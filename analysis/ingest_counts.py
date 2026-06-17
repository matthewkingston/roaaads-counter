"""
Ingest raw walking traffic count CSVs into a processed sessions file.

Reads all CSVs from data/counts/, extracts one record per session, and
writes data/counts_processed.json. Idempotent: existing session_ids are
skipped (and already-snapped sessions are skipped) so new data files can
be added and this script re-run safely.

Each session is matched to exactly one road link in the simulation network
using a GPS-accuracy-weighted least-squares fit, then the observer's
direction of travel is used to assign directed "with" and "against" links.

Usage:
  python3 analysis/ingest_counts.py
"""

import csv, glob, json, math, os
from collections import defaultdict
from datetime import datetime

import osmnx as ox
import pytz
from pyproj import Transformer
from shapely.geometry import Point

COUNTS_DIR       = "data/counts"
PROCESSED_FILE   = "data/counts_processed.json"
CONS_GRAPH       = "simulation/newtownards_consolidated.graphml"
HOURLY_FRACTIONS = "analysis/hourly_fractions.csv"

LOCAL_TZ = pytz.timezone("Europe/London")

# Sessions where GPS snap would land on the wrong link (e.g. observer on the
# parallel carriageway of a dual one-way road).  Maps session_id → forced
# directed links.  Takes priority over auto-snap; idempotent on re-ingest.
MANUAL_LINK_OVERRIDES = {
    # A20 Kempe Stones eastbound: observer was on the westbound carriageway
    # and could not access the eastbound side.
    "e644eae2": {"link_with": [8, 7], "link_against": [7, 8]},
    "760b0c8e": {"link_with": [8, 7], "link_against": [7, 8]},
}


def parse_iso(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def _u(n):
    """Jeffreys-prior Poisson uncertainty: sqrt(n + 0.5); None-safe."""
    return round(math.sqrt(n + 0.5), 3) if n is not None else None


# ── Load hourly fraction profile ─────────────────────────────────────────────

hourly_fracs = {}  # (day_of_week, hour) → (mean_fraction, std_fraction)
with open(HOURLY_FRACTIONS, newline="") as f:
    for row in csv.DictReader(f):
        dow  = int(row["day_of_week"])
        hour = int(row["hour"].split(":")[0])
        hourly_fracs[(dow, hour)] = (float(row["mean_fraction"]), float(row["std_fraction"]))


# ── Load all raw CSVs ────────────────────────────────────────────────────────

# sessions[sid]["gps"] = [(lat, lng, accuracy_m), ...]
sessions = defaultdict(lambda: {"meta": None, "gps": [], "cars": [], "files": set()})

csv_files = sorted(glob.glob(f"{COUNTS_DIR}/*.csv"))
if not csv_files:
    print(f"No CSV files found in {COUNTS_DIR}/")
    raise SystemExit(1)

seen_rows = set()  # (session_id, event_type, timestamp) — deduplicates sessions across files

for path in csv_files:
    filename = os.path.basename(path)
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            sid = row["session_id"]
            sessions[sid]["files"].add(filename)
            if sessions[sid]["meta"] is None:
                sessions[sid]["meta"] = {
                    "label":     row["session_label"],
                    "mode":      row["session_mode"],
                    "start_utc": row["session_start"],
                    "end_utc":   row["session_end"],
                }
            row_key = (sid, row["event_type"], row["timestamp"])
            if row_key in seen_rows:
                continue
            seen_rows.add(row_key)
            if row["event_type"] == "gps_track":
                sessions[sid]["gps"].append((
                    float(row["lat"]),
                    float(row["lng"]),
                    float(row["gps_accuracy_m"]),
                ))
            elif row["event_type"] == "car":
                sessions[sid]["cars"].append(row["direction"])

print(f"Read {len(csv_files)} file(s), found {len(sessions)} unique session(s)")

# ── Load existing processed file ─────────────────────────────────────────────

if os.path.exists(PROCESSED_FILE):
    with open(PROCESSED_FILE) as f:
        processed = json.load(f)
else:
    processed = {"sessions": {}}

# ── Load graph and build undirected edge geometry list ───────────────────────

needs_snap = [
    sid for sid, data in sessions.items()
    if sid not in processed["sessions"] or "matched_link_with" not in processed["sessions"][sid]
]

if needs_snap:
    print(f"Loading graph for link snapping ({len(needs_snap)} session(s) to snap) …")
    G = ox.load_graphml(CONS_GRAPH)

    to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)

    # Deduplicated undirected edges: store actual (u, v) so geometry direction
    # matches best_u/best_v in snap_to_link (fixes direction flip on one-way
    # roads where u > v).
    seen_pairs = set()
    edge_geoms = []
    for u, v, edata in G.edges(data=True):
        if "geometry" not in edata:
            continue  # virtual edges (e.g. stub nodes) have no geometry
        pair = (min(u, v), max(u, v))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        edge_geoms.append((u, v, edata["geometry"]))

    print(f"  {len(edge_geoms)} unique undirected links")

    def snap_to_link(gps_pts):
        """
        gps_pts: list of (lat, lng, accuracy_m)
        Returns (link_with, link_against, rmse_m) where link_with/against are [u, v].
        """
        utm_pts = []
        weights = []
        for lat, lng, acc in gps_pts:
            e, n = to_utm.transform(lng, lat)
            utm_pts.append(Point(e, n))
            weights.append(1.0 / (acc * acc))

        total_w = sum(weights)

        # Find best undirected edge by weighted SSE
        best_u, best_v, best_geom = None, None, None
        best_wsse = float("inf")
        for u, v, geom in edge_geoms:
            wsse = sum(w * pt.distance(geom) ** 2 for pt, w in zip(utm_pts, weights))
            if wsse < best_wsse:
                best_wsse = wsse
                best_u, best_v, best_geom = u, v, geom

        rmse_m = math.sqrt(best_wsse / total_w) if total_w > 0 else None

        # Determine observer direction: vector from first to last UTM point
        obs_dx = utm_pts[-1].x - utm_pts[0].x
        obs_dy = utm_pts[-1].y - utm_pts[0].y

        # Edge direction: from linestring start to end
        coords = list(best_geom.coords)
        edge_dx = coords[-1][0] - coords[0][0]
        edge_dy = coords[-1][1] - coords[0][1]

        dot = obs_dx * edge_dx + obs_dy * edge_dy
        if dot >= 0:
            link_with    = [best_u, best_v]
            link_against = [best_v, best_u]
        else:
            link_with    = [best_v, best_u]
            link_against = [best_u, best_v]

        return link_with, link_against, round(rmse_m, 1) if rmse_m is not None else None

# ── Process new sessions and snap ────────────────────────────────────────────

new_count         = 0
snap_count        = 0
uncertainty_count = 0
aadt_count        = 0

for sid, data in sessions.items():
    is_new = sid not in processed["sessions"]

    if is_new:
        meta = data["meta"]
        start = parse_iso(meta["start_utc"])
        end   = parse_iso(meta["end_utc"])
        duration_s = (end - start).total_seconds()

        with_count    = sum(1 for d in data["cars"] if d == "with")
        against_count = sum(1 for d in data["cars"] if d == "against")
        total_count   = with_count + against_count

        if meta["mode"] == "single":
            # Determine which direction was being recorded from car events.
            # If no cars observed we can't tell → discard all counts (null).
            dirs = {d for d in data["cars"]}
            if not dirs:
                with_count = against_count = total_count = None
            elif dirs == {"with"}:
                against_count = total_count = None
            elif dirs == {"against"}:
                with_count = total_count = None
            else:
                # Both directions present in single mode → ambiguous, discard
                with_count = against_count = total_count = None

        gps_pts = data["gps"]
        if gps_pts:
            centroid_lat = sum(p[0] for p in gps_pts) / len(gps_pts)
            centroid_lng = sum(p[1] for p in gps_pts) / len(gps_pts)
        else:
            centroid_lat = centroid_lng = None

        processed["sessions"][sid] = {
            "session_id":    sid,
            "label":         meta["label"],
            "mode":          meta["mode"],
            "start_utc":     meta["start_utc"],
            "end_utc":       meta["end_utc"],
            "duration_s":    round(duration_s, 1),
            "with_count":    with_count,
            "against_count": against_count,
            "total_count":   total_count,
            "centroid_lat":  round(centroid_lat, 6) if centroid_lat is not None else None,
            "centroid_lng":  round(centroid_lng, 6) if centroid_lng is not None else None,
            "source_files":  sorted(data["files"]),
        }
        new_count += 1

    # Jeffreys-prior Poisson uncertainty: sqrt(N + 0.5); null if count is null.
    # Gives non-zero uncertainty for N=0; converges to sqrt(N) for large N.
    if processed["sessions"][sid].get("uncertainty_method") != "jeffreys":
        rec = processed["sessions"][sid]
        rec["with_uncertainty"]    = _u(rec["with_count"])
        rec["against_uncertainty"] = _u(rec["against_count"])
        rec["total_uncertainty"]   = _u(rec["total_count"])
        rec["uncertainty_method"]  = "jeffreys"
        uncertainty_count += 1

    # Snap if new or previously processed without snapping
    if "matched_link_with" not in processed["sessions"][sid]:
        if sid in MANUAL_LINK_OVERRIDES:
            ov = MANUAL_LINK_OVERRIDES[sid]
            processed["sessions"][sid]["matched_link_with"]    = ov["link_with"]
            processed["sessions"][sid]["matched_link_against"] = ov["link_against"]
            processed["sessions"][sid]["match_rmse_m"]         = None
            processed["sessions"][sid]["match_method"]         = "manual"
            print(f"  {sid}  MANUAL link {ov['link_with'][0]}→{ov['link_with'][1]} "
                  f"(against: {ov['link_against'][0]}→{ov['link_against'][1]})")
        else:
            gps_pts = data["gps"]
            if gps_pts:
                link_with, link_against, rmse_m = snap_to_link(gps_pts)
                processed["sessions"][sid]["matched_link_with"]    = link_with
                processed["sessions"][sid]["matched_link_against"] = link_against
                processed["sessions"][sid]["match_rmse_m"]         = rmse_m
                print(f"  {sid}  link {link_with[0]}→{link_with[1]} "
                      f"(against: {link_against[0]}→{link_against[1]})  "
                      f"RMSE={rmse_m}m")
        # Validate that every non-null count maps to a real directed edge.
        # A count on a non-existent edge silently contributes zero to the model.
        rec = processed["sessions"][sid]
        lw, la = rec.get("matched_link_with"), rec.get("matched_link_against")
        if lw and la:
            if rec.get("with_count") is not None and not G.has_edge(lw[0], lw[1]):
                raise ValueError(
                    f"\nSession {sid}: with_count={rec['with_count']} assigned to "
                    f"{lw[0]}→{lw[1]} but that directed edge does not exist in the network. "
                    f"Add an entry to MANUAL_LINK_OVERRIDES with the correct link."
                )
            if rec.get("against_count") is not None and not G.has_edge(la[0], la[1]):
                raise ValueError(
                    f"\nSession {sid}: against_count={rec['against_count']} assigned to "
                    f"{la[0]}→{la[1]} but that directed edge does not exist in the network. "
                    f"Add an entry to MANUAL_LINK_OVERRIDES with the correct link."
                )
        snap_count += 1

    # AADT estimation via hourly fraction profile
    if processed["sessions"][sid].get("aadt_method") != "hourly_fraction_v3":
        rec  = processed["sessions"][sid]
        T    = rec["duration_s"]
        dt_l = parse_iso(rec["start_utc"]).astimezone(LOCAL_TZ)
        mean_f, std_f = hourly_fracs[(dt_l.weekday(), dt_l.hour)]

        def _aadt(n, sigma_n):
            if n is None:
                return None, None
            k      = 3600.0 / (T * mean_f)
            # Use Jeffreys posterior mean (n + 0.5) for both the point estimate
            # and the fraction-uncertainty term, so the estimate is consistent
            # with the Jeffreys uncertainty (sigma_n = sqrt(n + 0.5)).
            # For n >> 1 the 0.5 correction is negligible; for n=0 it shifts
            # the point estimate from 0 to ~0.5/T rather than anchoring at 0.
            n_eff  = n + 0.5
            aadt   = round(n_eff * k)
            sigma  = round(k * math.sqrt(sigma_n ** 2 + (n_eff * std_f / mean_f) ** 2))
            return aadt, sigma

        rec["with_aadt"],    rec["with_aadt_uncertainty"]    = _aadt(rec["with_count"],    rec["with_uncertainty"])
        rec["against_aadt"], rec["against_aadt_uncertainty"] = _aadt(rec["against_count"], rec["against_uncertainty"])
        rec["total_aadt"],   rec["total_aadt_uncertainty"]   = _aadt(rec["total_count"],   rec["total_uncertainty"])
        rec["hourly_fraction"]     = round(mean_f, 6)
        rec["hourly_fraction_std"] = round(std_f,  6)
        rec["time_slot"]           = [dt_l.weekday(), dt_l.hour]
        rec["frac_rel_std"]        = round(std_f / mean_f, 6)
        rec["aadt_method"]         = "hourly_fraction_v3"
        aadt_count += 1

# ── Save ─────────────────────────────────────────────────────────────────────

with open(PROCESSED_FILE, "w") as f:
    json.dump(processed, f, indent=2)

already_present = len(processed["sessions"]) - new_count
total = len(processed["sessions"])
print(f"{new_count} new session(s) added, {already_present} already present, "
      f"{snap_count} snapped, {uncertainty_count} uncertainty, {aadt_count} AADT estimated, {total} total")
print(f"Saved → {PROCESSED_FILE}")
