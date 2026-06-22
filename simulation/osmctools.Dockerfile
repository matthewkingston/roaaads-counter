# Minimal image providing osmconvert + osmfilter (osmctools) for build_network.py.
# build_network.py uses them to stream a small bbox + drivable-highway extract out
# of the full NI .osm.pbf (the same snapshot OSRM is built from) at ~0.5 GB RAM.
#
# osmctools is used in preference to osmium-tool here because osmium sizes its
# referenced-node id-set by OSM's max node id (~12e9), needing ~2-3+ GB regardless
# of extract area — which OOMs modest machines. osmconvert/osmfilter size memory by
# the node count actually in the region instead.
#
# Build once (build_network.py auto-builds it on first run if absent):
#   docker build -f simulation/osmctools.Dockerfile -t osmctools-roaaads simulation/
FROM debian:stable-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends osmctools ca-certificates \
 && rm -rf /var/lib/apt/lists/*
