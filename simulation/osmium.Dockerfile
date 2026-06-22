# Minimal image providing the `osmium` CLI (osmium-tool) for build_network.py.
# build_network.py uses it to stream a small bbox + drivable-highway extract out
# of the full NI .osm.pbf (the same snapshot OSRM is built from) without loading
# the whole 400 MB file into memory.
#
# Build once:
#   docker build -f simulation/osmium.Dockerfile -t osmium-roaaads simulation/
FROM debian:stable-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends osmium-tool ca-certificates \
 && rm -rf /var/lib/apt/lists/*
ENTRYPOINT ["osmium"]
