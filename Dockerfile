FROM python:3.12-slim AS base
RUN apt-get update && apt-get install -y --no-install-recommends sqlite3 p7zip-full rsync && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir gemmi numpy spglib
WORKDIR /w

# cctbx for dev/verify_cctbx.py, in its own environment: docker build --target dev -t hsrdb-dev .
FROM base AS dev
COPY --from=mambaorg/micromamba:2.9.0 /bin/micromamba /usr/local/bin/micromamba
RUN micromamba create -y -p /opt/cctbx -c conda-forge python=3.12 cctbx-base=2026.7 gemmi numpy spglib \
 && micromamba clean -afy

FROM base
