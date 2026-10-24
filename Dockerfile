# Environment for builds and checks; the repository is mounted, not copied, so the image serves every checkout:
#   docker build -t hsrdb .
#   docker run --rm -v "$PWD:/w" hsrdb python -m hsrdb_tools add cod/ --hs4 COD_HS4x.hsrdb --hs3 COD_HS3x.hsrdb
FROM python:3.12.14-slim-trixie AS base
RUN apt-get update && apt-get install -y --no-install-recommends sqlite3 p7zip-full rsync && rm -rf /var/lib/apt/lists/*
COPY requirements-release.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements-release.txt
WORKDIR /w

# cctbx for dev/verify_cctbx.py, in its own environment: docker build --target dev -t hsrdb-dev .
FROM base AS dev
COPY --from=mambaorg/micromamba:2.9.0 /bin/micromamba /usr/local/bin/micromamba
# the libraries at the versions of requirements-release.txt, so that the cross-check runs the release's numerics
RUN micromamba create -y -p /opt/cctbx -c conda-forge python=3.12 cctbx-base=2026.7 gemmi=0.7.5 numpy=2.5.3 spglib=2.7.0 \
 && micromamba clean -afy

FROM base
