#!/usr/bin/env bash
# Zero-config launcher for people who have Python 3 but nothing else installed.
# Creates a local virtual environment the first time, installs the package into
# it, then forwards all arguments to the `crispra-screen` command.
#
#   ./run.sh qc    --fastq my_library_R1.fastq.gz --outdir qc_out
#   ./run.sh score --treatment treat_R1.fastq.gz --use-bundled-mock --outdir score_out
#
# Re-runs reuse the environment, so only the first call pays the install cost.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
venv="$here/.venv"

if [ ! -x "$venv/bin/crispra-screen" ]; then
  echo "[run.sh] setting up a local environment in .venv (first run only) ..." >&2
  # --system-site-packages reuses numpy/pandas/scipy/matplotlib if the Python
  # you launched already has them (e.g. Anaconda), so only the small package
  # itself is fetched; a plain internet connection covers the rest.
  python3 -m venv --system-site-packages "$venv"
  "$venv/bin/pip" install --quiet --upgrade pip || true
  "$venv/bin/pip" install --quiet "$here"
fi

exec "$venv/bin/crispra-screen" "$@"
