#!/usr/bin/env bash
# Convert local frozen data/ and benchmark/ trees on a workstation or a Slurm node.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if [[ ! -f "$repo/rarecell/single_agent_data.py" ]]; then
    # sbatch runs a spooled copy of the script. Recover the submission checkout.
    repo=${SLURM_SUBMIT_DIR:?submit this script from the MultiGPTRareCell checkout}
fi
cd "$repo"
py=${RARECELL_PYTHON:-"$repo/.venv/bin/python"}
data=${DATA_ROOT:-"$repo/../data"}
bench=${BENCHMARK_ROOT:-"$repo/../LLMRareCell/benchmark"}
out=${DATASET_OUT:-"$repo/work/single-agent-datasets"}
if [[ $# == 0 ]]; then set -- all; fi
exec "$py" -m rarecell.single_agent_data --data-root "$data" --benchmark-root "$bench" \
    --all-views --out "$out" "$@"
