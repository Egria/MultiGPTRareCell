#!/bin/bash
#SBATCH --job-name=rarecell
#SBATCH --export=NIL
set -e
cd "${SLURM_SUBMIT_DIR:?Submit from the MultiGPTRareCell repository root}"
if [[ ! -x .venv/bin/python || ! -f rarecell/r/annotation.R ]]; then
    echo "Submit from the MultiGPTRareCell checkout containing .venv/bin/python" >&2
    exit 2
fi
exec .venv/bin/python -m rarecell "$@"
