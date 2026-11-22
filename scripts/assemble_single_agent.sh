#!/usr/bin/env bash
# Import existing frozen releases; build any missing releases with the original preparer.
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
project=${SINGLE_AGENT_PROJECT:-"$(dirname "$repo")/LLMRareCell"}
out=${DATASET_OUT:-"$repo/work/single-agent-datasets"}
py=${RARECELL_PYTHON:-"$repo/.venv/bin/python"}
prep_py=${PREP_PYTHON:-"$project/.venv/bin/python"}
code=${SINGLE_AGENT_CODE:-"$repo/work/single-agent-preparer"}
pin=eed0c7a2d07f6772fbc85c0cffea2e56244c30e2
if [ "$#" -eq 0 ]; then
  set -- villani pgc cellsius pbmc_tcr baron cellsius_spheroid gse156490 airway pbmc68k
fi
if [ ! -x "$py" ]; then echo "Set RARECELL_PYTHON to the MultiGPTRareCell environment's Python." >&2; exit 1; fi
build() {
  local name=$1
  local args=()
  if [ ! -x "$prep_py" ]; then
    echo "Set PREP_PYTHON to an environment with single-agent module/preparation/requirements.txt installed." >&2
    return 1
  fi
  if [ ! -d "$code" ]; then
    mkdir -p "$(dirname "$code")" || return
    git clone --no-checkout https://github.com/Egria/LLMRareCell.git "$code" || return
    git -C "$code" checkout --detach "$pin" || return
  fi
  if [ "$(git -C "$code" rev-parse HEAD)" != "$pin" ] || [ -n "$(git -C "$code" status --porcelain --untracked-files=no -- module/preparation module/datasets.py)" ]; then
    echo "Use a clean single-agent preparer checkout at $pin via SINGLE_AGENT_CODE." >&2; return 1
  fi
  if [ "$name" = cellsius_spheroid ] && [ -n "${SPHEROID_EXPORT:-}" ]; then
    args+=(--spheroid-export "$SPHEROID_EXPORT")
  fi
  if [ "${INCLUDE_EXTRA_VIEWS:-0}" = 1 ]; then args+=(--include-extra-views); fi
  "$prep_py" "$code/module/preparation/prepare.py" prepare --projectdir "$project" --dataset "$name" "${args[@]}"
}
cd "$repo"
failed=0
for requested in "$@"; do
  case "$requested" in
    pgc) name=human_pgc ;;
    spheroid) name=cellsius_spheroid ;;
    villani|human_pgc|cellsius|pbmc_tcr|baron|cellsius_spheroid|gse156490|airway|pbmc68k) name=$requested ;;
    *) echo "Unknown dataset: $requested" >&2; failed=1; continue ;;
  esac
  if [ ! -f "$project/benchmark/$name/release.lock.json" ]; then
    if ! build "$name"; then echo "Preparation failed for $name" >&2; failed=1; continue; fi
  fi
  opts=()
  if [ "${INCLUDE_EXTRA_VIEWS:-0}" = 1 ]; then opts+=(--all-views); fi
  "$py" -m rarecell.single_agent_data --project "$project" --out "$out" "${opts[@]}" "$name" || failed=1
done
exit "$failed"
