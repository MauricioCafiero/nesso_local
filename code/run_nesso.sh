#!/bin/zsh
# Run one Nesso-1 prediction (protein sequence + ligand -> binding affinity).
# Results go to outputs/predictions/<input name>/affinity.json.
# usage: zsh code/run_nesso.sh code/tutorial_examples/smiles.yaml [extra nesso options]
root="$(cd "$(dirname "$0")/.." && pwd)"
caffeinate -w $$ &   # keep the Mac awake until the run finishes
# xet chunked transfer stalled on flaky wifi; plain HTTPS LFS resumes properly
export HF_HUB_DISABLE_XET=1
"$root/.venv/bin/nesso" predict "$@" --out_dir "$root/outputs"