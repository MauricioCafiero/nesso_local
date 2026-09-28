# nesso_local

A local setup for running [Nesso-1](https://github.com/recursionpharma/nesso) (Recursion
Pharmaceuticals) on Apple silicon: helper scripts, the upstream tutorial inputs, and the
upstream CPU-only test suite, arranged so a prediction is a single command.

Nesso-1 is a coarse-grained cofolding model that predicts binding affinity. It takes a
protein sequence plus a ligand (SMILES, CCD code, or SDF) and returns an affinity scalar.
It folds the pair internally to get there, so no experimental structure, docking grid, or
pocket definition is needed — sequence in, number out. It is Apache-2.0 licensed; the
method and its evaluations are described in the
[technical report](https://www.biorxiv.org/content/10.64898/2026.08.01.742196v1).

This repository contains no model code; `nesso` itself is installed as a dependency.

## What it's good for

The model is small (165 MB of weights) and there is no structure-preparation step, which
makes it practical for volume on modest hardware: the tutorial example takes roughly three
minutes on an M-series laptop once the cache is warm. That suits

* ranking a set of candidate ligands against one target, especially a target with no
  solved structure,
* comparing analogs within a series to see which direction is worth pursuing,
* a cheap first-pass filter ahead of docking, free-energy work, or assays.

What it will not do:

* **Give you a pose.** A run writes `affinity.json` and nothing else — no coordinates. It
  is an affinity predictor, not a structure predictor, and the internal fold is not
  exported.
* **Aim at a specific pocket.** Pocket conditioning and structural templating are listed
  upstream as not yet implemented, so you cannot point it at one site of a multi-site
  protein.
* **Replace a measurement.** Treat the output as a ranking signal, and check
  `entropy_crop_pl` before trusting any single number (see below).

## Install

Python 3.10–3.13 (upstream constraint is `>=3.10,<3.14`). Tested on 3.12.

```sh
python -m venv .venv
.venv/bin/pip install "nesso @ git+https://github.com/recursionpharma/nesso.git"
```

Nesso is not on PyPI yet, so it installs from git; that pulls in its pinned runtime
dependencies. `requirements.txt` records the resolved versions for reference.

The first prediction downloads about 3.1 GB of weights into `.cache/`: the CCD dictionary
(`ccd.pkl`, 413 MB), the model checkpoint (165 MB), and ESM-2
(`facebook/esm2_t33_650M_UR50D`, ~2.5 GB). That directory is gitignored but should not be
deleted — everything after the first run reuses it.

## Run

```sh
zsh code/run_nesso.sh code/tutorial_examples/smiles.yaml
```

Extra flags pass straight through to the CLI, e.g. `--save_metadata`, `--override`,
`--accelerator cpu`. The wrapper pins `--out_dir outputs`, holds off idle sleep with
`caffeinate` for the length of the run, and sets `HF_HUB_DISABLE_XET=1` so the weight
downloads use plain HTTPS LFS, which resumes cleanly on a weak connection.

Results land in `outputs/predictions/<record_id>/affinity.json`, where `record_id` is the
stem of the input YAML. Existing predictions are reused unless you pass `--override`, so
give each protein–ligand pair its own filename.

## Screening many ligands against one protein

```sh
.venv/bin/python code/screen_nesso.py target.fasta ligands.smi --name mytarget
```

The ligand file is one SMILES per line with an optional name in the second column; `#`
comments and blank lines are ignored, unparseable SMILES are reported and skipped, and
repeated names are made unique. The script writes one YAML per ligand and hands the whole
directory to `nesso predict` in a single call, so the model loads once and the protein's
ESM embedding is computed once for the entire set rather than per ligand. Extra flags
(`--accelerator cpu`, `--override`, …) pass straight through.

Everything lands in `outputs/screens/<name>/`: the generated `inputs/`, the per-ligand
`predictions/`, and `results.csv` ranked strongest first, with `affinity_pred_value`,
pIC50, the ensemble spread, the binder probability, and `entropy_crop_pl` for each ligand.
The top ten are printed at the end.

A screen resumes cleanly. Ligands that already have an `affinity.json` are skipped on a
re-run, so an interrupted screen continues where it stopped, and `--collect-only` rebuilds
`results.csv` from whatever has finished without running any prediction. For a long screen,
launch it detached so it survives the terminal:

```sh
.venv/bin/python -c "import subprocess; subprocess.Popen(['.venv/bin/python','code/screen_nesso.py','target.fasta','ligands.smi','--name','mytarget'], stdout=open('outputs/screen.log','ab'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)"
```

### Target size is what sets the runtime

Cost is dominated by the protein, not the ligand: pair representations scale roughly with
the square of the token count, with larger temporaries inside triangle attention, and the
cuEquivariance kernels that would cut this are CUDA-only. Measured on an 8 GB M-series
laptop: ~20-residue peptides run about **2 s per ligand**, while a 384-residue protein takes
**~170 s for a single ligand** and can exhaust RAM and swap outright — the process pages out
and wedges rather than merely running slowly.

So a hundred-ligand peptide screen is a few minutes, and the same screen against a real
protein is hours, if it fits at all. On a memory-constrained machine, watch
`sysctl vm.swapusage` on a first run, and if a large target is the goal, trade accuracy for
headroom with `--accelerator cpu`, `--recycling_steps 1` (the default of 5 means six trunk
passes), `--no_refine_protein_inference`, or a smaller `--refine_protein_tokens_budget`.
Calibrate any cheapened setting against a few known ligands before trusting a whole screen.
A stalled run loses nothing — kill it and re-run to resume.

Note that only one ligand per YAML is the `binder`. `code/tutorial_examples/multi_ligand.yaml`
shows two ligands in one input, but that is a cofactor setup — a single affinity is
predicted for the designated binder, not one per ligand.

## Reading the output

* `affinity_pred_value` — log10(IC50 / µM), so **lower is stronger**: −3 ≈ 1 nM, 0 ≈ 1 µM,
  +2 ≈ 100 µM. Convert with pIC50 = 6 − value.
* `affinity_pred_value1` / `affinity_pred_value2` — the two ensemble members. Their spread
  is a cheap uncertainty estimate.
* `affinity_probability_binary` — a separate binder/non-binder classifier head.
* `entropy_crop_pl` — structural confidence. At 0.0 the model failed to place the ligand
  and the affinity should not be trusted.

## Layout

```
code/run_nesso.sh          prediction wrapper
code/tutorial_examples/    upstream tutorial YAMLs + extract_features.py
tutorial -> code/tutorial_examples   symlink the upstream tests expect
tests/                     upstream CPU-only test suite (no weights, no network)
outputs/                   prediction results (gitignored)
.cache/                    Hugging Face model cache (gitignored)
```

The tutorial inputs and tests were copied from upstream commit `6c72f66` on 2026-09-28;
see `code/tutorial_examples/README.md` for what each example does.

## Tests

```sh
.venv/bin/pytest tests
```

The suite runs on CPU and needs no weights. Seven CCD-dependent tests skip until a first
`nesso predict` has cached `ccd.pkl`; after that all 28 pass.

## License

The scripts here are MIT — see `LICENSE`. Nesso-1 itself is Apache-2.0; see the
[upstream repository](https://github.com/recursionpharma/nesso) for its license and
third-party notices.
