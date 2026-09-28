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

* separating scaffolds against one target, especially a target with no solved structure,
* a cheap first-pass filter ahead of docking, free-energy work, or assays.

Benchmarking on a sulfotransferase pair (see [Validation](#validation-sult1a3-vs-sult1a1))
found it reliable for that coarse separation and unreliable for finer distinctions, so
ranking close analogues is deliberately not on this list.

What it will not do:

* **Give you a pose.** A run writes `affinity.json` and nothing else — no coordinates. It
  is an affinity predictor, not a structure predictor, and the internal fold is not
  exported.
* **Aim at a specific pocket.** Pocket conditioning and structural templating are listed
  upstream as not yet implemented, so you cannot point it at one site of a multi-site
  protein.
* **Rank close analogues, or resolve isoform selectivity.** Measured on SULT1A3/SULT1A1:
  it inverted the published ranking of three cresol isomers and returned near-identical
  values for two enzymes that differ at the pocket, including the substitutions known to
  switch their substrate preference. See [Validation](#validation-sult1a3-vs-sult1a1).
* **Replace a measurement.** Treat the output as a ranking signal. Note that
  `entropy_crop_pl` read 0.44-0.61 throughout that benchmark, i.e. it looked confident
  while the rankings were wrong -- it reports placement, not accuracy.

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
protein is hours, if it fits at all. Measured on an 8 GB machine, the 384-residue tutorial
target needs about **5.7 GB** and stalls the process in uninterruptible wait rather than
failing cleanly. A 295-residue target needs about **4.0 GB** and finishes in ~3 minutes, so
on 8 GB the practical ceiling is roughly 300 residues, with no margin to spare.

Peak memory is fixed by the full sequence length and cannot be tuned down after the fact:
the pair representation is built at full length and the pocket crop only happens after the
first Pairformer pass, so `--refine_protein_tokens_budget`, `--affinity_protein_cutoff`
and `--recycling_steps` reduce *time*, not peak memory. What actually helps:

* **A shorter sequence.** Pair memory grows with the square of the token count, so passing
  only the domain or region around the binding site is the one large reduction available.
* **`--accelerator cpu`.** Same footprint, but CPU tensors are pageable and compressible
  where MPS allocations are wired and neither, so a tight machine swaps and finishes
  slowly instead of wedging.
* **`PYTORCH_MPS_HIGH_WATERMARK_RATIO`** to cap the MPS allocator, turning a silent stall
  into a clean out-of-memory error.

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

## Validation: SULT1A3 vs SULT1A1

Nesso was benchmarked on a system with independent reference data: the sulfotransferase
pair SULT1A3 and SULT1A1. The two enzymes are **92.9% identical** (21 substitutions over
295 residues) but differ sharply in substrate preference. SULT1A3 carries a charged back
pocket — Glu146 pairs with the protonated amine of catecholamines — while SULT1A1 has no
charged residue there apart from the catalytic lysine and so favours hydrophobic phenols.
Three of the 21 substitutions (D86A, E89I, E146A) are exactly the triple mutation that
experimentally converts SULT1A3's dopamine kinetics to SULT1A1-like, and E146A alone raises
the dopamine Km eightfold (~0.9 log units).

That makes the pair a controlled test: same ligands, near-identical sequences, one
remodelled pocket.

Six ligands were screened against both isoforms — dopamine, L-DOPA, paracetamol, and o-,
m- and p-cresol. The same ligands were docked into both with AutoDock Vina, at the
crystallographic site in each case (2A3R, with L-dopamine bound, and 1LS6, with
p-nitrophenol), using the tooling in `dock_assist`.

### Nesso does not distinguish the isoforms

| ligand | SULT1A3 | SULT1A1 | delta |
|---|---|---|---|
| dopamine | 1.266 | 1.219 | −0.047 |
| L-DOPA | 1.266 | 1.438 | +0.172 |
| o-cresol | 1.797 | 1.812 | +0.016 |
| m-cresol | 1.859 | 1.812 | −0.047 |
| paracetamol | 2.281 | 2.250 | −0.031 |
| p-cresol | 2.406 | 2.391 | −0.016 |

Five of six deltas fall within ±0.05. Inference runs in `bf16-mixed` and every output sits
on the bf16 grid, which steps by ~0.008 in this range, so those shifts are 2–6 steps — the
model's own output resolution. Where experiment gives ~0.9 log units weaker for dopamine,
nesso gives 0.047 in the wrong direction.

This is not a case of the model scoring the wrong region. Running with `--save_metadata`
and decoding `pocket_mask` shows the pocket contains D86, K106, H108 and E146 along with
all 11 dopamine atoms, and 15 of the 21 substitutions fall inside it. The model has the
relevant residues in view and is unmoved by changing them. The likely reason is dilution:
the crop keeps 115 of 295 residues, so a two-residue change in charge character is a small
perturbation to what the affinity head sees. That crop cannot be tightened from the CLI,
because it happens after the first full-length Pairformer pass.

### Docking recovers both the ranking and the selectivity

| ligand | Vina 1A3 | Vina 1A1 | delta (kcal/mol) |
|---|---|---|---|
| paracetamol | −6.2 | −3.1 | +3.1 |
| dopamine | −5.6 | −3.3 | +2.3 |
| L-DOPA | −5.6 | −1.7 | +3.9 |
| p-cresol | −5.4 | −4.0 | +1.4 |
| m-cresol | −5.3 | −4.1 | +1.2 |
| o-cresol | −5.2 | −4.0 | +1.2 |

On SULT1A3 this reproduces published MP2//DFT interaction energies for the same active
site (acetaminophen −47.93 kcal/mol strongest, p-cresol −41.25, o-cresol −34.69 weakest):
Vina puts paracetamol first and o-cresol last. Nesso inverts both, ranking o-cresol second
and paracetamol fifth.

Across the isoform swap the catecholamines lose 2–4 kcal/mol while the cresols lose only
~1.2, and the ranking inverts on SULT1A1 — the cresols become the best binders, ahead of
dopamine and paracetamol. That is the expected consequence of replacing a charged back
pocket with a hydrophobic one.

### Caveats

* The two dockings use different crystal structures, so the uniform part of the 1A3→1A1
  shift may be a systematic offset between them rather than selectivity. The robust result
  is the ranking inversion *within* SULT1A1, which involves one structure and one box.
* L-DOPA's −1.7 on SULT1A1 is a large outlier for the biggest, zwitterionic ligand and may
  be a docking artefact.
* Km is a substrate kinetic parameter while nesso predicts an inhibition-style affinity, so
  the observables are not identical — though no reasonable mismatch turns +0.9 into −0.05.
* Vina scores fragments this small within a ~1 kcal/mol spread, and neither method resolves
  m- from p-cresol. Both are coarse.

### What this means in practice

The two methods get reached for in the same situation — rank some ligands against a target,
quickly, without setting up anything expensive. On cost they are comparable, and nesso is
the more convenient of the two since it needs only a sequence. On this system the
physics-based method is the more accurate one: it recovers the reference ranking and the
isoform selectivity, and nesso does neither.

Nesso separated catechols from monophenols correctly on SULT1A3, which is the coarse
discrimination it should get right. Treat its output as scaffold-level triage, and do not
use it to rank close analogues or to reason about isoform selectivity without checking
against a structure-based method.

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
