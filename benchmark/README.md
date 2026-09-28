# SULT1A3 / SULT1A1 benchmark data

Raw results behind the [Validation](../README.md#validation-sult1a3-vs-sult1a1) section.
Six ligands against two sulfotransferases that are 92.9% identical but differ in substrate
preference, scored by Nesso-1 and by AutoDock Vina.

## Files

| file | contents |
|---|---|
| `ligands.smi` | the six ligands, SMILES + name |
| `sult1a3.fasta` | SULT1A3 sequence, 295 aa (from the boltz_local example input) |
| `sult1a1.fasta` | SULT1A1 sequence, 295 aa (UniProt P50225) |
| `nesso_sult1a3.csv` | Nesso results, ranked strongest first |
| `nesso_sult1a1.csv` | as above for SULT1A1 |
| `affinity/<target>/<ligand>.json` | raw Nesso output, full precision |
| `pocket_<target>.txt` | residues Nesso scored: name, count, comma-separated numbers |
| `vina_scores.csv` | AutoDock Vina scores for both targets |

## How these were produced

Nesso, once per target:

```sh
zsh code/run_screen.sh benchmark/sult1a3.fasta benchmark/ligands.smi sult1a3
zsh code/run_screen.sh benchmark/sult1a1.fasta benchmark/ligands.smi sult1a1
```

Vina, via the `dock_assist` repo, docking into the crystallographic site of each structure
— 2A3R for SULT1A3 (L-dopamine bound) and 1LS6 for SULT1A1 (p-nitrophenol). Box centres
are the centroid of the catalytic co-crystal ligand, chosen as the copy within 5 A of the
PAP cofactor; 1LS6 holds a second p-nitrophenol 9.1 A away at a non-catalytic site, which
is excluded. Cubic box 22 A, exhaustiveness 16, seed 42.

| target | structure | ligand | box centre |
|---|---|---|---|
| SULT1A3 | 2A3R | L-dopamine, chain A | (50.536, 117.105, −1.369) |
| SULT1A1 | 1LS6 | p-nitrophenol, res 3001 | (17.350, 104.843, 57.613) |

## Notes on reading the data

* Nesso values are `log10(IC50 / uM)`, lower is stronger. Vina scores are kcal/mol, also
  lower is stronger. The two are not on a common scale; only rankings compare.
* Nesso inference is `bf16-mixed` and the outputs sit on the bf16 grid, which steps by
  ~0.008 near these values. Differences below ~0.01 are not meaningful, and exact ties
  between similar ligands are expected rather than suspicious.
* `n_pocket_residues` is populated for one ligand per target. The crop is a property of
  the protein and barely varies across these ligands — 115 residues on SULT1A3, 110 on
  SULT1A1 — so one per target is enough to show what was scored.
* The `predictions.safetensors` files holding the full pair representations are ~28 MB
  each and are not tracked here; they stay under `outputs/`.
