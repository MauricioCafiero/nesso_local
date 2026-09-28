#!/usr/bin/env python
"""Screen many ligands against one protein with Nesso-1.

Builds one input YAML per ligand, runs them as a single `nesso predict` call on
the directory (one model load, and the protein's ESM embedding computed once for
the whole set), then collects the per-ligand affinity.json files into a table
ranked strongest first.

usage:
    python code/screen_nesso.py target.fasta ligands.smi --name mytarget
    python code/screen_nesso.py target.fasta ligands.smi --name mytarget --collect-only
    ... [extra flags passed through to nesso predict, e.g. --accelerator cpu]

The ligand file is one SMILES per line, with an optional name in the second
column; blank lines and lines starting with # are ignored. The protein is a
FASTA file or a plain file holding just the sequence.

Results land in outputs/screens/<name>/:
    inputs/         generated YAMLs, one per ligand
    predictions/    nesso output, one directory per ligand
    results.csv     ranked table

Re-running resumes: nesso skips any ligand that already has an affinity.json
(pass --override to force recomputation), so an interrupted screen picks up
where it stopped.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import yaml
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Screen many ligands against one protein with Nesso-1.",
        epilog="Unrecognised flags are passed through to `nesso predict`.",
    )
    parser.add_argument("protein", type=Path, help="FASTA file, or a file holding just the sequence")
    parser.add_argument("ligands", type=Path, help="SMILES file: 'SMILES [name]' per line")
    parser.add_argument("--name", help="screen name; defaults to the ligand file's stem")
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="where to put the screen; defaults to outputs/screens/<name>",
    )
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="skip prediction, just rebuild results.csv from whatever has finished",
    )
    args, passthrough = parser.parse_known_args()

    name = args.name or args.ligands.stem
    screen_dir = args.out_dir or ROOT / "outputs" / "screens" / name
    inputs_dir = screen_dir / "inputs"

    # Keep the Mac awake for the length of the screen; dies with this process.
    subprocess.Popen(["caffeinate", "-w", str(os.getpid())])

    if not args.collect_only:
        text = args.protein.read_text()
        if text.lstrip().startswith(">"):
            text = "".join(l for l in text.splitlines() if not l.startswith(">"))
        sequence = "".join(text.split()).upper()
        if not sequence:
            sys.exit(f"No sequence found in {args.protein}")

        # A bad SMILES anywhere in the list would take down the whole batch, and
        # vendor/ChEMBL exports routinely carry a few, so drop them up front.
        RDLogger.DisableLog("rdApp.*")
        ligands, seen = [], set()
        for n, line in enumerate(args.ligands.read_text().splitlines(), 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            smiles = fields[0]
            if Chem.MolFromSmiles(smiles) is None:
                print(f"  skipping line {n}: unparseable SMILES {smiles!r}")
                continue
            # The record id is the YAML stem, so it must be unique and path-safe.
            label = fields[1] if len(fields) > 1 else f"lig{n:04d}"
            label = re.sub(r"[^A-Za-z0-9._-]", "_", label)
            if label in seen:
                label = f"{label}_{n}"
            seen.add(label)
            ligands.append((label, smiles))

        if not ligands:
            sys.exit(f"No usable ligands in {args.ligands}")

        inputs_dir.mkdir(parents=True, exist_ok=True)
        for label, smiles in ligands:
            spec = {
                "sequences": [
                    {"protein": {"id": "A", "sequence": sequence}},
                    {"ligand": {"id": "B", "smiles": smiles}},
                ],
                "properties": [{"affinity": {"binder": "B"}}],
            }
            (inputs_dir / f"{label}.yaml").write_text(yaml.safe_dump(spec, sort_keys=False))

        print(f"{len(ligands)} ligands -> {inputs_dir}")

        env = dict(os.environ)
        env["HF_HUB_DISABLE_XET"] = "1"  # xet transfer stalls on weak wifi; LFS resumes
        result = subprocess.run(
            [str(ROOT / ".venv" / "bin" / "nesso"), "predict", str(inputs_dir),
             "--out_dir", str(screen_dir), "--save_metadata", *passthrough],
            env=env,
        )
        if result.returncode != 0:
            print("nesso predict failed; collecting whatever finished", file=sys.stderr)

    # The pocket the model actually scored. A prediction against the wrong
    # residues is meaningless however confident it looks, so record it.
    prot_text = args.protein.read_text()
    if prot_text.lstrip().startswith(">"):
        prot_text = "".join(l for l in prot_text.splitlines() if not l.startswith(">"))
    n_prot = len("".join(prot_text.split()))

    pockets = {}
    for meta in sorted((screen_dir / "predictions").glob("*/predictions.safetensors")):
        try:
            from safetensors.torch import safe_open
            with safe_open(meta, "pt") as fh:
                mask = fh.get_tensor("pocket_mask").tolist()
        except Exception as exc:
            print(f"  could not read pocket for {meta.parent.name}: {exc}")
            continue
        res = [i + 1 for i, v in enumerate(mask) if v and (not n_prot or i < n_prot)]
        pockets[meta.parent.name] = res

    if pockets:
        with (screen_dir / "pockets.txt").open("w") as fh:
            for name, res in sorted(pockets.items()):
                fh.write(f"{name}\t{len(res)}\t{','.join(map(str, res))}\n")

    rows = []
    for path in sorted((screen_dir / "predictions").glob("*/affinity.json")):
        d = json.loads(path.read_text())
        value = d["affinity_pred_value"]
        rows.append({
            "ligand": path.parent.name,
            "affinity_pred_value": round(value, 4),
            "pIC50": round(6 - value, 2),
            "ensemble_spread": round(abs(d["affinity_pred_value1"] - d["affinity_pred_value2"]), 4),
            "binder_probability": round(d["affinity_probability_binary"], 4),
            "entropy_crop_pl": round(d["entropy_crop_pl"], 4),
            "n_pocket_residues": len(pockets.get(path.parent.name, [])) or "",
        })

    if not rows:
        sys.exit(f"No predictions found under {screen_dir / 'predictions'}")

    rows.sort(key=lambda r: r["affinity_pred_value"])  # lower is stronger
    results = screen_dir / "results.csv"
    with results.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{len(rows)} results -> {results}")
    print(f"{'ligand':<24} {'value':>8} {'pIC50':>7} {'entropy_pl':>11}")
    for r in rows[:10]:
        print(f"{r['ligand']:<24} {r['affinity_pred_value']:>8} {r['pIC50']:>7} {r['entropy_crop_pl']:>11}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
