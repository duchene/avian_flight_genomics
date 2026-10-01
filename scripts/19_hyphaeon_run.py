#!/usr/bin/env python3
"""19_hyphaeon_run.py

Run the per-gene HyphAeon phenotype screen across all retained genes and all three
phenotypes.

================================================================================
THIS SCRIPT TARGETS ONE SPECIFIC MACHINE: an Intel (x86_64) Mac.
================================================================================

PyTorch's last macOS x86_64 build is 2.2.2, which predates NumPy 2. HyphAeon's default
dependency resolution installs NumPy 2, and the run then dies at tensor.numpy() with
"Numpy is not available". Downgrading NumPy afterwards segfaults, because scipy's wheel
is built against NumPy 2. The entire stack must therefore be pinned to the NumPy 1 era
BEFORE hyphaeon is installed:

    python3 -m venv hvenv
    hvenv/bin/pip install --upgrade pip
    hvenv/bin/pip install "numpy==1.26.4" "scipy==1.11.4" "pandas==2.1.4" \
                          "matplotlib==3.8.4" "torch==2.2.2"
    hvenv/bin/pip install --no-deps hyphaeon hyphaeon-core
    hvenv/bin/pip install biopython tn93 safetensors huggingface-hub tqdm pyyaml \
                          click truststore

Run this script with that venv's interpreter, not the system Python.

On Apple Silicon or Linux none of this pinning is needed, because torch is current there.
If the screen is moved to such a machine, delete the version check below and install
hyphaeon normally. Nothing else in this script is platform-specific.

See HYPHAEON.md for the full notes.

--------------------------------------------------------------------------------
WHY A PYTHON DRIVER RATHER THAN THE CLI

Measured on a real gene: 1.5 s of analysis but 5.3 s per process, because importing torch
and loading the model dominate. Calling the CLI once per gene for 14,957 genes across
3 phenotypes would spend roughly 65 hours on process overhead alone. This driver imports
once and loops.

In:  ../output/hyphaeon_input/manifest.csv, ../output/hyphaeon_input/phenotypes.csv
Out: ../output/hyphaeon_out/<phenotype>/<gene>.json
     ../output/hyphaeon_out/progress_<phenotype>.log
"""

import argparse
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# The exceptiongroup backport in this dependency tree segfaults inside its own exception
# hook while formatting errors, which replaces every real message with a traceback from
# exceptiongroup/_formatting.py and then crashes. Take the hook back before anything can
# raise, or debugging is impossible.
try:
    import exceptiongroup  # noqa: F401
except ImportError:
    pass
sys.excepthook = sys.__excepthook__

HERE = os.path.dirname(os.path.abspath(__file__))
IN_DIR = os.path.join(HERE, "..", "output", "hyphaeon_input")
OUT_ROOT = os.path.join(HERE, "..", "output", "hyphaeon_out")

PHENOTYPES = {
    # name           column          continuous?
    "loading": ("loading_index", True),
    "flightloss": ("flight_loss", False),
    "migration": ("migratory", False),
}


def check_environment():
    """Fail loudly and usefully rather than at tensor.numpy() a minute in."""
    import platform

    import numpy as np
    import torch

    numpy_major = int(np.__version__.split(".")[0])
    msg = []
    if platform.machine() == "x86_64" and platform.system() == "Darwin":
        if numpy_major >= 2:
            msg.append(
                f"NumPy {np.__version__} with torch {torch.__version__} on an Intel Mac. "
                "torch 2.2.2 predates NumPy 2 and tensor.numpy() will fail. "
                "Rebuild the venv with numpy<2 pinned first; see the header of this file."
            )
    try:
        torch.zeros(2).numpy()
    except Exception as exc:  # pragma: no cover
        msg.append(f"torch cannot convert tensors to NumPy arrays: {exc}")
    if msg:
        sys.exit("[!] Environment problem:\n    " + "\n    ".join(msg))
    print(f"[*] torch {torch.__version__}, numpy {np.__version__}, "
          f"{platform.machine()} {platform.system()}: OK")


def load_manifest(limit=None):
    import csv

    path = os.path.join(IN_DIR, "manifest.csv")
    if not os.path.exists(path):
        sys.exit(f"[!] {path} not found. Run 18_hyphaeon_prep.Rscript first.")
    rows = [r for r in csv.DictReader(open(path)) if r["keep"].upper() == "TRUE"]
    if limit:
        rows = rows[:limit]
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-p", "--phenotype", choices=list(PHENOTYPES), required=True,
                    help="which phenotype to screen")
    ap.add_argument("--permulations", type=int, default=1000,
                    help="Brownian-motion permulations for the RERconverge-style gene "
                         "p-value (default 1000)")
    ap.add_argument("--alpha", type=float, default=0.05, help="site-level FDR threshold")
    ap.add_argument("--min-taxa", type=int, default=4, help="minimum taxa per site")
    ap.add_argument("--limit", type=int, default=None, help="run only the first N genes")
    ap.add_argument("--resume", action="store_true",
                    help="skip genes whose output JSON already exists")
    args = ap.parse_args()

    check_environment()

    from hyphaeon.phenotype import run_phenotype_association

    trait_col, continuous = PHENOTYPES[args.phenotype]
    pheno_file = os.path.join(IN_DIR, "phenotypes.csv")
    out_dir = os.path.join(OUT_ROOT, args.phenotype)
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(OUT_ROOT, f"progress_{args.phenotype}.log")

    genes = load_manifest(args.limit)
    print(f"[*] phenotype '{args.phenotype}' (column {trait_col}, "
          f"{'continuous' if continuous else 'binary'})")
    print(f"[*] {len(genes)} genes, {args.permulations} permulations each")
    print(f"[*] writing to {out_dir}")

    done = failed = skipped = 0
    t_start = time.time()
    with open(log_path, "a") as log:
        for i, g in enumerate(genes, 1):
            dest = os.path.join(out_dir, g["gene"] + ".json")
            if args.resume and os.path.exists(dest):
                skipped += 1
                continue
            t0 = time.time()
            try:
                res = run_phenotype_association(
                    alignment_path=g["alignment"],
                    tree_path=g["tree"],
                    weights_path=None,
                    variant=None,
                    preset=None,
                    foreground=None,
                    background=None,
                    phenotype_file=pheno_file,
                    trait_col=trait_col,
                    species_col="species",
                    continuous=continuous,
                    permulations=args.permulations,
                    min_taxa_per_site=args.min_taxa,
                    alpha=args.alpha,
                    cpu=True,
                    use_tn93=False,
                    n_permutations=10000,
                    max_perm_p=None,
                )
                res["gene"] = g["gene"]
                with open(dest, "w") as fh:
                    json.dump(res, fh)
                done += 1
                status = f"p={res.get('gene_p_value_perm')}"
            except Exception as exc:
                failed += 1
                status = f"FAILED {type(exc).__name__}: {exc}"
            log.write(f"{g['gene']}\t{time.time() - t0:.2f}\t{status}\n")
            log.flush()

            if i % 50 == 0 or i == len(genes):
                el = time.time() - t_start
                rate = (done + failed) / el if el else 0
                left = (len(genes) - i) / rate / 3600 if rate else float("nan")
                print(f"  {i}/{len(genes)}  done={done} failed={failed} "
                      f"skipped={skipped}  {rate:.2f} genes/s  ETA {left:.1f} h",
                      flush=True)

    el = time.time() - t_start
    print(f"\n[*] finished in {el / 3600:.2f} h: {done} written, {failed} failed, "
          f"{skipped} skipped")
    print(f"[*] per-gene log: {log_path}")
    if failed:
        print("[!] rerun with --resume to retry the failures")
    print("[*] next: scripts/20_hyphaeon_collect.Rscript")


if __name__ == "__main__":
    main()
