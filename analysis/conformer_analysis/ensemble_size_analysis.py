#!/usr/bin/env python3
"""
Ensemble-size analysis for gEDES-HM (JCIM ci-2026-02521e, Reviewer 1 comment 6).

For each complex, counts how many of the 101 gEDES conformers are "productive",
i.e. yield at least one pose with RMSD_BM <= threshold. The probability that a
randomly drawn sub-ensemble of size k contains at least one productive conformer
is then exact (hypergeometric), so no Monte Carlo is needed:

    P(success | k) = 1 - C(101 - m, k) / C(101, k)

Run from the directory holding the per-complex HADDOCK3 output
(inside this archive: docking/; see analysis/README.md)

    python3 ensemble_size_analysis.py                 # threshold 2.5 A
    python3 ensemble_size_analysis.py --threshold 3.0
"""

import argparse
import glob
import os
import sys

import pandas as pd
from scipy.special import comb

KS = [1, 2, 3, 5, 10, 15, 25, 50, 75, 101]


def per_complex(run_dir, threshold):
    """Return (n_conformers, n_productive, min_rmsd_bm) for one run directory."""
    merged = os.path.join(run_dir, "merged_rmsds.tsv")
    p2c = os.path.join(run_dir, "pose2conf.tsv")
    if not (os.path.isfile(merged) and os.path.isfile(p2c)):
        return None

    df = pd.read_csv(merged, sep="\t")
    conf = pd.read_csv(p2c, sep="\t", header=None, names=["pose", "conformer"])

    if "pose" not in df.columns or "combined_bs_lig_rmsd" not in df.columns:
        sys.exit(f"unexpected columns in {merged}: {list(df.columns)}")

    df = df.merge(conf, on="pose", how="inner")
    if df.empty:
        sys.exit(f"join produced no rows for {run_dir} - check pose id formatting")

    best = df.groupby("conformer")["combined_bs_lig_rmsd"].min()
    return len(best), int((best <= threshold).sum()), float(best.min())


def p_success(n, m, k):
    """Probability that k conformers drawn without replacement include >=1 productive."""
    if m == 0:
        return 0.0
    if k >= n - m + 1:
        return 1.0
    return 1.0 - comb(n - m, k, exact=True) / comb(n, k, exact=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=2.5,
                    help="RMSD_BM success threshold in Angstrom (default 2.5)")
    ap.add_argument("--glob", default="haddock_*",
                    help="run directory pattern (default haddock_*)")
    ap.add_argument("-o", "--out-prefix", default="ensemble_size",
                    help="prefix for output tsv files")
    args = ap.parse_args()

    rows = []
    for d in sorted(glob.glob(args.glob)):
        if not os.path.isdir(d):
            continue
        res = per_complex(d, args.threshold)
        if res is None:
            print(f"skipping {d} (missing merged_rmsds.tsv or pose2conf.tsv)", file=sys.stderr)
            continue
        n, m, best = res
        rows.append({"complex": d, "n_conformers": n, "productive": m,
                     "fraction_productive": m / n, "best_rmsd_bm": best})

    if not rows:
        sys.exit("no complexes processed")

    per = pd.DataFrame(rows).sort_values("productive", ascending=False)

    curve = []
    for k in KS:
        probs = [p_success(r["n_conformers"], r["productive"], k) for _, r in per.iterrows()]
        curve.append({"k": k,
                      "expected_successes": sum(probs),
                      "expected_rate": sum(probs) / len(probs)})
    curve = pd.DataFrame(curve)

    tag = f"{args.threshold:g}A".replace(".", "p")
    per_f = f"{args.out_prefix}_per_complex_{tag}.tsv"
    cur_f = f"{args.out_prefix}_curve_{tag}.tsv"
    per.to_csv(per_f, sep="\t", index=False, float_format="%.4f")
    curve.to_csv(cur_f, sep="\t", index=False, float_format="%.4f")

    n_tot = len(per)
    print(f"\nthreshold: {args.threshold} A over {n_tot} complexes\n")
    print("productive conformers per complex (out of 101)")
    print(per.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nexpected success rate vs ensemble size (out of {n_tot})")
    print(curve.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nwritten: {per_f}, {cur_f}")
    print("sanity check: expected_successes at k=101 must equal the number of "
          "complexes with productive > 0")


if __name__ == "__main__":
    main()
