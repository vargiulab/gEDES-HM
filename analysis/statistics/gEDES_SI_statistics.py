#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Statistical analysis of gEDES homology-model docking
#
# Purpose
#   Quantify which structural quantity is associated with the binding-mode
#   accuracy of HADDOCK3 ensemble docking into a gEDES homology model
#   (native-like when RMSD_BM <= 2.5 Angstrom), across 22 adenylate-kinase
#   homologs (26 ligand poses, counting the two nucleotides of bisubstrate
#   complexes separately).
#
#   Quantities tested:
#       rmsd_BS    binding-site RMSD of the best docked complex vs reference
#       rmsd_prot  global protein-model RMSD of the best complex vs reference
#       rmsd_lig   ligand RMSD after final refinement
#       Fnat       fraction of native contacts
#       bs_id      binding-site sequence identity to the reference (P69441)
#       prot_id    overall sequence identity to the reference
#
# Method rationale: the sample is small (n = 26 poses / 22 targets), contains
# influential points (e.g. target 11, global RMSD 6.7 A), and the associations
# are not assumed linear or Gaussian. We therefore use the Spearman rank
# correlation throughout, which tests for monotonic rather than strictly linear
# association and is robust to outliers. Because the multi-ligand complexes
# contribute correlated poses, every association is reported at both the pose
# level (n = 26) and the conservative target level (n = 22, per-target means).
#
# Note that rmsd_BS and rmsd_BM are not independent: the binding-site heavy
# atoms that define rmsd_BS also enter rmsd_BM together with the ligand atoms,
# so their correlation is inflated by shared atoms and is reported as a
# consistency check rather than as an independent association.
#
# Dependencies: pandas, scipy  (Python >= 3.9)
# Usage:        python gEDES_SI_statistics.py
# Outputs:      console summary
# -----------------------------------------------------------------------------

import os

import pandas as pd
from scipy.stats import spearmanr

THRESHOLD = 2.5          # native-like binding-mode cutoff (Angstrom)
PREDICTORS = ["rmsd_BS", "rmsd_prot", "rmsd_lig", "Fnat", "bs_id", "prot_id"]


DATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "gEDES_docking_data.csv")


def load(path=None):
    """Read the per-complex table.

    The default input is gEDES_docking_data.csv, shipped beside this script:
    one row per complex (26 rows, the two nucleotides of the bisubstrate
    complexes counted separately) with the columns

        sys, label, prot_id, bs_id, rmsd_BM, rmsd_lig, Fnat, rmsd_BS, rmsd_prot

    A tab-separated table using the BM_new / BS_new / lig_new / prot_new /
    Fnat_new names is also accepted and is renamed on the fly.
    """
    path = path or DATA_FILE
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    df = df.rename(columns={"BM_new": "rmsd_BM", "BS_new": "rmsd_BS",
                            "lig_new": "rmsd_lig", "prot_new": "rmsd_prot",
                            "Fnat_new": "Fnat"})
    df["native_like"] = df["rmsd_BM"] <= THRESHOLD
    return df


def target_level(df):
    """Conservative aggregation: one independent row per target (per-target mean)."""
    agg = (df.groupby("sys")
             .agg(rmsd_BM=("rmsd_BM", "mean"), rmsd_BS=("rmsd_BS", "mean"),
                  rmsd_prot=("rmsd_prot", "mean"), rmsd_lig=("rmsd_lig", "mean"),
                  Fnat=("Fnat", "mean"), bs_id=("bs_id", "first"),
                  prot_id=("prot_id", "first"))
             .reset_index())
    agg["native_like"] = agg["rmsd_BM"] <= THRESHOLD
    return agg


def spearman_table(df, agg):
    print("Spearman correlation of RMSD_BM with each predictor")
    print(f"  {'predictor':10s}  {'pose level (n=26)':24s}  {'target level (n=22)':24s}")
    for c in PREDICTORS:
        r1, p1 = spearmanr(df["rmsd_BM"], df[c])
        r2, p2 = spearmanr(agg["rmsd_BM"], agg[c])
        print(f"  {c:10s}  rho={r1:+.2f}  p={p1:.2e}      rho={r2:+.2f}  p={p2:.2e}")


def separation(df):
    pmax = df.loc[df.native_like, "rmsd_BS"].max()
    fmin = df.loc[~df.native_like, "rmsd_BS"].min()
    print("\nSeparating interval on RMSD_BS")
    print(f"  native-like max RMSD_BS = {pmax:.3f} A ; failed min RMSD_BS = {fmin:.3f} A")
    print(f"  -> every recovered complex lies below {pmax:.3f} A and both unrecovered "
          f"ones above {fmin:.3f} A")
    print("  With only two failure events this is a descriptive boundary, not a "
          "criterion applied a priori.")


def main():
    df = load()
    agg = target_level(df)
    print(f"n poses = {len(df)}, n targets = {len(agg)}, "
          f"native-like = {df.native_like.sum()}/{len(df)} ({100*df.native_like.mean():.1f}%)\n")
    spearman_table(df, agg)
    separation(df)


if __name__ == "__main__":
    main()
