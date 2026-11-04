#!/usr/bin/env python3
"""Both collective-coordinate descriptions of the ADK transition, for the gEDES
library, the two E. coli references, the 22 homologue targets and the seven
experimental intermediates. This is the source of Table S3 and of Figure 2.

Distances  centre-of-mass separations over all atoms, with the domains of
           Table S3 (LID 118-160, NMP 30-67, CORE 1-29, 68-117, 161-214).
Angles     theta1 = CORE-[hinge]-LID, theta2 = hinge-[CORE]-NMP, as defined by
           Li, Liu & Ji, Biophys J 2015, 109, 647, over backbone centres of
           mass with LID 123-155, NMP 50-59, hinge 161-165 and
           CORE 1-8, 79-85, 104-110, 190-198.

Everything is in E. coli (1AKE/4AKE) numbering; the homologue and intermediate
structures were renumbered onto that frame beforehand with map_to_reference.py,
and are stored here as the *_ecnum.pdb files under analysis/structures/.

Inputs   analysis/structures/ecoli/{1ake,4ake}.pdb
         analysis/structures/homologues_ecnum/T<n>_ecnum.pdb        (n = 1..22)
         analysis/structures/intermediates_ecnum/<PDBID>_ecnum.pdb  (7 files)
         proteins_structural_models/gEDES/reference_gEDES_templates/<i>.pdb.gz
                                                                    (i = 0..100)
Output   adk_coordinates.tsv  (one row per structure)

Usage    python compute_adk_coordinates.py [-o adk_coordinates.tsv]

Paths are resolved by analysis/common/adk_paths.py and contain no absolute
component; set ADK_ARCHIVE_ROOT if you move the archive apart.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "common"))
import numpy as np
from adk_geom import read_pdb, coords
from adk_paths import ECOLI, HOMOLOGUES, INTERMEDIATES, GEDES_TEMPLATES, require

COLS = ["label", "role", "dLID", "dNMP", "theta1", "theta2",
        "cov_dist_LID", "cov_ang_LID", "cov_ang_CORE"]

INTERM_LID = ("1ZIN", "1ZIO", "1ZIP")   # partly open LID, NMP closed
INTERM_NMP = ("1AK2", "2AK2", "1DVR", "2C9Y")   # LID closed, NMP open


def row(label, role, path, chain=None):
    c = coords(read_pdb(require(path, f"structure for {label}"), chain=chain))
    out = dict(label=label, role=role)
    out.update({k: c[k] for k in COLS[2:]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--out", default="adk_coordinates.tsv")
    args = ap.parse_args()

    rows = [row("4AKE", "reference", os.path.join(ECOLI, "4ake.pdb")),
            row("1AKE", "reference", os.path.join(ECOLI, "1ake.pdb"))]
    for pid in INTERM_LID:
        rows.append(row(pid, "interm_lid",
                        os.path.join(INTERMEDIATES, f"{pid}_ecnum.pdb")))
    for pid in INTERM_NMP:
        rows.append(row(pid, "interm_nmp",
                        os.path.join(INTERMEDIATES, f"{pid}_ecnum.pdb")))
    for i in range(1, 23):
        p = os.path.join(HOMOLOGUES, f"T{i}_ecnum.pdb")
        if os.path.exists(p):
            rows.append(row(f"T{i:02d}", "homologue", p))
    for i in range(101):
        p = os.path.join(GEDES_TEMPLATES, f"{i}.pdb.gz")
        if not os.path.exists(p):                     # tolerate a decompressed set
            p = os.path.join(GEDES_TEMPLATES, f"{i}.pdb")
        rows.append(row(f"gEDES_{i:03d}", "gedes", p))

    with open(args.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                               for c in COLS) + "\n")
    print(f"wrote {args.out}, {len(rows)} structures")


if __name__ == "__main__":
    main()
