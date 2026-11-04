#!/usr/bin/env python3
"""
map_to_reference.py

Align a homologue structure's sequence to the E. coli reference by
Needleman-Wunsch (BLOSUM62, affine gaps), renumber it into reference
numbering, and emit the list of reference residues that have a one-to-one
structural partner.

Needed because the Aquifex aeolicus TSE construct (PDB 3SR0, 201 residues)
cannot be mapped onto E. coli numbering (214 residues) by a rigid offset:
there are insertions and deletions, and a shift that looks right in the P-loop
will silently misalign the LID.

Dependencies: numpy only.

Usage
-----
    python3 map_to_reference.py \
        --reference pdb/1ake.pdb \
        --target    tse_msmd.pdb \
        --out       tse_msmd_ecnum.pdb \
        --columns-out common_columns.txt \
        --report    alignment_report.txt

Handles multi-model PDBs: the sequence is taken from model 1 and the
renumbering applied to every model. Target residues with no aligned partner
are dropped, so the output contains only residues that map one-to-one.
"""

import argparse
import sys
import numpy as np

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
    "ASH": "D", "GLH": "E", "LYN": "K", "ACE": "X", "NME": "X",
    "NHE": "X",
}

_B62_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62_ROWS = """
 4-1-2-2 0-1-1 0-2-1-1-1-1-2-1 1 0-3-2 0-2-1 0-4
-1 5 0-2-3 1 0-2 0-3-2 2-1-3-2-1-1-3-2-3-1 0-1-4
-2 0 6 1-3 0 0 0 1-3-3 0-2-3-2 1 0-4-2-3 3 0-1-4
-2-2 1 6-3 0 2-1-1-3-4-1-3-3-1 0-1-4-3-3 4 1-1-4
 0-3-3-3 9-3-4-3-3-1-1-3-1-2-3-1-1-2-2-1-3-3-2-4
-1 1 0 0-3 5 2-2 0-3-2 1 0-3-1 0-1-2-1-2 0 3-1-4
-1 0 0 2-4 2 5-2 0-3-3 1-2-3-1 0-1-3-2-2 1 4-1-4
 0-2 0-1-3-2-2 6-2-4-4-2-3-3-2 0-2-2-3-3-1-2-1-4
-2 0 1-1-3 0 0-2 8-3-3-1-2-1-2-1-2-2 2-3 0 0-1-4
-1-3-3-3-1-3-3-4-3 4 2-3 1 0-3-2-1-3-1 3-3-3-1-4
-1-2-3-4-1-2-3-4-3 2 4-2 2 0-3-2-1-2-1 1-4-3-1-4
-1 2 0-1-3 1 1-2-1-3-2 5-1-3-1 0-1-3-2-2 0 1-1-4
-1-1-2-3-1 0-2-3-2 1 2-1 5 0-2-1-1-1-1 1-3-1-1-4
-2-3-3-3-2-3-3-3-1 0 0-3 0 6-4-2-2 1 3-1-3-3-1-4
-1-2-2-1-3-1-1-2-2-3-3-1-2-4 7-1-1-4-3-2-2-1-2-4
 1-1 1 0-1 0 0 0-1-2-2 0-1-2-1 4 1-3-2-2 0 0 0-4
 0-1 0-1-1-1-1-2-2-1-1-1-1-2-1 1 5-2-2 0-1-1 0-4
-3-3-4-4-2-2-3-2-2-3-2-3-1 1-4-3-211 2-3-4-3-2-4
-2-2-2-3-2-1-2-3 2-1-1-2-1 3-3-2-2 2 7-1-3-2-1-4
 0-3-3-3-1-2-2-3-3 3 1-2 1-1-2-2 0-3-1 4-3-2-1-4
-2-1 3 4-3 0 1-1 0-3-4 0-3-3-2 0-1-4-3-3 4 1-1-4
-1 0 0 1-3 3 4-2 0-3-3 1-1-3-1 0-1-3-2-2 1 4-1-4
 0-1-1-1-2-1-1-1-1-1-1-1-1-1-2 0 0-2-1-1-1-1-1-4
-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4-4 1
"""


def blosum62():
    rows = [r for r in _B62_ROWS.strip("\n").splitlines()]
    m = {}
    for i, row in enumerate(rows):
        vals = [int(row[j:j + 2]) for j in range(0, 2 * len(_B62_ORDER), 2)]
        for j, v in enumerate(vals):
            m[(_B62_ORDER[i], _B62_ORDER[j])] = v
    return m


B62 = blosum62()


def score(a, b):
    return B62.get((a, b), B62.get((b, a), -1))


def nw_affine(s1, s2, gap_open=-11, gap_extend=-1):
    """Global alignment, Gotoh. Returns list of (i, j) index pairs, 0-based,
    for aligned positions (gaps excluded)."""
    n, m = len(s1), len(s2)
    NEG = -1e9
    M = np.full((n + 1, m + 1), NEG)
    X = np.full((n + 1, m + 1), NEG)   # gap in s2
    Y = np.full((n + 1, m + 1), NEG)   # gap in s1
    ptrM = np.zeros((n + 1, m + 1), dtype=np.int8)
    ptrX = np.zeros((n + 1, m + 1), dtype=np.int8)
    ptrY = np.zeros((n + 1, m + 1), dtype=np.int8)

    M[0, 0] = 0.0
    for i in range(1, n + 1):
        X[i, 0] = gap_open + (i - 1) * gap_extend
        ptrX[i, 0] = 1
    for j in range(1, m + 1):
        Y[0, j] = gap_open + (j - 1) * gap_extend
        ptrY[0, j] = 2

    for i in range(1, n + 1):
        a = s1[i - 1]
        for j in range(1, m + 1):
            b = s2[j - 1]
            sc = score(a, b)
            prev = (M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
            k = int(np.argmax(prev))
            M[i, j] = prev[k] + sc
            ptrM[i, j] = k

            cand = (M[i - 1, j] + gap_open, X[i - 1, j] + gap_extend,
                    Y[i - 1, j] + gap_open)
            k = int(np.argmax(cand))
            X[i, j] = cand[k]
            ptrX[i, j] = k

            cand = (M[i, j - 1] + gap_open, X[i, j - 1] + gap_open,
                    Y[i, j - 1] + gap_extend)
            k = int(np.argmax(cand))
            Y[i, j] = cand[k]
            ptrY[i, j] = k

    state = int(np.argmax((M[n, m], X[n, m], Y[n, m])))
    i, j = n, m
    pairs = []
    while i > 0 or j > 0:
        if state == 0:
            pairs.append((i - 1, j - 1))
            state = int(ptrM[i, j])
            i, j = i - 1, j - 1
        elif state == 1:
            state = int(ptrX[i, j])
            i -= 1
        else:
            state = int(ptrY[i, j])
            j -= 1
    pairs.reverse()
    return pairs


def read_ca_seq(path, chain="-"):
    """First model only. Returns (resnums, one_letter_string)."""
    nums, seq = [], []
    want = None
    with open(path) as fh:
        for line in fh:
            if line[:6] == "ENDMDL":
                break
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue
            ch = line[21]
            if chain == "-":
                if want is None:
                    want = ch
                if ch != want:
                    continue
            elif ch != chain:
                continue
            if line[26] != " ":
                continue
            rn = line[17:20].strip()
            one = THREE2ONE.get(rn, "X")
            if one == "X":
                continue                      # caps and unknowns
            nums.append(int(line[22:26]))
            seq.append(one)
    return nums, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--columns-out", default="common_columns.txt")
    ap.add_argument("--report", default="alignment_report.txt")
    ap.add_argument("--ref-chain", default="-")
    ap.add_argument("--target-chain", default="-")
    args = ap.parse_args()

    rnums, rseq = read_ca_seq(args.reference, args.ref_chain)
    tnums, tseq = read_ca_seq(args.target, args.target_chain)
    sys.stderr.write("reference %d residues, target %d residues\n"
                     % (len(rseq), len(tseq)))

    pairs = nw_affine(rseq, tseq)
    ident = sum(1 for i, j in pairs if rseq[i] == tseq[j])
    pct = 100.0 * ident / max(1, len(pairs))
    sys.stderr.write("aligned %d positions, %d identical (%.1f%%)\n"
                     % (len(pairs), ident, pct))

    # target residue number -> reference residue number
    t2r = {tnums[j]: rnums[i] for i, j in pairs}
    common = sorted(rnums[i] for i, j in pairs)

    # ---- report -------------------------------------------------------
    with open(args.report, "w") as fh:
        fh.write("reference : %s (%d residues)\n" % (args.reference, len(rseq)))
        fh.write("target    : %s (%d residues)\n" % (args.target, len(tseq)))
        fh.write("aligned   : %d positions, %d identical (%.1f%%)\n\n"
                 % (len(pairs), ident, pct))
        fh.write("reference residues WITHOUT a partner (%d):\n"
                 % (len(rnums) - len(pairs)))
        miss = sorted(set(rnums) - set(common))
        fh.write("  %s\n\n" % (miss if miss else "none"))
        fh.write("ref_resid ref_aa target_resid target_aa identical\n")
        for i, j in pairs:
            fh.write("%9d %6s %12d %9s %9s\n"
                     % (rnums[i], rseq[i], tnums[j], tseq[j],
                        "yes" if rseq[i] == tseq[j] else "no"))

    open(args.columns_out, "w").write(
        " ".join(str(r) for r in common) + "\n")

    # ---- renumber every model ------------------------------------------
    nwritten = 0
    with open(args.target) as fin, open(args.out, "w") as fout:
        for line in fin:
            rec = line[:6]
            if rec in ("MODEL ", "ENDMDL", "TER   ", "END   "):
                fout.write(line)
                if rec == "ENDMDL":
                    nwritten += 1
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue
            try:
                old = int(line[22:26])
            except ValueError:
                continue
            if old not in t2r:
                continue
            fout.write(line[:22] + "%4d" % t2r[old] + line[26:])
    sys.stderr.write("wrote %s (%d models), %s, %s\n"
                     % (args.out, max(nwritten, 1), args.columns_out, args.report))


if __name__ == "__main__":
    main()
