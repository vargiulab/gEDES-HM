#!/usr/bin/env python3
"""
 unk_symswap.py — C2-symmetry-swapped reference for the HADDOCK AP5 ligand
(resname UNK). Same logic as ap5_symswap.py; only the atom-name map differs
because HADDOCK renamed/reordered every atom. C2 axis runs through the
central phosphate P3. Mapping derived from connectivity (traced from
coordinates), not from the labels.

Usage:  python ap5_symswap_haddock.py UNK_ref1.pdb UNK_ref2.pdb
"""
import sys, math

# ---- C2 swap map (involution) ----------------------------------------
# adenine half-1  <->  adenine half-2  (by chemical position N9/C4/C8/N7/
# C5/N1/C2/N3/C6/N6-amine)
SWAP = {
    'N1': 'N6',  'N6': 'N1',     # N9 glycosidic
    'C8': 'C16', 'C16': 'C8',    # C8
    'N2': 'N7',  'N7': 'N2',     # N7
    'C9': 'C17', 'C17': 'C9',    # C5
    'C12': 'C20','C20': 'C12',   # C4
    'C10': 'C18','C18': 'C10',   # C6
    'N3': 'N8',  'N8': 'N3',     # N6 (exocyclic amine)
    'N4': 'N9',  'N9': 'N4',     # N1 (ring)
    'C11': 'C19','C19': 'C11',   # C2
    'N5': 'N10', 'N10': 'N5',    # N3 (ring)
    # ribose-1 <-> ribose-2 (C1'/C2'/O2'/C3'/O3'/C4'/O4'ring/C5'/O5')
    'C7': 'C15', 'C15': 'C7',    # C1'
    'C13': 'C1', 'C1': 'C13',    # C2'
    'O20': 'O1', 'O1': 'O20',    # O2'
    'C14': 'C2', 'C2': 'C14',    # C3'
    'O21': 'O2', 'O2': 'O21',    # O3'
    'C6': 'C3',  'C3': 'C6',     # C4'
    'O19': 'O22','O22': 'O19',   # O4' (ring)
    'C5': 'C4',  'C4': 'C5',     # C5'
    'O18': 'O3', 'O3': 'O18',    # O5'
    # phosphates: P5<->P1, P4<->P2, P3 fixed
    'P5': 'P1',  'P1': 'P5',
    'P4': 'P2',  'P2': 'P4',
    'P3': 'P3',
    # bridging O: O15(P5-P4)<->O6(P2-P1); O12(P4-P3)<->O9(P3-P2)
    'O15': 'O6', 'O6': 'O15',
    'O12': 'O9', 'O9': 'O12',
    # non-bridging O (label-based deterministic pairing)
    'O16': 'O4', 'O4': 'O16',    # P5 / P1
    'O17': 'O5', 'O5': 'O17',
    'O13': 'O7', 'O7': 'O13',    # P4 / P2
    'O14': 'O8', 'O8': 'O14',
    'O10': 'O11','O11': 'O10',   # P3 central
}


def parse_atoms(path):
    lines, by_name, xyz = [], {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                name = line[12:16].strip()
                lines.append(line)
                by_name[name] = line
                xyz[name] = (float(line[30:38]), float(line[38:46]),
                             float(line[46:54]))
    return lines, by_name, xyz


def verify(by_name, xyz):
    # involution
    for a, b in SWAP.items():
        assert SWAP.get(b) == a, f"map not involutive at {a}->{b}"
    # all atoms covered exactly once
    names = set(by_name)
    assert names == set(SWAP), (
        f"name mismatch: only-in-pdb={names - set(SWAP)} "
        f"only-in-map={set(SWAP) - names}")
    # bond-graph preservation: build bonds by distance, check swap maps
    # every bonded pair to a bonded pair
    ns = list(names)
    def d(i, j):
        return math.dist(xyz[i], xyz[j])
    bonds = {(i, j) for x, i in enumerate(ns) for j in ns[x + 1:]
             if d(i, j) < 1.85}
    bset = {frozenset(b) for b in bonds}
    for i, j in bonds:
        si, sj = SWAP[i], SWAP[j]
        assert frozenset((si, sj)) in bset, (
            f"bond {i}-{j} not preserved under swap ({si}-{sj})")
    return len(bonds)


def write_swapped(in_path, out_path):
    lines, by_name, xyz = parse_atoms(in_path)
    nb = verify(by_name, xyz)
    with open(out_path, 'w') as f:
        for line in lines:
            name = line[12:16].strip()
            src = by_name[SWAP[name]]
            f.write(line[:30] + src[30:54] + line[54:])
        f.write("END\n")
    print(f"OK: {len(by_name)} atoms, {nb} bonds, "
          f"map verified involutive + bond-preserving")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("usage: unk_symswap.py <in.pdb> <out.pdb>")
    write_swapped(sys.argv[1], sys.argv[2])