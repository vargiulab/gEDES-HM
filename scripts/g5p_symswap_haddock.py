#!/usr/bin/env python3
"""
Usage:  python gp5_symswap_haddock.py UNK_ref1.pdb UNK_ref2.pdb
"""
import sys
import numpy as np

# covalent radii (Angstrom); bond if dist < BOND_TOL * (r_i + r_j)
COV = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05,
       'F': 0.57, 'CL': 1.02, 'BR': 1.20}
BOND_TOL = 1.3


def parse(path):
    names, elems, xyz, lines = [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                lines.append(line)
                names.append(line[12:16].strip())
                el = line[76:78].strip().upper()
                if not el:  # fall back to atom name
                    el = ''.join(c for c in line[12:16] if c.isalpha())[:2].upper()
                elems.append(el)
                xyz.append((float(line[30:38]), float(line[38:46]),
                            float(line[46:54])))
    return names, elems, np.array(xyz), lines


def bond_graph(elems, xyz):
    n = len(elems)
    adj = [[] for _ in range(n)]
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    for i in range(n):
        ri = COV.get(elems[i], 0.77)
        for j in range(i + 1, n):
            cut = BOND_TOL * (ri + COV.get(elems[j], 0.77))
            if 0.4 < d[i, j] < cut:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def automorphisms(elems, adj, limit=20000):
    """Backtracking search for element-respecting graph automorphisms."""
    n = len(elems)
    deg = [len(a) for a in adj]
    sig = [(elems[i], deg[i]) for i in range(n)]
    # candidate targets per node, pruned by (element, degree)
    cand = [[j for j in range(n) if sig[j] == sig[i]] for i in range(n)]
    order = sorted(range(n), key=lambda i: len(cand[i]))
    maps = []
    perm = [-1] * n
    used = [False] * n

    def bt(k):
        if len(maps) >= limit:
            return
        if k == n:
            maps.append(perm.copy())
            return
        i = order[k]
        for j in cand[i]:
            if used[j]:
                continue
            ok = True
            for nb in adj[i]:
                pn = perm[nb]
                if pn != -1 and pn not in adj[j]:
                    ok = False
                    break
            if not ok:
                continue
            # neighbor-count consistency
            if len(adj[i]) != len(adj[j]):
                continue
            perm[i] = j
            used[j] = True
            bt(k + 1)
            perm[i] = -1
            used[j] = False

    bt(0)
    return maps


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1, 1, d])
    R = V @ D @ Wt
    diff = Pc @ R - Qc
    return np.sqrt((diff ** 2).sum() / len(P))


def pick_c2(elems, xyz, autos):
    """Among non-identity automorphisms, choose the one whose permuted
    coordinates best superpose on the original (the geometric C2)."""
    n = len(elems)
    best, best_rmsd, best_moved = None, np.inf, 0
    for p in autos:
        if p == list(range(n)):
            continue
        moved = sum(1 for i in range(n) if p[i] != i)
        r = kabsch_rmsd(xyz, xyz[p])
        # prefer the deep, well-superposing involution (the true C2)
        is_inv = all(p[p[i]] == i for i in range(n))
        score = r + (0 if is_inv else 5.0) + (0 if moved > n // 2 else 5.0)
        if score < best_rmsd:
            best_rmsd, best, best_moved = score, p, moved
    return best, best_moved


def main(inp, outp):
    names, elems, xyz, lines = parse(inp)
    n = len(names)
    adj = bond_graph(elems, xyz)
    nb = sum(len(a) for a in adj) // 2
    autos = automorphisms(elems, adj)
    if not autos:
        raise SystemExit("no automorphisms found (check bond detection)")
    perm, moved = pick_c2(elems, xyz, autos)
    if perm is None:
        raise SystemExit("only identity automorphism found - not symmetric?")
    inv_ok = all(perm[perm[i]] == i for i in range(n))
    fixed = [names[i] for i in range(n) if perm[i] == i]
    # write ref2: atom i keeps its line/name, takes coords of perm[i]
    with open(outp, 'w') as f:
        for i, line in enumerate(lines):
            src = lines[perm[i]]
            f.write(line[:30] + src[30:54] + line[54:])
        f.write("END\n")
    rms = kabsch_rmsd(xyz, xyz[perm])
    print(f"OK: {n} atoms, {nb} bonds, {len(autos)} automorphism(s)")
    print(f"    chosen C2: involution={inv_ok}, atoms moved={moved}, "
          f"fixed={fixed if fixed else 'none'}")
    print(f"    self-superposition RMSD under C2 = {rms:.3f} A "
          f"(structural asymmetry; small = good)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("usage: ap5_symswap_haddock.py <in.pdb> <out.pdb>")
    main(sys.argv[1], sys.argv[2])