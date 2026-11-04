"""Collective coordinates of the ADK open-to-closed transition.

Distances (Table S3 / Figure 2A): mass-weighted centres of mass over ALL atoms,
    LID 118-160, NMP 30-67, CORE 1-29, 68-117, 161-214.
Angles (Li, Liu & Ji, Biophys J 2015, 109, 647 / Figure 2B): over BACKBONE atoms,
    LID 123-155, NMP 50-59, HINGE 161-165, CORE 1-8, 79-85, 104-110, 190-198.
    theta1 = angle CORE-[HINGE]-LID, theta2 = angle HINGE-[CORE]-NMP.
All residue ranges in E. coli (1AKE/4AKE) numbering.
"""
import gzip
import numpy as np

MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'S': 32.06, 'P': 30.974, 'M': 24.305}

D_LID, D_NMP = [(118, 160)], [(30, 67)]
D_CORE = [(1, 29), (68, 117), (161, 214)]

LID, NMP, HINGE = [(123, 155)], [(50, 59)], [(161, 165)]
CORE = [(1, 8), (79, 85), (104, 110), (190, 198)]

BACKBONE = ('N', 'CA', 'C', 'O')

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
       'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
       'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
       'HID':'H','HIE':'H','HIP':'H','MSE':'M','SEC':'U','PYL':'O'}


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_pdb(path, chain=None, model=None):
    """Return {resnum: {atomname: (xyz, element)}} for one model."""
    res, cur = {}, 0
    with _open(path) as fh:
        for ln in fh:
            if ln.startswith('MODEL'):
                cur = int(ln.split()[1])
                continue
            if ln.startswith('ENDMDL'):
                if model is not None and cur == model:
                    break
                continue
            if not ln.startswith('ATOM'):
                continue
            if model is not None and cur != model:
                continue
            if ln[16] not in (' ', 'A'):          # altloc
                continue
            if chain and ln[21] != chain:
                continue
            name = ln[12:16].strip()
            elem = ln[76:78].strip() or name[0]
            try:
                num = int(ln[22:26])
            except ValueError:
                continue
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            res.setdefault(num, {})[name] = (xyz, elem, ln[17:20].strip())
    return res


def iter_models(path, chain=None):
    """Yield (model_number, residues) for every MODEL in a multi-model PDB."""
    res, cur = {}, None
    with _open(path) as fh:
        for ln in fh:
            if ln.startswith('MODEL'):
                cur = int(ln.split()[1]); res = {}
                continue
            if ln.startswith('ENDMDL'):
                yield cur, res
                res = {}
                continue
            if not ln.startswith('ATOM'):
                continue
            if ln[16] not in (' ', 'A'):
                continue
            if chain and ln[21] != chain:
                continue
            name = ln[12:16].strip()
            elem = ln[76:78].strip() or name[0]
            try:
                num = int(ln[22:26])
            except ValueError:
                continue
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            res.setdefault(num, {})[name] = (xyz, elem, ln[17:20].strip())
    if res:
        yield cur, res


def _com(res, ranges, backbone_only=False, keymap=None):
    xs, ws, have, tot = [], [], 0, 0
    for lo, hi in ranges:
        for r in range(lo, hi + 1):
            tot += 1
            k = keymap.get(r) if keymap is not None else r
            if k is None or k not in res:
                continue
            have += 1
            for nm, (xyz, elem, _rn) in res[k].items():
                if backbone_only and nm not in BACKBONE:
                    continue
                xs.append(xyz)
                ws.append(MASS.get(elem[0], 12.011))
    if not xs:
        return None, have / tot if tot else 0.0
    xs, ws = np.array(xs), np.array(ws)
    return (xs * ws[:, None]).sum(0) / ws.sum(), have / tot


def frac(res, ranges, keymap=None):
    tot = sum(hi - lo + 1 for lo, hi in ranges)
    have = 0
    for lo, hi in ranges:
        for r in range(lo, hi + 1):
            k = keymap.get(r) if keymap is not None else r
            if k is not None and k in res:
                have += 1
    return have / tot


def _angle(a, vertex, b):
    if a is None or vertex is None or b is None:
        return float('nan')
    u, v = a - vertex, b - vertex
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def coords(res, keymap=None):
    """Return dict with dLID, dNMP, theta1, theta2 and coverage fractions."""
    core_d, _ = _com(res, D_CORE, keymap=keymap)
    lid_d, cov_dlid = _com(res, D_LID, keymap=keymap)
    nmp_d, _ = _com(res, D_NMP, keymap=keymap)
    core_a, _ = _com(res, CORE, backbone_only=True, keymap=keymap)
    lid_a, _ = _com(res, LID, backbone_only=True, keymap=keymap)
    nmp_a, _ = _com(res, NMP, backbone_only=True, keymap=keymap)
    hng_a, _ = _com(res, HINGE, backbone_only=True, keymap=keymap)
    return dict(
        dLID=float(np.linalg.norm(lid_d - core_d)) if lid_d is not None and core_d is not None else float('nan'),
        dNMP=float(np.linalg.norm(nmp_d - core_d)) if nmp_d is not None and core_d is not None else float('nan'),
        theta1=_angle(core_a, hng_a, lid_a),
        theta2=_angle(hng_a, core_a, nmp_a),
        cov_dist_LID=cov_dlid,
        cov_ang_LID=frac(res, LID, keymap),
        cov_ang_CORE=frac(res, CORE, keymap),
        cov_ang_NMP=frac(res, NMP, keymap),
        cov_ang_HINGE=frac(res, HINGE, keymap),
    )
