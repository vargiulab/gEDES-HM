"""E. coli (1AKE) residue number -> homology-model residue number, per target,
taken from the MODELLER .pir alignment actually used to build the models."""
import glob
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, "common"))
from adk_geom import iter_models, AA3
from adk_paths import GEDES_MODELS as ZEN


def parse_pir(path):
    entries, name, hdr, seq = [], None, None, []
    for ln in open(path):
        ln = ln.rstrip('\n')
        if ln.startswith('>P1;'):
            if name: entries.append((name, hdr, ''.join(seq)))
            name, hdr, seq = ln[4:].strip(), None, []
        elif hdr is None and name is not None:
            hdr = ln
        elif name is not None:
            seq.append(ln.strip())
    if name: entries.append((name, hdr, ''.join(seq)))
    return [(n, h, s.replace('*', '')) for n, h, s in entries]


def keymap_for(target):
    """{ecoli_resnum: model_resnum} plus the target sequence."""
    pirs = glob.glob(f"{ZEN}/{target}/*.pir")
    assert len(pirs) == 1, (target, pirs)
    ents = parse_pir(pirs[0])
    tmpl = [e for e in ents if e[1].startswith('structure')]
    tgt  = [e for e in ents if e[1].startswith('sequence')]
    assert len(tmpl) == 1 and len(tgt) == 1, (target, [e[0] for e in ents])
    tseq, qseq = tmpl[0][2], tgt[0][2]
    assert len(tseq) == len(qseq), (target, len(tseq), len(qseq))
    # template range from the header: structure:ADK:<start>:<ch>:<end>:...
    f = tmpl[0][1].split(':')
    tstart = int(f[2])
    ti, qi, km = tstart - 1, 0, {}
    for a, b in zip(tseq, qseq):
        if a != '-': ti += 1
        if b != '-': qi += 1
        if a != '-' and b != '-':
            km[ti] = qi
    return km, qseq


def check(target, km, qseq):
    """Fraction of mapped positions whose model residue name matches the
    target sequence letter, over the first model."""
    path = f"{ZEN}/{target}/concat_selected_clean.pdb.gz"
    _, res = next(iter_models(path))
    ok = tot = 0
    present = 0
    for e, m in km.items():
        if m not in res: continue
        present += 1
        rn = next(iter(res[m].values()))[2]
        tot += 1
        if AA3.get(rn, 'X') == qseq[m - 1]: ok += 1
    return ok, tot, present, len(km), sorted(res)[0], sorted(res)[-1]


if __name__ == '__main__':
    print(f"{'tgt':>4} {'mapped':>7} {'inmodel':>8} {'nameOK':>8} {'range':>12}")
    for i in range(1, 23):
        km, q = keymap_for(str(i))
        ok, tot, present, nmap, lo, hi = check(str(i), km, q)
        flag = '' if (tot and ok == tot) else '  <-- MISMATCH'
        print(f"{i:>4} {nmap:>7} {present:>8} {ok}/{tot:<8} {lo}-{hi}{flag}")
