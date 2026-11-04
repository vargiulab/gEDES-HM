#!/usr/bin/env python
"""
apo homology modeling (single 4AKE template) -- lean pipeline + safe trim (v14)

= v13 (light refinement, template-faithful, fast) PLUS terminal trimming done
the safe way:

  1. Trim terminal target overhangs from the alignment so MODELLER never builds
     the spaghetti N/C tail (only columns where the TEMPLATE is gapped; internal
     gaps untouched).
  2. After building, shift output residue numbers back by the number of
     N-terminal residues that were trimmed, so the model carries FULL-LENGTH
     target numbering. Your cpptraj dual-reference RMSD masks keep working
     unchanged -- this is what broke in v12 and is now fixed.

Nothing else about the geometry changes: same light refinement that gave you
good BS / complex RMSD. The trim only removes residues that had zero template
support and were excluded from RMSD anyway.

Assumes the target's full-length numbering starts at TARGET_FIRST_RESNUM
(UniProt-style, usually 1). Adjust if your target starts elsewhere.
"""

import os
import glob
import shutil

from modeller import Environ, Alignment, log
from modeller.automodel import AutoModel, refine, assess

# Single apo template (PDB 4AKE), relative to this system folder in the archive.
# The archived template is gzipped; decompress it in place
# ('gunzip ../reference_apo_template_4ake.pdb.gz') before running.
TEMPLATE_DIR      = '../'
TEMPLATES         = ['../reference_apo_template_4ake.pdb']
N_TEMPLATES       = len(TEMPLATES)
RAW_PIR           = '9_adk_Q8IJV6.pir'   # YOUR alignment, used as is (no realign)
KNOWN             = 'ADK'
TARGET            = 'Q8IJV6'
TARGET_FIRST_RESNUM = 1                    # first residue number of full target

ALN_PIR           = 'aln_trimmed.pir'
ALN_PAP           = 'aln_trimmed.pap'
MODELS_PER_TMPL   = 10
MD_LEVEL          = refine.slow            # try very_fast/none for more speed

log.verbose()


# ----------------------------------------------------------------------------
def _read_pir(path):
    with open(path) as fh:
        text = fh.read()
    recs, order = {}, []
    for block in text.split('>P1;')[1:]:
        lines = block.splitlines()
        code = lines[0].strip()
        recs[code] = [lines[1], ''.join(lines[2:]).replace('*', '')]
        order.append(code)
    return recs, order


def _write_pir(path, recs, order):
    with open(path, 'w') as fh:
        for code in order:
            desc, seq = recs[code]
            fh.write('>P1;%s\n%s\n' % (code, desc))
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')
            fh.seek(fh.tell() - 1)
            fh.write('*\n')


def trim_termini(in_pir, out_pir):
    """Strip leading/trailing columns where the TEMPLATE is gapped.

    Returns (n_trim_N, n_trim_C): target residues removed at each terminus.
    """
    recs, order = _read_pir(in_pir)
    T, Q = recs[KNOWN][1], recs[TARGET][1]
    assert len(T) == len(Q), 'aligned lengths differ; check the .pir'
    L = len(T)

    start = 0
    while start < L and T[start] == '-':
        start += 1
    end = L
    while end > start and T[end - 1] == '-':
        end -= 1

    n_trim_N = sum(1 for c in Q[:start] if c != '-')
    n_trim_C = sum(1 for c in Q[end:]  if c != '-')

    recs[KNOWN][1]  = T[start:end]
    recs[TARGET][1] = Q[start:end]
    _write_pir(out_pir, recs, order)
    print('Trimmed %d N-terminal and %d C-terminal target residue(s).'
          % (n_trim_N, n_trim_C))
    return n_trim_N, n_trim_C


def renumber_pdb(pdb_file, offset):
    """Shift every residue number by offset (PDB cols 23-26), in place."""
    if offset == 0:
        return
    with open(pdb_file) as fh:
        lines = fh.readlines()
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM', 'TER')) and len(line) >= 26:
            try:
                new = int(line[22:26]) + offset
                line = line[:22] + '%4d' % new + line[26:]
            except ValueError:
                pass
        out.append(line)
    with open(pdb_file, 'w') as fh:
        fh.writelines(out)


def audit_internal_gaps(pir_file, min_len=10):
    recs, _ = _read_pir(pir_file)
    T, Q = recs[KNOWN][1], recs[TARGET][1]
    tgt_res, in_gap, glen, start = 0, False, 0, 0
    for t, q in zip(T, Q):
        if q != '-':
            tgt_res += 1
            if in_gap:
                if glen >= min_len:
                    print('  WARNING: %d-residue internal gap after target '
                          'residue %d. Real missing subdomain or alignment '
                          'error? Check the .pap.' % (glen, start))
                in_gap, glen = False, 0
        elif t != '-':
            if not in_gap:
                in_gap, start, glen = True, tgt_res, 0
            glen += 1


# ----------------------------------------------------------------------------
def main():
    env = Environ()
    env.io.atom_files_directory = ['.', TEMPLATE_DIR]

    n_trim_N, _ = trim_termini(RAW_PIR, ALN_PIR)
    offset = (TARGET_FIRST_RESNUM - 1) + n_trim_N   # add to model numbers

    aln = Alignment(env)
    aln.append(file=ALN_PIR, alignment_format='PIR')
    aln.write(file=ALN_PAP, alignment_format='PAP')
    print('Inspect %s. Output models will be renumbered by +%d so pocket '
          'residues keep full-length numbering.' % (ALN_PAP, offset))
    audit_internal_gaps(ALN_PIR)

    for i, template_file in enumerate(TEMPLATES):
        if not os.path.exists(template_file):
            print('Template %s not found. Skipping.' % template_file)
            continue
        print('\n=== Template %d: %s ===' % (i, template_file))
        shutil.copy(template_file, 'ADK.pdb')

        a = AutoModel(env,
                      alnfile=ALN_PIR,
                      knowns=KNOWN,
                      sequence=TARGET,
                      assess_methods=(assess.DOPE, assess.normalized_dope))
        a.starting_model = 1
        a.ending_model = MODELS_PER_TMPL
        a.md_level = MD_LEVEL
        a.make()

        # renumber every output model back to full-length numbering
        for pdb in glob.glob('%s.*.pdb' % TARGET):
            renumber_pdb(pdb, offset)

        ok = [x for x in a.outputs if x['failure'] is None]
        if ok:
            ok.sort(key=lambda x: x['Normalized DOPE score'])
            print('  z-DOPE ranking (best first):')
            for x in ok:
                print('    %-30s zDOPE=%.3f' % (x['name'],
                                                x['Normalized DOPE score']))

        result_folder = 'results_%d' % i
        os.makedirs(result_folder, exist_ok=True)
        for f in glob.glob('%s.*' % TARGET):
            shutil.move(f, os.path.join(result_folder, os.path.basename(f)))
        shutil.copy('ADK.pdb',
                    os.path.join(result_folder, os.path.basename(template_file)))
        os.remove('ADK.pdb')


if __name__ == '__main__':
    main()