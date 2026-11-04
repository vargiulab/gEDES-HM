#!/usr/bin/env python
"""
gEDES homology modeling -- lean pipeline + safe trim (v15)

= v14 PLUS:
  1. TEMPLATE_PATH can now be EITHER
       (a) a directory containing numbered cluster templates 0.pdb, 1.pdb, ...
           (old behaviour, used when N_TEMPLATES > 1 and the path is a dir), OR
       (b) a single PDB file to use directly as the one-and-only template
           (new: fixes "Template .../foo.pdb/0.pdb not found" when you just
           want to model off one structure, no clustering involved).
     Detection is automatic based on os.path.isfile()/isdir().
  2. All output now goes into a single results/ folder (results/template_%d/
     sub-folders are only created when there is genuinely more than one
     cluster template; for the single-template case everything lands flat
     in results/).

Everything else (terminal trimming, renumbering back to full-length
numbering, DOPE ranking, internal-gap audit) is unchanged from v14.

Assumes the target's full-length numbering starts at TARGET_FIRST_RESNUM
(UniProt-style, usually 1). Adjust if your target starts elsewhere.
"""

import os
import glob
import shutil

from modeller import Environ, Alignment, log
from modeller.automodel import AutoModel, refine, assess

# ----------------------------------------------------------------------------
# EITHER a directory of numbered cluster templates (0.pdb, 1.pdb, ...)
# OR a single PDB file to use directly as the one template.
TEMPLATE_PATH     = '/data/HomologyModeling_Project2/1ake_as_template/ref/1ake_A_prep_nolig.pdb'
N_TEMPLATES       = 1                      # only used when TEMPLATE_PATH is a directory
RAW_PIR           = '20_adk_A0A452CSM2.pir'   # YOUR alignment, used as is (no realign)
KNOWN             = 'ADK'
TARGET            = 'A0A452CSM2'
TARGET_FIRST_RESNUM = 1                    # first residue number of full target

ALN_PIR           = 'aln_trimmed.pir'
ALN_PAP           = 'aln_trimmed.pap'
MODELS_PER_TMPL   = 10
MD_LEVEL          = refine.slow            # try very_fast/none for more speed

RESULTS_ROOT      = 'results'              # everything goes under here now

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


def resolve_templates(template_path, n_templates):
    """Return a list of (index, pdb_filepath, atom_files_directory) tuples.

    - If template_path is a FILE: single-template mode. Returns one entry,
      index 0, pointing at that exact file. The atom_files_directory is the
      file's parent directory (so MODELLER/knowns lookups still work).
    - If template_path is a DIRECTORY: cluster mode (old behaviour). Returns
      up to n_templates entries for 0.pdb, 1.pdb, ... that actually exist.
    """
    if os.path.isfile(template_path):
        return [(0, template_path)], os.path.dirname(template_path) or '.'

    if os.path.isdir(template_path):
        found = []
        for i in range(n_templates):
            candidate = os.path.join(template_path, '%d.pdb' % i)
            if os.path.exists(candidate):
                found.append((i, candidate))
            else:
                print('Template %s not found. Skipping.' % candidate)
        return found, template_path

    raise FileNotFoundError(
        'TEMPLATE_PATH %r is neither an existing file nor a directory.'
        % template_path)


# ----------------------------------------------------------------------------
def main():
    templates, atom_dir = resolve_templates(TEMPLATE_PATH, N_TEMPLATES)
    if not templates:
        print('No templates found under %s -- nothing to do.' % TEMPLATE_PATH)
        return

    single_template_mode = len(templates) == 1

    env = Environ()
    env.io.atom_files_directory = ['.', atom_dir]

    n_trim_N, _ = trim_termini(RAW_PIR, ALN_PIR)
    offset = (TARGET_FIRST_RESNUM - 1) + n_trim_N   # add to model numbers

    aln = Alignment(env)
    aln.append(file=ALN_PIR, alignment_format='PIR')
    aln.write(file=ALN_PAP, alignment_format='PAP')
    print('Inspect %s. Output models will be renumbered by +%d so pocket '
          'residues keep full-length numbering.' % (ALN_PAP, offset))
    audit_internal_gaps(ALN_PIR)

    os.makedirs(RESULTS_ROOT, exist_ok=True)

    for i, template_file in templates:
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

        for pdb in glob.glob('%s.*.pdb' % TARGET):
            renumber_pdb(pdb, offset)

        ok = [x for x in a.outputs if x['failure'] is None]
        if ok:
            ok.sort(key=lambda x: x['Normalized DOPE score'])
            print('  z-DOPE ranking (best first):')
            for x in ok:
                print('    %-30s zDOPE=%.3f' % (x['name'],
                                                x['Normalized DOPE score']))

        # single template -> flat results/, multiple -> results/template_i/
        result_folder = RESULTS_ROOT if single_template_mode else \
            os.path.join(RESULTS_ROOT, 'template_%d' % i)
        os.makedirs(result_folder, exist_ok=True)

        for f in glob.glob('%s.*' % TARGET):
            shutil.move(f, os.path.join(result_folder, os.path.basename(f)))
        shutil.copy('ADK.pdb',
                    os.path.join(result_folder, os.path.basename(template_file)))
        os.remove('ADK.pdb')

    print('\nAll done. Results in %s/' % RESULTS_ROOT)


if __name__ == '__main__':
    main()

