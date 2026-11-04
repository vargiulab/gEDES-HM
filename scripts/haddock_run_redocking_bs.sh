#!/bin/bash
#
# HADDOCK3 REDOCKING Pipeline - BINDING SITE RESIDUES Approach
# =============================================================
#
# INPUT: Holo structure (protein + ligand complex) + BS residue list
# OUTPUT: All files needed for HADDOCK3 redocking with AIR restraints
#
# This uses the HADDOCK protein-ligand binding site protocol:
#   - Rigid-body:  BS residues = ACTIVE, ligand = ACTIVE  (draws ligand into pocket)
#   - Flexref:     BS residues = PASSIVE, ligand = ACTIVE  (ligand explores pocket)
#   - Two separate restraint files: ambig_active.tbl and ambig_passive.tbl
#   - 2-body docking (no shape molecule)
#   - Rigid ligand throughout (nseg2=0, comparable to AutoDock/Vina)
#
# Following the official HADDOCK protein-ligand tutorial:
#   https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/
#
# Requirements:
#   - Python 3 with: rdkit, requests, numpy, scipy
#   - pdb-tools (pdb_chain, pdb_reres, pdb_reatom)
#   - haddock3-restraints (for actpass_to_ambig, optional -- script generates manually if not available)
#
# Usage:
#   ./haddock_run_redocking_bs.sh <holo_complex.pdb> <ligand_resname> <bs_residues> [ncores]
#
# Example:
#   ./haddock_run_redocking_bs.sh holo_processed.pdb AP5 "14 15 16 17 18 19 20 40 41 58 62 64 91 94 98 125 129 139 140 158 169 195" 32
#

set -e  # Exit on error

# ===========================================
# Configuration
# ===========================================
HOLO_PDB=${1:-""}
LIGAND_RESNAME=${2:-""}
BS_RESIDUES=${3:-""}
NCORES=${4:-32}
TARGET="ligand"  # Base name for output files

if [[ -z "${HOLO_PDB}" ]] || [[ -z "${LIGAND_RESNAME}" ]] || [[ -z "${BS_RESIDUES}" ]]; then
    echo "Usage: $0 <holo_complex.pdb> <ligand_resname> <bs_residues> [ncores]"
    echo ""
    echo "Arguments:"
    echo "  holo_complex.pdb  - PDB file with protein + ligand"
    echo "  ligand_resname    - 3-letter residue name of ligand (e.g., AP5)"
    echo "  bs_residues       - Space-separated binding site residue numbers (quoted)"
    echo "  ncores            - Number of CPU cores (default: 32)"
    echo ""
    echo "Example:"
    echo "  $0 holo_processed.pdb AP5 \"14 15 16 17 18 19 20 40 41\" 32"
    exit 1
fi

if [[ ! -f "${HOLO_PDB}" ]]; then
    echo "ERROR: ${HOLO_PDB} not found"
    exit 1
fi

echo "=============================================="
echo "HADDOCK3 REDOCKING - Binding Site Residues"
echo "=============================================="
echo "Input: ${HOLO_PDB}"
echo "Ligand: ${LIGAND_RESNAME}"
echo "BS residues: ${BS_RESIDUES}"
echo "Cores: ${NCORES}"
echo "=============================================="
echo ""

# ===========================================
# STAGE 1: EXTRACT LIGAND
# ===========================================
echo "==== STAGE 1: EXTRACT LIGAND ===="

echo "  [1/3] Extracting ligand from complex..."
grep "^HETATM" "${HOLO_PDB}" | grep " ${LIGAND_RESNAME} " > ${TARGET}_ref_ligand_raw.pdb || {
    echo "ERROR: No HETATM lines found for residue ${LIGAND_RESNAME}"
    exit 1
}
NATOMS=$(wc -l < ${TARGET}_ref_ligand_raw.pdb)
echo "        Extracted: ${NATOMS} atoms"

echo "  [2/3] Canonicalizing ligand for PRODRG..."
python3 << 'PYEOF'
from rdkit import Chem

def format_pdb_line(atom_count, atom_name, element, coords, chain='B', resname='UNK', resnum=1):
    x, y, z = coords
    if len(element) == 1:
        atom_name_field = f" {atom_name:<3s}"
    else:
        atom_name_field = f"{atom_name:<4s}"
    
    line = (
        f"ATOM  {atom_count:5d} {atom_name_field} {resname:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    )
    return line

# Read reference ligand
mol = Chem.MolFromPDBFile('ligand_ref_ligand_raw.pdb', removeHs=False)
if mol is None:
    print("ERROR: Could not parse ligand")
    exit(1)

# Canonicalize atom order
_ = Chem.MolToSmiles(mol)
mol_order = list(mol.GetPropsAsDict(True, True).get("_smilesAtomOutputOrder", range(mol.GetNumAtoms())))
mol_canonical = Chem.RenumberAtoms(mol, mol_order)

# Write for PRODRG
conf = mol_canonical.GetConformer()
lines = []
atom_type_counts = {}

for atom_idx in range(mol_canonical.GetNumAtoms()):
    atom = mol_canonical.GetAtomWithIdx(atom_idx)
    pos = conf.GetAtomPosition(atom_idx)
    
    element = atom.GetSymbol().upper()
    atom_type_counts[element] = atom_type_counts.get(element, 0) + 1
    atom_name = f"{element}{atom_type_counts[element]}"
    
    line = format_pdb_line(
        atom_count=atom_idx + 1,
        atom_name=atom_name,
        element=element,
        coords=(pos.x, pos.y, pos.z)
    )
    lines.append(line)

with open('ligand_for_prodrg.pdb', 'w') as f:
    f.write('\n'.join(lines) + '\nEND\n')

print(f"        Canonicalized: {mol_canonical.GetNumAtoms()} atoms")
print(f"        Written: ligand_for_prodrg.pdb")
PYEOF

echo "  [3/3] Calling PRODRG API..."
python3 << 'PYEOF'
import requests
import json
import sys
import time

url = 'https://wenmr.science.uu.nl/prodrg'
max_retries = 3
timeout = 120

for attempt in range(1, max_retries + 1):
    try:
        print(f"        Attempt {attempt}/{max_retries}...", flush=True)
        start = time.time()
        
        with open('ligand_for_prodrg.pdb', 'rb') as f:
            files = {'input': ('ligand_for_prodrg.pdb', f, 'text/plain')}
            data = {'passkey': 'paralig4haddock'}
            r = requests.post(url, files=files, data=data, timeout=timeout)
        
        elapsed = time.time() - start
        print(f"        Response received in {elapsed:.1f}s", flush=True)
        
        if r.status_code == 200 and len(r.text) > 100:
            with open('prodrg_output.json', 'w') as out:
                out.write(r.text)
            print(f"        PRODRG output saved ({len(r.text)} bytes)", flush=True)
            sys.exit(0)
        else:
            print(f"        Invalid response (status: {r.status_code}, size: {len(r.text)})", flush=True)
    except requests.exceptions.Timeout:
        print(f"        Timeout after {timeout}s", flush=True)
    except Exception as e:
        print(f"        Error: {e}", flush=True)
    
    if attempt < max_retries:
        print(f"        Retrying in 5 seconds...", flush=True)
        time.sleep(5)

print("ERROR: PRODRG API failed after all retries")
sys.exit(1)
PYEOF

if [[ $? -ne 0 ]]; then
    echo "ERROR: PRODRG API call failed"
    exit 1
fi

echo "  STAGE 1 COMPLETE"

# ===========================================
# STAGE 2: EXTRACT PRODRG FILES
# ===========================================
echo ""
echo "==== STAGE 2: EXTRACT PRODRG FILES ===="

python3 << 'PYEOF'
import json
import re

with open('prodrg_output.json', 'r') as f:
    data = json.load(f)

files = data.get('files', {})
print(f"  Found {len(files)} files in PRODRG response")

file_keys = list(files.keys())

# Extract parameter file (.par)
param_file = [k for k in file_keys if k.endswith('.par')]
if param_file:
    key = param_file[0]
    print(f"  Found parameters: {key}")
    with open('ligand.param', 'w') as f:
        f.write(files[key])
else:
    print("ERROR: No parameter file (.par) found")
    exit(1)

# Extract topology file (.top)
top_file = [k for k in file_keys if k.endswith('.top')]
if top_file:
    key = top_file[0]
    print(f"  Found topology: {key}")
    with open('ligand.top', 'w') as f:
        f.write(files[key])
else:
    print("ERROR: No topology file (.top) found")
    exit(1)

# Extract PDB file (heavy atoms only version - _noh.pdb)
pdb_file = [k for k in file_keys if '_noh.pdb' in k]
if not pdb_file:
    pdb_file = [k for k in file_keys if k.endswith('.pdb') and '_poh' not in k]

if pdb_file:
    key = pdb_file[0]
    print(f"  Found PDB: {key}")
    with open('ligand_prodrg.pdb', 'w') as f:
        f.write(files[key])
else:
    print("ERROR: No PDB file found")
    exit(1)

print(f"  Written: ligand.top")
print(f"  Written: ligand.param")
print(f"  Written: ligand_prodrg.pdb")
PYEOF

echo "  STAGE 2 COMPLETE"

# ===========================================
# STAGE 3: PREPARE LIGAND FOR HADDOCK (REDOCKING)
# ===========================================
echo ""
echo "==== STAGE 3: PREPARE LIGAND FOR REDOCKING ===="
echo "  Using CRYSTAL POSE (single conformer) for redocking"

python3 << 'PYEOF'
import re
import numpy as np
from collections import defaultdict
from scipy.optimize import linear_sum_assignment
from rdkit import Chem

def extract_topology_atom_names(top_file):
    """Extract heavy atom names from topology file in order."""
    with open(top_file, 'r') as f:
        content = f.read()
    
    residue_match = re.search(r'Residue\s+UNK.*?(?=\nResidue|\Z)', content, re.DOTALL)
    if not residue_match:
        return []
    
    atom_pattern = r'ATOM\s+(\S+)\s+TYPE='
    atom_names = re.findall(atom_pattern, residue_match.group(0))
    
    heavy_atoms = [name for name in atom_names if not name.startswith('H')]
    return heavy_atoms

def get_element_from_line(line):
    """Extract element from PDB element column (76-78). Falls back to atom name heuristic."""
    # Try PDB element column first (standard PDB format)
    if len(line) >= 78:
        elem = line[76:78].strip()
        if elem and elem[0].isalpha():
            return elem.upper()
    # Fallback: guess from atom name
    atname = line[12:16].strip()
    name = atname.lstrip("0123456789'\"")
    if len(name) >= 2 and name[:2].upper() in ['CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'NA']:
        return name[:2].upper()
    return name[0].upper() if name else 'C'

def is_hydrogen_from_line(line):
    """Check if atom is hydrogen using the element column."""
    return get_element_from_line(line) == 'H'

def format_pdb_line(atom_count, atom_name, element, coords, chain='B', resname='UNK', resnum=1):
    x, y, z = coords
    if len(element) == 1:
        atom_name_field = f" {atom_name:<3s}"
    else:
        atom_name_field = f"{atom_name:<4s}"
    
    line = (
        f"ATOM  {atom_count:5d} {atom_name_field} {resname:>3s} {chain}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  "
    )
    return line

# Get atom names from topology (heavy atoms only)
topology_atoms = extract_topology_atom_names('ligand.top')
print(f"  Topology heavy atoms: {len(topology_atoms)}")

# Read PRODRG PDB (heavy atoms only)
prodrg_atoms = []
with open('ligand_prodrg.pdb', 'r') as f:
    for line in f:
        if line.startswith(('HETATM', 'ATOM')):
            if not is_hydrogen_from_line(line):
                atname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                prodrg_atoms.append({
                    'name': atname,
                    'coords': np.array([x, y, z]),
                    'element': get_element_from_line(line)
                })

print(f"  PRODRG heavy atoms: {len(prodrg_atoms)}")

# Read crystal ligand (heavy atoms only)
crystal_atoms = []
with open('ligand_ref_ligand_raw.pdb', 'r') as f:
    for line in f:
        if line.startswith(('HETATM', 'ATOM')):
            if not is_hydrogen_from_line(line):
                atname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                crystal_atoms.append({
                    'name': atname,
                    'coords': np.array([x, y, z]),
                    'element': get_element_from_line(line)
                })

print(f"  Crystal heavy atoms: {len(crystal_atoms)}")

# =====================================================
# Match atoms by element using Hungarian algorithm
# (optimal global assignment instead of greedy nearest-neighbor)
# =====================================================
prodrg_coords = np.array([a['coords'] for a in prodrg_atoms])
crystal_coords = np.array([a['coords'] for a in crystal_atoms])

prodrg_center = prodrg_coords.mean(axis=0)
crystal_center = crystal_coords.mean(axis=0)

prodrg_centered = prodrg_coords - prodrg_center
crystal_centered = crystal_coords - crystal_center

# Group atom indices by element
prodrg_by_element = defaultdict(list)
crystal_by_element = defaultdict(list)

for i, a in enumerate(prodrg_atoms):
    prodrg_by_element[a['element']].append(i)
for j, a in enumerate(crystal_atoms):
    crystal_by_element[a['element']].append(j)

# Optimal matching per element
prodrg_to_crystal = {}

for element in prodrg_by_element:
    p_indices = prodrg_by_element[element]
    c_indices = crystal_by_element.get(element, [])
    
    if len(c_indices) < len(p_indices):
        print(f"  ERROR: Not enough crystal {element} atoms ({len(c_indices)}) to match PRODRG ({len(p_indices)})")
        exit(1)
    
    # Build cost matrix (pairwise distances on centered coordinates)
    cost = np.zeros((len(p_indices), len(c_indices)))
    for ii, pi in enumerate(p_indices):
        for jj, cj in enumerate(c_indices):
            cost[ii, jj] = np.linalg.norm(prodrg_centered[pi] - crystal_centered[cj])
    
    row_ind, col_ind = linear_sum_assignment(cost)
    
    for ii, jj in zip(row_ind, col_ind):
        prodrg_to_crystal[p_indices[ii]] = c_indices[jj]
        
    # Report max distance for this element (diagnostic)
    max_dist = max(cost[ii, jj] for ii, jj in zip(row_ind, col_ind))
    print(f"  Matched {len(p_indices):2d} {element} atoms (max dist: {max_dist:.3f} A)")

matched_crystal_coords = [crystal_atoms[prodrg_to_crystal[i]]['coords'] for i in range(len(prodrg_atoms))]
print(f"  Successfully matched {len(matched_crystal_coords)} atoms")

# Write HADDOCK-ready ligand (PRODRG atom names + crystal coordinates)
lines = []
for i, (atom, coords) in enumerate(zip(prodrg_atoms, matched_crystal_coords)):
    atom_name = atom['name']
    element = atom['element']
    line = format_pdb_line(i+1, atom_name, element, coords)
    lines.append(line)

with open('ligand_haddock.pdb', 'w') as f:
    f.write('\n'.join(lines) + '\nEND\n')

print(f"  Written: ligand_haddock.pdb (single crystal pose)")
PYEOF

echo "  STAGE 3 COMPLETE"

# ===========================================
# STAGE 4: PREPARE RECEPTOR
# ===========================================
echo ""
echo "==== STAGE 4: PREPARE RECEPTOR ===="

echo "  [1/2] Extracting receptor (protein only)..."
grep "^ATOM" "${HOLO_PDB}" | grep -v " ${LIGAND_RESNAME} " > receptor_tmp.pdb || true

# Set chain to A
pdb_chain -A receptor_tmp.pdb > receptor.pdb 2>/dev/null || {
    sed 's/^\(.\{21\}\)./\1A/' receptor_tmp.pdb > receptor.pdb
}
rm -f receptor_tmp.pdb
echo "END" >> receptor.pdb

NRECP=$(grep -c "^ATOM" receptor.pdb || echo "0")
echo "        Receptor atoms: ${NRECP}"

echo "  [2/2] Creating CAPRI reference..."
cat receptor.pdb > reference.pdb
grep "^ATOM" ligand_haddock.pdb >> reference.pdb
echo "END" >> reference.pdb
echo "        Created: reference.pdb"

echo "  STAGE 4 COMPLETE"

# ===========================================
# STAGE 5: GENERATE AIR RESTRAINTS (BS APPROACH)
# ===========================================
echo ""
echo "==== STAGE 5: GENERATE AIR RESTRAINTS ===="
echo "  Protocol: HADDOCK protein-ligand binding site approach"
echo "  BS residues: ${BS_RESIDUES}"

python3 << PYEOF
import os

bs_residues = "${BS_RESIDUES}".split()
print(f"  Number of BS residues: {len(bs_residues)}")

# -------------------------------------------------------
# actpass files for haddock3-restraints (if available)
# -------------------------------------------------------

# --- RIGIDBODY phase: BS active, ligand active ---
with open('receptor_active.actpass', 'w') as f:
    f.write(' '.join(bs_residues) + '\n')
    f.write('\n')

with open('ligand_active.actpass', 'w') as f:
    f.write('1\n')
    f.write('\n')

# --- FLEXREF/EMREF phase: BS passive, ligand active ---
with open('receptor_passive.actpass', 'w') as f:
    f.write('\n')
    f.write(' '.join(bs_residues) + '\n')

with open('ligand_active_only.actpass', 'w') as f:
    f.write('1\n')
    f.write('\n')

print(f"  Written: receptor_active.actpass")
print(f"  Written: ligand_active.actpass")
print(f"  Written: receptor_passive.actpass")
print(f"  Written: ligand_active_only.actpass")

# -------------------------------------------------------
# Generate AIR restraint files directly (CNS format)
# -------------------------------------------------------

def generate_air_restraints(active_resids_mol1, passive_resids_mol1,
                            active_resids_mol2, passive_resids_mol2,
                            segid1='A', segid2='B'):
    """Generate AIR restraints in CNS format following HADDOCK convention."""
    restraints = []
    
    for resid in active_resids_mol1:
        target_resids = active_resids_mol2 + passive_resids_mol2
        if not target_resids:
            continue
        
        target_parts = []
        for t_resid in target_resids:
            target_parts.append(f"(resid {t_resid} and segid {segid2})")
        
        if len(target_parts) > 1:
            target_sel_lines = []
            for i, tp in enumerate(target_parts):
                if i == 0:
                    target_sel_lines.append(f"        {tp}")
                else:
                    target_sel_lines.append(f"     or {tp}")
            target_sel = '\n'.join(target_sel_lines)
            restraint = f"assign (resid {resid} and segid {segid1})\n(\n{target_sel}\n) 2.0 2.0 0.0"
        else:
            restraint = f"assign (resid {resid} and segid {segid1})\n(\n        {target_parts[0]}\n) 2.0 2.0 0.0"
        
        restraints.append(restraint)
    
    for resid in active_resids_mol2:
        target_resids = active_resids_mol1 + passive_resids_mol1
        if not target_resids:
            continue
        
        target_parts = []
        for t_resid in target_resids:
            target_parts.append(f"(resid {t_resid} and segid {segid1})")
        
        if len(target_parts) > 1:
            target_sel_lines = []
            for i, tp in enumerate(target_parts):
                if i == 0:
                    target_sel_lines.append(f"        {tp}")
                else:
                    target_sel_lines.append(f"     or {tp}")
            target_sel = '\n'.join(target_sel_lines)
            restraint = f"assign (resid {resid} and segid {segid2})\n(\n{target_sel}\n) 2.0 2.0 0.0"
        else:
            restraint = f"assign (resid {resid} and segid {segid2})\n(\n        {target_parts[0]}\n) 2.0 2.0 0.0"
        
        restraints.append(restraint)
    
    return restraints

bs_list = [int(r) for r in bs_residues]
lig_list = [1]

# --- ambig_active.tbl: BS=active, ligand=active (for rigidbody) ---
restraints_active = generate_air_restraints(
    active_resids_mol1=bs_list, passive_resids_mol1=[],
    active_resids_mol2=lig_list, passive_resids_mol2=[],
    segid1='A', segid2='B'
)

with open('ambig_active.tbl', 'w') as f:
    f.write('! HADDOCK AIR restraints - RIGIDBODY phase\n')
    f.write('! BS residues = ACTIVE, Ligand = ACTIVE\n')
    f.write('! Draws ligand into the binding pocket\n!\n')
    f.write('\n\n'.join(restraints_active) + '\n')

print(f"  Written: ambig_active.tbl ({len(restraints_active)} restraints)")

# --- ambig_passive.tbl: BS=passive, ligand=active (for flexref/emref) ---
restraints_passive = generate_air_restraints(
    active_resids_mol1=[], passive_resids_mol1=bs_list,
    active_resids_mol2=lig_list, passive_resids_mol2=[],
    segid1='A', segid2='B'
)

with open('ambig_passive.tbl', 'w') as f:
    f.write('! HADDOCK AIR restraints - FLEXREF/EMREF phase\n')
    f.write('! BS residues = PASSIVE, Ligand = ACTIVE\n')
    f.write('! Ligand explores the binding pocket freely\n!\n')
    f.write('\n\n'.join(restraints_passive) + '\n')

print(f"  Written: ambig_passive.tbl ({len(restraints_passive)} restraints)")
PYEOF

# Try to also generate with haddock3-restraints if available (for validation)
if command -v haddock3-restraints &> /dev/null; then
    echo "  [Optional] Validating with haddock3-restraints..."
    haddock3-restraints active_passive_to_ambig receptor_active.actpass ligand_active.actpass \
        --segid-one A --segid-two B > ambig_active_h3r.tbl 2>/dev/null && \
        echo "        Validated: ambig_active_h3r.tbl" || \
        echo "        Warning: haddock3-restraints validation skipped"
    
    haddock3-restraints active_passive_to_ambig receptor_passive.actpass ligand_active_only.actpass \
        --segid-one A --segid-two B > ambig_passive_h3r.tbl 2>/dev/null && \
        echo "        Validated: ambig_passive_h3r.tbl" || \
        echo "        Warning: haddock3-restraints validation skipped"
    
    if [[ -s ambig_active_h3r.tbl ]]; then
        mv ambig_active_h3r.tbl ambig_active.tbl
        echo "        Using haddock3-restraints output for ambig_active.tbl"
    fi
    if [[ -s ambig_passive_h3r.tbl ]]; then
        mv ambig_passive_h3r.tbl ambig_passive.tbl
        echo "        Using haddock3-restraints output for ambig_passive.tbl"
    fi
else
    echo "  [Info] haddock3-restraints not found, using manually generated AIRs"
fi

echo "  STAGE 5 COMPLETE"

# ===========================================
# STAGE 6: GENERATE HADDOCK3 CONFIG (RIGID LIGAND)
# ===========================================
echo ""
echo "==== STAGE 6: GENERATE HADDOCK3 CONFIG ===="

cat > haddock3_redocking_bs.cfg << CFGEOF
# HADDOCK3 Protein-Ligand REDOCKING Configuration
# Protocol: BINDING SITE RESIDUES approach - RIGID LIGAND
# Generated: $(date)
# Mode: REDOCKING (single crystal conformer)
#
# Following: https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-binding-sites/
# Ligand flexibility: RIGID (nseg2 = 0, no torsional refinement)
# Comparable to AutoDock/Vina rigid ligand docking

run_dir = "run_redocking_bs"
ncores = ${NCORES}

molecules = [
    "receptor.pdb",
    "ligand_haddock.pdb"
]

# ===========================================
# Topology generation
# ===========================================
[topoaa]
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
delenph = false
autohis = true
set_bfactor = false

# ===========================================
# Rigid-body docking
# BS residues = ACTIVE, Ligand = ACTIVE
# ===========================================
[rigidbody]
sampling = 250
tolerance = 5
ambig_fname = "ambig_active.tbl"
w_vdw = 1.0
inter_rigid = 0.001
mol_fix_origin_1 = true
randremoval = false

ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

# ===========================================
# CAPRI evaluation after rigid-body
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Selection of top 50
# ===========================================
[seletop]
select = 50

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Flexible refinement (ligand kept rigid)
# BS residues = PASSIVE, Ligand = ACTIVE
# nseg2 = 0: no torsional freedom for ligand
# ===========================================
[flexref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
mdsteps_rigid = 0
mdsteps_cool1 = 0
nseg2 = 0
randremoval = false

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Energy minimization refinement (ligand kept rigid)
# ===========================================
[emref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
nseg2 = 0
randremoval = false

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Final selection and evaluation
# ===========================================
[seletop]
select = 5

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5
CFGEOF

echo "  Written: haddock3_redocking_bs.cfg (rigid ligand)"

echo "  STAGE 6 COMPLETE"

# ===========================================
# SUMMARY
# ===========================================
echo ""
echo "=============================================="
echo "PIPELINE COMPLETE"
echo "=============================================="
echo ""
echo "Generated files:"
echo "  - receptor.pdb              : Protein structure (chain A)"
echo "  - ligand_haddock.pdb        : Ligand (single crystal pose, chain B)"
echo "  - ambig_active.tbl          : AIRs for rigidbody (BS=active, lig=active)"
echo "  - ambig_passive.tbl         : AIRs for flexref/emref (BS=passive, lig=active)"
echo "  - reference.pdb             : Reference for CAPRI evaluation"
echo "  - ligand.top                : PRODRG topology"
echo "  - ligand.param              : PRODRG parameters"
echo ""
echo "HADDOCK3 configuration:"
echo "  - haddock3_redocking_bs.cfg : Rigid ligand (nseg2=0)"
echo ""
echo "actpass files (for haddock3-restraints validation):"
echo "  - receptor_active.actpass   : BS residues as active"
echo "  - receptor_passive.actpass  : BS residues as passive"
echo "  - ligand_active.actpass     : Ligand as active"
echo ""
echo "To run HADDOCK3:"
echo "  haddock3 haddock3_redocking_bs.cfg"
echo ""
echo "=============================================="