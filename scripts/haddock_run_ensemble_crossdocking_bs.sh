#!/bin/bash
#
# HADDOCK3 ENSEMBLE CROSS-DOCKING: All Homology Models in ONE run
# Protocol: BINDING SITE RESIDUES approach - RIGID LIGAND
# ================================================================
#
# Difference vs. the per-model (one-by-one) script:
#   - All receptor homology models (same sequence/length, different
#     conformations, selected by DOPE) are concatenated into a single
#     multi-MODEL ensemble PDB and docked together as molecule 1.
#   - HADDOCK distributes sampling across the ensemble members.
#   - Pipeline: rigidbody -> seletop (by score) -> flexref -> emref -> caprieval.
#   - emref (energy-minimization refinement) relieves the w_vdw=0 clashes that
#     flexref (mdsteps=0) leaves, giving meaningful (negative) HADDOCK scores.
#   - NO HADDOCK clustering module is run (neither iL-RMSD nor FCC): emref is
#     the final stage. Poses are kept individually, their true RMSD is computed
#     afterwards with cpptraj (symmetry-corrected), and pose clustering is done
#     with the Amber/cpptraj workflow - this is what is reported in the paper.
#
# Everything else (restraints, cmrest/cmtight, w_vdw=0, tolerance=20,
# semi-flexible BS segments, rigid ligand) is identical to the
# one-by-one script.
#
# The [caprieval] blocks after rigidbody / seletop / flexref are analysis
# checkpoints only (they neither filter nor re-rank the models); the reported
# evaluation is the last one, after emref.
#
# Resulting step folders (zero-padded in recent HADDOCK3, e.g. 00_topoaa):
#   0_topoaa 1_rigidbody 2_caprieval 3_seletop 4_caprieval
#   5_flexref 6_caprieval 7_emref 8_caprieval   <- last one = reported evaluation
#
# Usage:
#   ./haddock_run_ensemble_crossdocking_bs.sh <models_dir> <redocking_dir> [ncores]
#
# Example:
#   ./haddock_run_ensemble_crossdocking_bs.sh selected_best_structures redocking 32
#

set -e

MODELS_DIR=${1:-"selected_best_structures"}
REDOCKING_DIR=${2:-"redocking"}
NCORES=${3:-32}
OUTDIR="ensemble_crossdocking_bs"

if [[ ! -d "${MODELS_DIR}" ]]; then
    echo "ERROR: Models directory ${MODELS_DIR} not found"
    exit 1
fi

# Check required files from redocking
for f in ligand_haddock.pdb ligand.top ligand.param ambig_active.tbl reference.pdb; do
    if [[ ! -f "${REDOCKING_DIR}/${f}" ]]; then
        echo "ERROR: ${REDOCKING_DIR}/${f} not found"
        echo "       Run redocking setup first!"
        exit 1
    fi
done

# Count models and derive scaled sampling
NMODELS=$(ls "${MODELS_DIR}"/*.pdb 2>/dev/null | wc -l | tr -d ' ')
if [[ "${NMODELS}" -eq 0 ]]; then
    echo "ERROR: No .pdb models found in ${MODELS_DIR}"
    exit 1
fi
# Rigidbody sampling = 250 per conformer, exactly like the one-by-one runs.
# No iL-RMSD matrix anymore -> no 10k cap -> full sampling is fine.
SAMPLING=$((250 * NMODELS))

# --- Stage-1 (AFTER rigidbody): SCORE-based funnel to flexref (no iL-RMSD) ---
# Keep 50 per conformer (50*N), matching one-by-one's 250 rigid -> 50 flex.
# NOTE: this is a GLOBAL top-(50*N) by score, not a guaranteed 50-per-conformer
# (HADDOCK ensembles can't do per-conformer seletop) -- validate with cpptraj.
SELETOP_FLEX=$((50 * NMODELS))   # = 5050 for N=101 -> flexref

# --- Stage-2 (AFTER emref): none. No HADDOCK clustering module is run. ---
# Every emref pose is kept; clustering/ranking is done downstream with
# cpptraj (symmetry-corrected RMSD), as reported in the paper.

echo "=============================================="
echo "HADDOCK3 ENSEMBLE CROSS-DOCKING: BS Residues"
echo "Protocol: RIGID LIGAND (single ensemble run)"
echo "=============================================="
echo "Models dir:   ${MODELS_DIR} (N=${NMODELS} conformers)"
echo "Redocking:    ${REDOCKING_DIR}"
echo "Output:       ${OUTDIR}"
echo "Cores:        ${NCORES}"
echo "----------------------------------------------"
echo "  Rigidbody:   sampling = ${SAMPLING}"
echo "  w_vdw=0.0, tolerance=20, cmrest=true, cmtight=true"
echo "  Stage-1 (post-rigidbody): seletop = ${SELETOP_FLEX} (by HADDOCK score) -> flexref"
echo "  Refinement: emref (EM) after flexref -> meaningful (negative) scores"
echo "  Stage-2 (post-emref):     none (no HADDOCK clustering module)"
echo "  NOTE: all emref poses kept; clustering + true RMSD via cpptraj (symmetry-corrected)"
echo "=============================================="
echo ""

mkdir -p "${OUTDIR}"

# ===========================================
# Extract BS residues from ambig_active.tbl
# (identical to one-by-one script)
# ===========================================
echo "==== EXTRACTING BS RESIDUES FROM RESTRAINTS ===="
# (portable: GNU-only 'grep -oP ... \K' replaced by two BSD/GNU-safe greps)
BS_RESIDUES=$(grep -o -i 'resi[[:space:]]*[0-9][0-9]*' "${REDOCKING_DIR}/ambig_active.tbl" | \
              grep -o '[0-9][0-9]*' | grep -v '^1$' | sort -un | tr '\n' ' ')
echo "  BS residues: ${BS_RESIDUES}"
if [[ -z "${BS_RESIDUES}" ]]; then
    echo "ERROR: No BS residues extracted from ambig_active.tbl"
    exit 1
fi

# ===========================================
# Generate semi-flexible segment definitions
# (identical to one-by-one script)
# ===========================================
echo "==== GENERATING SEMI-FLEXIBLE SEGMENTS ===="
SEMIFLEX_BLOCK=$(python3 << PYEOF
bs_residues = [${BS_RESIDUES// /, }]

ranges = []
for r in sorted(set(bs_residues)):
    ranges.append((r - 1, r + 1))

merged = [ranges[0]]
for start, end in ranges[1:]:
    prev_start, prev_end = merged[-1]
    if start <= prev_end + 1:
        merged[-1] = (prev_start, max(prev_end, end))
    else:
        merged.append((start, end))

lines = []
lines.append("# Semi-flexible segments for receptor (molecule 1)")
lines.append(f"# {len(merged)} merged segments from {len(bs_residues)} BS residues")
for i, (start, end) in enumerate(merged, start=1):
    contained = [r for r in sorted(set(bs_residues)) if start <= r <= end]
    lines.append(f"# Segment {i}: res {','.join(map(str, contained))}")
    lines.append(f"seg_sta_1_{i} = {start}")
    lines.append(f"seg_end_1_{i} = {end}")

print('\n'.join(lines))
PYEOF
)
echo "${SEMIFLEX_BLOCK}"
echo ""

# ===========================================
# BUILD THE RECEPTOR ENSEMBLE
# Each model -> ATOM-only, chain A; then all
# models concatenated into one multi-MODEL PDB.
# ===========================================
echo "==== BUILDING RECEPTOR ENSEMBLE ===="

WORKDIR="${OUTDIR}"
RECEPTOR_PARTS=()
idx=0
for MODEL_PDB in "${MODELS_DIR}"/*.pdb; do
    idx=$((idx + 1))
    PART="${WORKDIR}/receptor_$(printf '%03d' ${idx}).pdb"

    grep "^ATOM" "${MODEL_PDB}" > "${WORKDIR}/_tmp.pdb" || true
    if command -v pdb_chain >/dev/null 2>&1; then
        pdb_chain -A "${WORKDIR}/_tmp.pdb" > "${PART}" 2>/dev/null || \
            sed 's/^\(.\{21\}\)./\1A/' "${WORKDIR}/_tmp.pdb" > "${PART}"
    else
        sed 's/^\(.\{21\}\)./\1A/' "${WORKDIR}/_tmp.pdb" > "${PART}"
    fi
    rm -f "${WORKDIR}/_tmp.pdb"
    RECEPTOR_PARTS+=("${PART}")
done

ENSEMBLE="${WORKDIR}/receptor_ensemble.pdb"
if command -v pdb_mkensemble >/dev/null 2>&1; then
    pdb_mkensemble "${RECEPTOR_PARTS[@]}" > "${ENSEMBLE}"
else
    # Fallback: wrap each model in MODEL/ENDMDL, drop END/TER between
    echo "  (pdb_mkensemble not found -- using python fallback)"
    python3 - "${ENSEMBLE}" "${RECEPTOR_PARTS[@]}" << 'PYEOF'
import sys
out_path = sys.argv[1]
parts = sys.argv[2:]
with open(out_path, "w") as out:
    for i, p in enumerate(parts, start=1):
        out.write(f"MODEL     {i:>4}\n")
        with open(p) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM", "ANISOU", "TER")):
                    out.write(line)
        out.write("ENDMDL\n")
    out.write("END\n")
PYEOF
fi

NENS=$(grep -c '^MODEL' "${ENSEMBLE}" || true)
echo "  Ensemble written: ${ENSEMBLE} (${NENS} models)"
if [[ "${NENS}" -ne "${NMODELS}" ]]; then
    echo "WARNING: ensemble model count (${NENS}) != input count (${NMODELS})"
fi

# ===========================================
# Copy redocking inputs into the run dir
# ===========================================
cp "${REDOCKING_DIR}/ligand_haddock.pdb" "${WORKDIR}/"
cp "${REDOCKING_DIR}/ligand.top"         "${WORKDIR}/"
cp "${REDOCKING_DIR}/ligand.param"       "${WORKDIR}/"
cp "${REDOCKING_DIR}/ambig_active.tbl"   "${WORKDIR}/"
cp "${REDOCKING_DIR}/reference.pdb"      "${WORKDIR}/"

# ===========================================
# WRITE THE SINGLE ENSEMBLE CONFIG
# ===========================================
cat > "${WORKDIR}/haddock3_ensemble_bs_rigid.cfg" << CFGEOF
# HADDOCK3 ENSEMBLE Cross-docking: ${NMODELS} homology models
# Protocol: Binding Site Residues approach - RIGID LIGAND
run_dir = "run_ensemble_crossdock_bs"
ncores = ${NCORES}
mode = "local"

molecules = [
    "receptor_ensemble.pdb",
    "ligand_haddock.pdb"
]

# ===========================================
# Topology generation
# (generated for every ensemble member)
# ===========================================
[topoaa]
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
delenph = false
autohis = true
set_bfactor = false

# ===========================================
# Rigid-body docking
# sampling scaled by N so each conformer is
# sampled ~250 times.
# ===========================================
[rigidbody]
sampling = ${SAMPLING}
tolerance = 20
ambig_fname = "ambig_active.tbl"
w_vdw = 0.0
inter_rigid = 0.001
mol_fix_origin_1 = true
cmrest = true
cmtight = true
randremoval = false
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# STAGE-1 (after rigidbody): SCORE-based funnel
# No iL-RMSD matrix here -> no 10k-model cap.
# Keep the top ${SELETOP_FLEX} rigidbody models by
# HADDOCK score and pass them to flexref.
# ===========================================
[seletop]
select = ${SELETOP_FLEX}

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Flexible refinement (rigid ligand, nseg2=0)
# semi-flexible BS side-chains only
# ===========================================
[flexref]
tolerance = 50
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"
mdsteps_rigid = 0
mdsteps_cool1 = 0
mdsteps_cool3 = 0
nseg2 = 0
randremoval = false
${SEMIFLEX_BLOCK}

[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5

# ===========================================
# Energy-minimization refinement (emref)
# Relieves the w_vdw=0 pocket-penetration clashes that
# flexref (mdsteps=0) leaves behind, so HADDOCK scores
# become meaningful (negative for good complexes).
# NOTE: EM cannot cross energy barriers -> severely
# clashed (non-native, high-RMSD) poses may stay
# positive; switch to [mdref] if you must relax them.
# Confirm params with 'haddock3-cfg -m emref'.
# ===========================================
[emref]
tolerance = 5
ambig_fname = "ambig_active.tbl"
ligand_param_fname = "ligand.param"
ligand_top_fname = "ligand.top"

# ===========================================
# FINAL CAPRI evaluation of the emref models
# (i-RMSD, l-RMSD, Fnat, DockQ vs. reference).
# emref is the last refinement stage.
# Pose clustering is done downstream with Amber/cpptraj
# on symmetry-corrected RMSD -- that is the protocol reported in the paper.
# ===========================================
[caprieval]
reference_fname = "reference.pdb"
irmsd_cutoff = 3.5
CFGEOF

echo ""
echo "  Config written: ${WORKDIR}/haddock3_ensemble_bs_rigid.cfg"

# ===========================================
# RUN (single ensemble run)
# ===========================================
echo ""
echo "==== RUNNING HADDOCK3 (ENSEMBLE) ===="
STARTTIME=$(date +%s)
if (cd "${WORKDIR}" && haddock3 haddock3_ensemble_bs_rigid.cfg > haddock3_ensemble.log 2>&1); then
    ELAPSED=$(( $(date +%s) - STARTTIME ))
    echo "  DONE (${ELAPSED}s)"
else
    ELAPSED=$(( $(date +%s) - STARTTIME ))
    echo "  FAILED (${ELAPSED}s) -- check ${WORKDIR}/haddock3_ensemble.log"
    exit 1
fi

echo ""
echo "=============================================="
echo "ENSEMBLE RUN COMPLETE"
echo "  Results: ${WORKDIR}/run_ensemble_crossdock_bs"
echo "  Ranked poses:  .../*_caprieval/capri_ss.tsv (last one = post-emref)"
echo "  Final models:  .../*_emref/ (no clustering; capri cluster columns are '-')"
echo "=============================================="
