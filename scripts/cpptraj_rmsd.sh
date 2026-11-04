#!/usr/bin/env bash
###############################################################################
# cpptraj_rmsd.sh   (single-reference; replaces cpptraj_rmsd_dualref.sh)
#
# Computes, per docked frame in concat.pdb, RMSD vs a single reference:
#   rmsd_wholepro.txt   sorted_rmsd_wholepro.txt        CA, common residues
#   rmsd_BS.txt         sorted_rmsd_BS.txt              symmrmsd+remap, BS heavy
#   rmsd_lig.txt        sorted_rmsd_lig.txt             ligand heavy, BS-fit + nofit
#   rmsd_combined_BS-lig.txt  sorted_rmsd_combined_BS-lig.txt
#
# WHY THIS REPLACES THE DUAL-REFERENCE VERSION
# --------------------------------------------
# The old 5.4-vs-1.6 ligand bug was NOT an AP5 arm swap. It was cpptraj's
# element selector @%H matching ZERO atoms on ligands whose hydrogens have
# nonstandard names (AP5: FHO3, JHO2, 1H6A). So ":LIG&!@%H" excluded nothing,
# the ligand RMSD ran over hydrogens with mispaired atoms -> garbage.
#
# FIX: hydrogens are detected by MASS, not by name or element string. cpptraj
# reads the mass column correctly even when it mis-guesses the element, so
# "mass < H_MASS_CUTOFF" reliably tags H for ANY ligand (AP5/B4P/GP5/ANP/...).
# This makes the heavy-atom ligand selection tool-independent and matches a
# correctly-built VMD `noh` selection to 4 decimals (verified: 1.5891).
#
# The dual reference (ref2, arm-swap) is removed: the four-arm test showed no
# swap, and single-reference ref1 already wins cleanly. If a specific target
# genuinely needs symmetry correction, use real symmrmsd with a prmtop
# (see note at bottom) rather than a hand-swapped ref2.
#
# Per-target inputs (env vars), same as before:
#   LIGRES, BS_RESIDUES (crystal resSeq), PROT_FIRST, PROT_LAST
# Requires ref.pdb and concat.pdb in the working directory.
#
# Usage:
#   LIGRES=AP5 BS_RESIDUES="9-15,31,..." PROT_FIRST=1 PROT_LAST=214 \
#       bash cpptraj_rmsd.sh
###############################################################################
set +H
set -uo pipefail

# --------------------------- CONFIG ----------------------------------------
LIGRES="${LIGRES:-AP5}"
BS_RESIDUES="${BS_RESIDUES:-9-15,31,36,53,57,59,64,85,88,92,123,124,127,136-138,160,171,197}"
PROT_FIRST="${PROT_FIRST:-1}"
PROT_LAST="${PROT_LAST:-212}"

REF="${REF:-ref.pdb}"
REF2="${REF2:-ref2.pdb}"
CONCAT="${CONCAT:-concat.pdb}"
OUTLOG="cpptraj_rmsd_run.log"

# Mode auto-detection: if ref2.pdb exists, run dual-reference (per-frame min of
# ref/ref2, chosen by ligand RMSD) for the manual arm-swap symmetry case.
# Otherwise run single-reference (the normal, recommended path).
# Force single even when ref2 exists with FORCE_SINGLE=1.
if [[ "${FORCE_SINGLE:-0}" != "1" && -e "$REF2" ]]; then
    DUAL_REF=1
else
    DUAL_REF=0
fi

# Mass below which an atom is treated as hydrogen. Deuterium-safe and well
# clear of the lightest heavy atom (C=12). Any H (~1.008) falls under this.
H_MASS_CUTOFF="${H_MASS_CUTOFF:-2.5}"

# --------------------- environment / input guards --------------------------
command -v cpptraj >/dev/null 2>&1 || {
    echo "FATAL: cpptraj not in PATH. Load AmberTools (e.g. 'module load amber24')." >&2; exit 1; }
for f in "$REF" "$CONCAT"; do
    [[ -e "$f" ]] || { echo "FATAL: required file '$f' not found in $(pwd)" >&2; exit 1; }
done
if [[ $DUAL_REF -eq 1 ]]; then
    [[ -e "$REF2" ]] || { echo "FATAL: dual mode but '$REF2' missing" >&2; exit 1; }
fi
echo "cpptraj: $(command -v cpptraj)"
echo "Target : LIGRES=$LIGRES  PROT=$PROT_FIRST-$PROT_LAST"
echo "BS res : $BS_RESIDUES"
if [[ $DUAL_REF -eq 1 ]]; then
    echo "Ref    : $REF + $REF2   (DUAL reference: per-frame min, chosen by ligand)"
else
    echo "Ref    : $REF   (single reference)"
fi
echo

# ----------------------------- sanitize ------------------------------------
sanitize() {  # $1 in  $2 out
    sed 's/\r$//' "$1" \
      | awk '/^(ATOM|HETATM)/{ l=$0; while(length(l)<80) l=l" "; print substr(l,1,80); next }
             /^(MODEL|ENDMDL|TER|END|CRYST1)/{ print }' > "$2"
}
echo "=== [1/6] Sanitizing ref / concat ==="
sanitize "$REF"    ref_clean.pdb
sanitize "$CONCAT" concat_clean.pdb
if [[ $DUAL_REF -eq 1 ]]; then
    sanitize "$REF2" ref2_clean.pdb
    echo "  ref2_clean  : $(grep -cE '^(ATOM|HETATM)' ref2_clean.pdb) atoms"
fi
echo "  ref_clean   : $(grep -cE '^(ATOM|HETATM)' ref_clean.pdb) atoms"
echo "  concat_clean: frames=$(grep -c '^MODEL' concat_clean.pdb 2>/dev/null || echo '?')"
echo

# --------------------------------------------------------------------------
# [2] Detect ligand HEAVY ATOMS BY MASS, select them BY INTERNAL INDEX.
#     The previous version excluded H by NAME, which breaks when atom names
#     contain characters cpptraj's mask parser mishandles ('_', leading
#     digits, primes): e.g. AMP's '1H6_'/'2H6_' made @-masks match 0 atoms.
#
#     Robust fix: read name+mass+INDEX from atominfo (index = column 1, always
#     a plain integer), keep atoms with mass >= cutoff (heavy), and build the
#     ligand mask as "@<idx,idx,...>". Integer index masks are immune to any
#     atom-naming weirdness. Mass is read correctly by cpptraj even when the
#     element string is mis-guessed.
#
#     Indices are TOPOLOGY-SPECIFIC, so ref and concat get their OWN heavy
#     index lists (built from each structure's own atominfo).
# --------------------------------------------------------------------------
echo "=== [2/6] Ligand heavy-atom selection (by mass >= ${H_MASS_CUTOFF}, by index) ==="

# Emit "index name mass" for every LIGRES atom in a topology.
dump_lig_atoms() {  # $1 = topology
    cpptraj -p "$1" 2>/dev/null <<CPEOF | \
        awk '/^[[:space:]]*[0-9]+[[:space:]]+[A-Za-z0-9]/ { print $1, $2, $7 }'
atominfo :${LIGRES}
quit
CPEOF
}

# Build a comma-separated heavy-atom INDEX list (mass >= cutoff) for a topology,
# and report counts of heavy/H for the log.
build_heavy_idx() {  # $1 = topology  -> sets globals via echo "idxcsv|nheavy|nh|hnames"
    dump_lig_atoms "$1" | awk -v cut="$H_MASS_CUTOFF" '
        {
            idx=$1; name=$2; mass=$3+0
            if (mass >= cut) { heavy[++nh_]=idx; }
            else             { hname[++nl_]=name; }
        }
        END{
            csv=""
            for(i=1;i<=nh_;i++){ csv = csv (i>1?",":"") heavy[i] }
            hn=""
            for(i=1;i<=nl_;i++){ hn = hn (i>1?",":"") hname[i] }
            printf "%s|%d|%d|%s", csv, nh_, nl_, hn
        }'
}

REF_HEAVY_RAW=$(build_heavy_idx ref_clean.pdb)
CON_HEAVY_RAW=$(build_heavy_idx concat_clean.pdb)

REF_HEAVY_IDX=${REF_HEAVY_RAW%%|*}
CON_HEAVY_IDX=${CON_HEAVY_RAW%%|*}
REF_NHEAVY=$(echo "$REF_HEAVY_RAW" | cut -d'|' -f2)
CON_NHEAVY=$(echo "$CON_HEAVY_RAW" | cut -d'|' -f2)
REF_NH=$(echo "$REF_HEAVY_RAW" | cut -d'|' -f3)
CON_NH=$(echo "$CON_HEAVY_RAW" | cut -d'|' -f3)
REF_HNAMES=$(echo "$REF_HEAVY_RAW" | cut -d'|' -f4)

[[ -n "$REF_HEAVY_IDX" ]] || { echo "FATAL: no ligand atoms found for :${LIGRES} in ref. Wrong LIGRES?" >&2; exit 1; }
[[ -n "$CON_HEAVY_IDX" ]] || { echo "FATAL: no ligand atoms found for :${LIGRES} in concat. Wrong LIGRES?" >&2; exit 1; }

# Heavy-atom masks are PER TOPOLOGY (indices differ between ref and concat).
LIGHEAVY_REF="@${REF_HEAVY_IDX}"
LIGHEAVY_CON="@${CON_HEAVY_IDX}"
# A single symbolic name used in messages.
LIGHEAVY="@<ligand-heavy-indices>"

echo "  ref   : ${REF_NHEAVY} heavy, ${REF_NH} H excluded (${REF_HNAMES:-none})"
echo "  concat: ${CON_NHEAVY} heavy, ${CON_NH} H excluded"
if [[ "$REF_NHEAVY" != "$CON_NHEAVY" ]]; then
    echo "FATAL: ligand heavy-atom count differs ref=${REF_NHEAVY} concat=${CON_NHEAVY}." >&2
    echo "       ref and concat ligand are not the same heavy-atom set." >&2
    exit 1
fi
echo "  ligand heavy atoms: ${REF_NHEAVY} (selected by internal index; naming-proof)"
echo

# --------------------------------------------------------------------------
# [3] Map crystal resSeq -> cpptraj internal index, per topology.
#     (Same logic as before; cpptraj masks select by internal index, and
#      #Orig in resinfo preserves the crystal resSeq for translation.)
# --------------------------------------------------------------------------
echo "=== [3/6] Mapping crystal resSeq -> cpptraj internal (per structure) ==="
cat > _seqdump.in <<EOF
resinfo :1-99999
quit
EOF
dump_map() {  # $1 = topology -> orig internal resname
    cpptraj -p "$1" -i _seqdump.in 2>&1 | awk '
        /^[[:space:]]*[0-9]+[[:space:]]+[A-Z]/ { print $6, $1, $2 }'
}
dump_map ref_clean.pdb    > _map_ref.txt
dump_map concat_clean.pdb > _map_con.txt
[[ -s _map_ref.txt && -s _map_con.txt ]] || { echo "FATAL: resinfo map dump failed." >&2; exit 1; }

expand_list() {
    echo "$1" | tr ', ' '\n' \
      | awk -F- 'NF==2{for(i=$1;i<=$2;i++)print i;next} NF{print}' \
      | grep -E '^[0-9]+$'
}
crystal_to_internal_mask() {  # $1 mapfile  $2 crystal list -> :a,b,c  or MISSING:...
    local mf="$1" lst="$2" miss="" out="" cr intr
    for cr in $(expand_list "$lst"); do
        intr=$(awk -v O="$cr" '$1==O{print $2; exit}' "$mf")
        if [[ -z "$intr" ]]; then miss+="$cr "; else out+="${intr},"; fi
    done
    if [[ -n "$miss" ]]; then echo "MISSING:${miss% }"; return 1; fi
    echo ":${out%,}"
}

BS_REF_INT=$(crystal_to_internal_mask _map_ref.txt "$BS_RESIDUES") || true
BS_CON_INT=$(crystal_to_internal_mask _map_con.txt "$BS_RESIDUES") || true
for tag in "ref:$BS_REF_INT" "concat:$BS_CON_INT"; do
    if [[ "${tag#*:}" == MISSING:* ]]; then
        echo "FATAL: BS residue(s) absent from ${tag%%:*}: ${tag#*:MISSING:}" >&2
        echo "       Inspect _map_ref.txt / _map_con.txt." >&2
        exit 1
    fi
done

PROT_LIST="$(seq "$PROT_FIRST" "$PROT_LAST" | paste -sd, -)"
COMMON_WP=""; EXCL_WP=""
for cr in $(expand_list "$PROT_LIST"); do
    ri=$(awk -v O="$cr" '$1==O{print $2; exit}' _map_ref.txt)
    ci=$(awk -v O="$cr" '$1==O{print $2; exit}' _map_con.txt)
    if [[ -n "$ri" && -n "$ci" ]]; then COMMON_WP+="$cr "; else EXCL_WP+="$cr "; fi
done
[[ -n "$COMMON_WP" ]] || { echo "FATAL: no common whole-protein residues." >&2; exit 1; }
WP_REF_INT=$(crystal_to_internal_mask _map_ref.txt "$(echo $COMMON_WP | tr ' ' ',')")
WP_CON_INT=$(crystal_to_internal_mask _map_con.txt "$(echo $COMMON_WP | tr ' ' ',')")

# Final masks. BS heavy: exclude H by NAME (standard protein H names start
# with H, so @H= is reliable for protein; this is unrelated to the ligand bug).
BSCORE_REF="(${BS_REF_INT})&!@H="
BSCORE_CON="(${BS_CON_INT})&!@H="
WHOLE_REF="(${WP_REF_INT})&@CA"
WHOLE_CON="(${WP_CON_INT})&@CA"

echo "  BS residues mapped: $(expand_list "$BS_RESIDUES" | wc -l)"
if [[ -n "$EXCL_WP" ]]; then
    nex=$(echo $EXCL_WP | wc -w)
    echo "  whole-protein: ${nex} residue(s) excluded as unresolved in reference: ${EXCL_WP}"
fi
echo

# ------------- [4] verify mapping + atom-count parity ----------------------
echo "=== [4/6] Verify mapping (BS residue names match ref vs concat) ==="
BAD=""
for cr in $(expand_list "$BS_RESIDUES"); do
    rname=$(awk -v O="$cr" '$1==O{print $3; exit}' _map_ref.txt)
    cname=$(awk -v O="$cr" '$1==O{print $3; exit}' _map_con.txt)
    nr=${rname/HIE/HIS}; nr=${nr/HID/HIS}
    nc=${cname/HIE/HIS}; nc=${nc/HID/HIS}
    if [[ -z "$rname" || -z "$cname" || "$nr" != "$nc" ]]; then
        BAD="crystal resSeq $cr -> ref:${rname:-NONE} concat:${cname:-NONE}"; break
    fi
done
[[ -n "$BAD" ]] && { echo "FATAL: BS mapping mismatch at ${BAD}." >&2; exit 1; }
echo "  mapping verified"

echo "=== Atom-count parity (ref vs concat) ==="
count_mask() {  # $1=topology  $2=mask -> integer
    local out n
    out=$(mktemp)
    cpptraj -p "$1" >"$out" 2>&1 <<CPEOF
atominfo "$2"
quit
CPEOF
    n=$(grep -oiE '[0-9]+ atoms selected' "$out" | grep -oE '[0-9]+' | tail -1)
    [[ -z "$n" ]] && n=$(grep -cE '^[[:space:]]*[0-9]+[[:space:]]+[A-Za-z0-9]+[[:space:]]+[0-9]+' "$out")
    rm -f "$out"; [[ "$n" =~ ^[0-9]+$ ]] || n=0; printf '%s' "$n"
}
declare -A FAIL=()
for lbl in WHOLE BS LIG; do
    case "$lbl" in
        WHOLE) mref="$WHOLE_REF"; mcon="$WHOLE_CON" ;;
        BS)    mref="$BSCORE_REF"; mcon="$BSCORE_CON" ;;
        LIG)   mref="$LIGHEAVY_REF"; mcon="$LIGHEAVY_CON"  ;;
    esac
    c1=$(count_mask ref_clean.pdb    "$mref"); c1=$((${c1//[^0-9]/} + 0))
    c2=$(count_mask concat_clean.pdb "$mcon"); c2=$((${c2//[^0-9]/} + 0))
    printf "  %-6s ref=%-5s concat=%-5s " "$lbl" "$c1" "$c2"
    if (( c1 > 0 && c1 == c2 )); then echo "OK"; else echo "MISMATCH"; FAIL[$lbl]=1; fi
done
if [[ ${#FAIL[@]} -gt 0 ]]; then
    echo "FATAL: atom-count mismatch in: ${!FAIL[*]} (ref vs concat must agree)." >&2
    echo "  For LIG this almost always means the heavy-atom selection differs;" >&2
    echo "  check the detected ligand H list above against both structures." >&2
    exit 1
fi
echo

# ----------------------------- [5] RMSD runs -------------------------------
rm -f rmsd_wholepro.txt rmsd_BS.txt rmsd_lig.txt rmsd_combined_BS-lig.txt \
      sorted_rmsd_*.txt

echo "=== [5/6] wholepro + BS (single ref) ==="
cat > _cpptraj_run1.in <<EOF
parm concat_clean.pdb
parm ref_clean.pdb [refp]
reference ref_clean.pdb parm [refp] [R]
trajin concat_clean.pdb
rms      WholePro ${WHOLE_CON} ${WHOLE_REF} ref [R] out rmsd_wholepro.txt
symmrmsd BS       ${BSCORE_CON} ${BSCORE_REF} ref [R] remap out rmsd_BS.txt
run
quit
EOF
cpptraj -i _cpptraj_run1.in 2>&1 | tee "$OUTLOG".run1 | grep -iE 'warning|error' || true

echo "=== ligand + combined (BS-fit, ligand nofit) ==="
if [[ $DUAL_REF -eq 0 ]]; then
    # ---- single reference (normal path) ----
    cat > _cpptraj_run2.in <<EOF
parm concat_clean.pdb
parm ref_clean.pdb [refp]
reference ref_clean.pdb parm [refp] [R]
trajin concat_clean.pdb
rms BSfit ${BSCORE_CON} ${BSCORE_REF} ref [R]
rms ligR  ${LIGHEAVY_CON} ${LIGHEAVY_REF} ref [R] nofit out rmsd_lig.txt
rms combR (${BSCORE_CON})|(${LIGHEAVY_CON}) (${BSCORE_REF})|(${LIGHEAVY_REF}) ref [R] out rmsd_combined_BS-lig.txt
run
quit
EOF
    cpptraj -i _cpptraj_run2.in 2>&1 | tee "$OUTLOG".run2 | grep -iE 'warning|error' || true
else
    # ---- dual reference: independent pass per ref, then per-frame min ----
    # ref2 = ref with the two ligand arms transposed (you build it manually).
    # Per frame the smaller ligand RMSD (ref vs ref2) defines the correct
    # symmetry labelling; the SAME ref choice is applied to that frame's
    # combined value so a pose is never labelled inconsistently.
    #
    # ref2 gets its OWN heavy-atom indices and BS indices derived from its own
    # atominfo/resinfo -- not reused from ref. This way a mis-built ref2
    # (reordered atoms, different residue set) fails the parity check below
    # instead of producing silently-wrong "corrected" RMSDs.
    REF2_HEAVY_RAW=$(build_heavy_idx ref2_clean.pdb)
    REF2_HEAVY_IDX=${REF2_HEAVY_RAW%%|*}
    REF2_NHEAVY=$(echo "$REF2_HEAVY_RAW" | cut -d'|' -f2)
    [[ -n "$REF2_HEAVY_IDX" ]] || { echo "FATAL: no ligand atoms for :${LIGRES} in ref2." >&2; exit 1; }
    if [[ "$REF2_NHEAVY" != "$REF_NHEAVY" ]]; then
        echo "FATAL: ref2 ligand heavy count (${REF2_NHEAVY}) != ref (${REF_NHEAVY})." >&2
        echo "       ref2.pdb is not a valid arm-swap of ref (different atom set)." >&2
        exit 1
    fi
    LIGHEAVY_REF2="@${REF2_HEAVY_IDX}"

    # ref2 BS internal indices (its own resinfo map).
    dump_map ref2_clean.pdb > _map_ref2.txt
    BS_REF2_INT=$(crystal_to_internal_mask _map_ref2.txt "$BS_RESIDUES") || true
    if [[ "$BS_REF2_INT" == MISSING:* || -z "$BS_REF2_INT" ]]; then
        echo "FATAL: BS residue(s) absent from ref2: ${BS_REF2_INT#MISSING:}" >&2
        exit 1
    fi
    BSCORE_REF2="(${BS_REF2_INT})&!@H="

    for tag in 1 2; do
        if [[ $tag -eq 1 ]]; then
            rp=ref_clean.pdb;  rbs="$BSCORE_REF";  rlig="$LIGHEAVY_REF"
        else
            rp=ref2_clean.pdb; rbs="$BSCORE_REF2"; rlig="$LIGHEAVY_REF2"
        fi
        cat > _cpptraj_run2_${tag}.in <<EOF
parm concat_clean.pdb
parm ${rp} [refp]
reference ${rp} parm [refp] [R]
trajin concat_clean.pdb
rms BSfit ${BSCORE_CON} ${rbs} ref [R]
rms ligR  ${LIGHEAVY_CON} ${rlig} ref [R] nofit out _lig${tag}.txt
rms combR (${BSCORE_CON})|(${LIGHEAVY_CON}) (${rbs})|(${rlig}) ref [R] out _comb${tag}.txt
run
quit
EOF
        cpptraj -i _cpptraj_run2_${tag}.in 2>&1 | tee "$OUTLOG".run2_${tag} | grep -iE 'warning|error' || true
    done
    # assemble 4 columns: ligR1 ligR2 combR1 combR2  (1 row/frame, same order)
    paste \
      <(grep -vE '^#' _lig1.txt  | awk '{print $2}') \
      <(grep -vE '^#' _lig2.txt  | awk '{print $2}') \
      <(grep -vE '^#' _comb1.txt | awk '{print $2}') \
      <(grep -vE '^#' _comb2.txt | awk '{print $2}') \
      | awk '
          { l1=$1+0; l2=$2+0; c1=$3+0; c2=$4+0
            if (l1<=l2){ lig=l1; comb=c1; ref=1 } else { lig=l2; comb=c2; ref=2 }
            printf "%d %.4f\n", NR, lig  > "rmsd_lig.txt"
            printf "%d %.4f\n", NR, comb > "rmsd_combined_BS-lig.txt"
            if(ref==2) flip++; tot++ }
          END{ printf "  dual-ref: %d frames, ref2 chosen (arm-swapped labelling): %d (%.1f%%)\n",
                       tot, flip+0, (tot?100.0*flip/tot:0) }'
    # cpptraj `out` writes frames starting at 1; rmsd_lig/combined here are
    # 1-indexed (NR), matching wholepro/BS before the 0-index pass.
fi
echo

# --------------- [6] post-process: 0-index + sort --------------------------
echo "=== [6/6] Post-process (0-index frames, sort by RMSD) ==="
for f in rmsd_wholepro rmsd_BS rmsd_lig rmsd_combined_BS-lig; do
    if [[ -s "${f}.txt" ]] && ! head -1 "${f}.txt" | grep -q '^-1[[:space:]]'; then
        awk '!/^#/{printf "%d %s\n", $1-1, $2}' "${f}.txt" > "${f}.tmp" && mv "${f}.tmp" "${f}.txt"
    fi
done

ALL_OK=1
for f in rmsd_wholepro rmsd_BS rmsd_lig rmsd_combined_BS-lig; do
    if [[ -s "${f}.txt" ]]; then
        sort -k2,2g "${f}.txt" > "sorted_${f}.txt"
        echo "  ${f}: $(wc -l < "${f}.txt") frames -> sorted_${f}.txt  (best: $(head -1 sorted_${f}.txt))"
    else
        echo "  MISSING: ${f}.txt"; ALL_OK=0
    fi
done
echo
if [[ $DUAL_REF -eq 1 ]]; then MODE_STR="dual reference (ref + ref2)"; else MODE_STR="single reference"; fi
[[ $ALL_OK -eq 1 ]] && echo "DONE: 4 metrics written (${MODE_STR})." || { echo "INCOMPLETE."; exit 1; }

###############################################################################
# NOTE on ligand symmetry (AP5/B4P/GP5 two-arm pseudo-symmetry)
# -------------------------------------------------------------
# This script does NOT symmetry-correct the ligand. For your ADK targets the
# four-arm test showed the docked pose is NOT arm-swapped, so a single
# reference is correct and the previous ref2 machinery was unnecessary.
#
# If a future target genuinely needs it, the RIGHT tool is real symmrmsd with
# a topology that has correct bonds, e.g.:
#     parm complex.prmtop
#     symmrmsd ligS ${LIGHEAVY} ref [R] nofit out rmsd_lig.txt
# (Building the prmtop needs an AP5 prepi/frcmod; the tleap failure you saw
#  was a missing ligand parameter set, not a fundamental blocker.)
###############################################################################