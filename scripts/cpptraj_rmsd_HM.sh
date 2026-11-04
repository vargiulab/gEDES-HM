#!/usr/bin/env bash
###############################################################################
# cpptraj_rmsd_HM.sh
#
# Single-reference RMSD pipeline for HOMOLOGY MODELS (no ligand).
#
# Computes WholePro (CA) and BS (heavy-atom) RMSD vs ref.pdb for two concat
# files in one invocation:
#   concat_all.pdb       -> rmsd_wholepro_all.txt       rmsd_BS_all.txt
#   concat_selected.pdb  -> rmsd_wholepro_selected.txt  rmsd_BS_selected.txt
# plus sorted_* counterparts.
#
# Inputs via env vars:
#   BS_RESIDUES, PROT_FIRST, PROT_LAST
#
# Requires in working dir:
#   ref.pdb, concat_all.pdb, concat_selected.pdb
#
# NOTE: sanitize() strips ALL hydrogens (protein + ligand) so heavy-atom
#       parity holds against heavy-atom-only concats (e.g. bioemu output).
###############################################################################

set +H
set -uo pipefail

# --------------------------- CONFIG ----------------------------------------
BS_RESIDUES="${BS_RESIDUES:-}"
PROT_FIRST="${PROT_FIRST:-1}"
PROT_LAST="${PROT_LAST:-99999}"

REF="${REF:-ref.pdb}"
CONCAT_ALL="${CONCAT_ALL:-concat_all.pdb}"
CONCAT_SEL="${CONCAT_SEL:-concat_selected.pdb}"
OUTLOG="cpptraj_rmsd_run.log"

# --------------------- environment / input guards --------------------------
command -v cpptraj >/dev/null 2>&1 || {
    echo "FATAL: cpptraj not in PATH. Load AmberTools." >&2; exit 1; }
for f in "$REF" "$CONCAT_ALL" "$CONCAT_SEL"; do
    [[ -e "$f" ]] || { echo "FATAL: required file '$f' not found in $(pwd)" >&2; exit 1; }
done
[[ -n "$BS_RESIDUES" ]] || { echo "FATAL: BS_RESIDUES is empty." >&2; exit 1; }

echo "cpptraj: $(command -v cpptraj)"
echo "Target : PROT=$PROT_FIRST-$PROT_LAST"
echo "BS res : $BS_RESIDUES"
echo "Ref    : $REF"
echo "Concat : $CONCAT_ALL , $CONCAT_SEL"
echo

# ----------------------------- sanitize ------------------------------------
# Pad records to 80 cols, drop alt-loc cruft, and STRIP ALL HYDROGENS.
# Hydrogen detection is by atom name (cols 13-16): matches H, HE, HZ1,
# 1HH1, 2HH2, HG, etc. via /^[0-9]*H/. This guarantees heavy-atom parity
# regardless of whether the element column or @H= wildcard behaves.
sanitize() {  # $1 in  $2 out
    sed 's/\r$//' "$1" \
      | awk '/^(ATOM|HETATM)/{
                 name=substr($0,13,4); gsub(/ /,"",name);
                 if (name ~ /^[0-9]*H/) next;        # drop hydrogens
                 l=$0; while(length(l)<80) l=l" "; print substr(l,1,80); next }
             /^(MODEL|ENDMDL|TER|END|CRYST1|CONECT)/{ print }' > "$2"
}

echo "=== [1/4] Sanitizing ref / concat files (stripping hydrogens) ==="
sanitize "$REF"        ref_clean.pdb
sanitize "$CONCAT_ALL" concat_all_clean.pdb
sanitize "$CONCAT_SEL" concat_selected_clean.pdb
echo "  ref_clean             : $(grep -cE '^(ATOM|HETATM)' ref_clean.pdb) atoms"
echo "  concat_all_clean      : $(grep -cE '^(ATOM|HETATM)' concat_all_clean.pdb) atoms"
echo "  concat_selected_clean : $(grep -cE '^(ATOM|HETATM)' concat_selected_clean.pdb) atoms"
# sanity: confirm no H survived
for f in ref_clean concat_all_clean concat_selected_clean; do
    h=$(grep -E '^(ATOM|HETATM)' "${f}.pdb" | awk '{n=substr($0,13,4); gsub(/ /,"",n); if(n ~ /^[0-9]*H/) c++} END{print c+0}')
    (( h == 0 )) || echo "  WARNING: ${f}.pdb still has ${h} H atom(s)"
done
echo

# --------------------------------------------------------------------------
# Map crystal resSeq -> cpptraj internal index (per topology)
# --------------------------------------------------------------------------
echo "=== [2/4] Mapping crystal resSeq -> cpptraj internal (per structure) ==="
cat > _seqdump.in <<EOF
resinfo :1-99999
quit
EOF

dump_map() {  # $1 = topology -> stdout: orig internal resname
    cpptraj -p "$1" -i _seqdump.in 2>&1 | awk '
        /^[[:space:]]*[0-9]+[[:space:]]+[A-Z]/ {
            internal=$1; name=$2; orig=$6;
            print orig, internal, name
        }'
}

dump_map ref_clean.pdb                > _map_ref.txt
dump_map concat_all_clean.pdb         > _map_all.txt
dump_map concat_selected_clean.pdb    > _map_sel.txt
for m in _map_ref.txt _map_all.txt _map_sel.txt; do
    [[ -s "$m" ]] || { echo "FATAL: resinfo map dump failed ($m)." >&2; exit 1; }
done

expand_list() {
    echo "$1" | tr ', ' '\n' \
      | awk -F- 'NF==2{for(i=$1;i<=$2;i++)print i;next} NF{print}' \
      | grep -E '^[0-9]+$'
}

crystal_to_internal_mask() {  # $1=mapfile  $2=crystal list
    local mf="$1" lst="$2" miss="" out="" cr intr
    for cr in $(expand_list "$lst"); do
        intr=$(awk -v O="$cr" '$1==O{print $2; exit}' "$mf")
        if [[ -z "$intr" ]]; then miss+="$cr "; else out+="${intr},"; fi
    done
    if [[ -n "$miss" ]]; then echo "MISSING:${miss% }"; return 1; fi
    echo ":${out%,}"
}

BS_REF_INT=$(crystal_to_internal_mask _map_ref.txt "$BS_RESIDUES") || true
BS_ALL_INT=$(crystal_to_internal_mask _map_all.txt "$BS_RESIDUES") || true
BS_SEL_INT=$(crystal_to_internal_mask _map_sel.txt "$BS_RESIDUES") || true
PROT_LIST="$(seq "$PROT_FIRST" "$PROT_LAST" | paste -sd, -)"

for tag in "ref:$BS_REF_INT" "concat_all:$BS_ALL_INT" "concat_selected:$BS_SEL_INT"; do
    if [[ "${tag#*:}" == MISSING:* ]]; then
        echo "FATAL: BS residue(s) absent from ${tag%%:*}: ${tag#*:MISSING:}" >&2
        exit 1
    fi
done

# whole-protein: residues that exist in ref AND in BOTH concat files
COMMON_WP=""
EXCL_WP=""
for cr in $(expand_list "$PROT_LIST"); do
    ri=$(awk  -v O="$cr" '$1==O{print $2; exit}' _map_ref.txt)
    ai=$(awk  -v O="$cr" '$1==O{print $2; exit}' _map_all.txt)
    si=$(awk  -v O="$cr" '$1==O{print $2; exit}' _map_sel.txt)
    if [[ -n "$ri" && -n "$ai" && -n "$si" ]]; then
        COMMON_WP+="$cr "
    else
        EXCL_WP+="$cr "
    fi
done
[[ -n "$COMMON_WP" ]] || { echo "FATAL: no common whole-protein residues across ref + both concats." >&2; exit 1; }
WP_REF_INT=$(crystal_to_internal_mask _map_ref.txt "$(echo $COMMON_WP | tr ' ' ',')")
WP_ALL_INT=$(crystal_to_internal_mask _map_all.txt "$(echo $COMMON_WP | tr ' ' ',')")
WP_SEL_INT=$(crystal_to_internal_mask _map_sel.txt "$(echo $COMMON_WP | tr ' ' ',')")

# Hydrogens already stripped at sanitize; &!@H= kept as a harmless guard.
BSCORE_REF="${BS_REF_INT}&!@H="
BSCORE_ALL="${BS_ALL_INT}&!@H="
BSCORE_SEL="${BS_SEL_INT}&!@H="
WHOLE_REF="${WP_REF_INT}&@CA"
WHOLE_ALL="${WP_ALL_INT}&@CA"
WHOLE_SEL="${WP_SEL_INT}&@CA"

echo "  BS residues mapped: $(expand_list "$BS_RESIDUES" | wc -l) (all present in ref + both concats)"
if [[ -n "$EXCL_WP" ]]; then
    nex=$(echo $EXCL_WP | wc -w)
    echo "  whole-protein: ${nex} residue(s) excluded as unresolved in at least one structure"
    echo "    crystal resSeq: ${EXCL_WP}"
else
    echo "  whole-protein: all PROT residues present in ref + both concats"
fi
echo

# ------------- verify mapping correctness ----------------------------------
echo "=== [3/4] Verifying mapping + atom-count parity ==="
BAD=""
for cr in $(expand_list "$BS_RESIDUES"); do
    rname=$(awk -v O="$cr" '$1==O{print $3; exit}' _map_ref.txt)
    aname=$(awk -v O="$cr" '$1==O{print $3; exit}' _map_all.txt)
    sname=$(awk -v O="$cr" '$1==O{print $3; exit}' _map_sel.txt)
    nr=${rname/HIE/HIS}; nr=${nr/HID/HIS}
    na=${aname/HIE/HIS}; na=${na/HID/HIS}
    ns=${sname/HIE/HIS}; ns=${ns/HID/HIS}
    if [[ -z "$rname" || -z "$aname" || -z "$sname" || "$nr" != "$na" || "$nr" != "$ns" ]]; then
        BAD="crystal resSeq $cr -> ref:${rname:-NONE}  all:${aname:-NONE}  sel:${sname:-NONE}"
        break
    fi
done
if [[ -n "$BAD" ]]; then
    echo "FATAL: mapping mismatch at ${BAD}." >&2
    exit 1
fi
echo "  mapping verified: every BS crystal residue resolves to the same residue name"

# Count atoms in a mask using cpptraj 'select' (stable 'Selected N atoms' line).
count_mask_cpptraj() {  # $1=topology  $2=mask -> integer
    local out n
    out=$(mktemp)
    cpptraj -p "$1" >"$out" 2>&1 <<CPEOF
select $2
quit
CPEOF
    n=$(grep -oiE 'Selected [0-9]+ atoms' "$out" | grep -oE '[0-9]+' | tail -1)
    [[ -z "$n" ]] && n=$(grep -oiE '[0-9]+ atoms selected' "$out" | grep -oE '[0-9]+' | tail -1)
    rm -f "$out"
    [[ "$n" =~ ^[0-9]+$ ]] || n=0
    printf '%s' "$n"
}

declare -A FAIL=()
# WHOLE
cR=$(count_mask_cpptraj ref_clean.pdb                "$WHOLE_REF")
cA=$(count_mask_cpptraj concat_all_clean.pdb         "$WHOLE_ALL")
cS=$(count_mask_cpptraj concat_selected_clean.pdb    "$WHOLE_SEL")
cR=$((cR+0)); cA=$((cA+0)); cS=$((cS+0))
printf "  %-6s ref=%-5s all=%-5s sel=%-5s " "WHOLE" "$cR" "$cA" "$cS"
if (( cR > 0 && cR == cA && cR == cS )); then echo "OK"; else echo "MISMATCH"; FAIL[WHOLE]=1; fi
# BS
cR=$(count_mask_cpptraj ref_clean.pdb                "$BSCORE_REF")
cA=$(count_mask_cpptraj concat_all_clean.pdb         "$BSCORE_ALL")
cS=$(count_mask_cpptraj concat_selected_clean.pdb    "$BSCORE_SEL")
cR=$((cR+0)); cA=$((cA+0)); cS=$((cS+0))
printf "  %-6s ref=%-5s all=%-5s sel=%-5s " "BS" "$cR" "$cA" "$cS"
if (( cR > 0 && cR == cA && cR == cS )); then echo "OK"; else echo "MISMATCH"; FAIL[BS]=1; fi

if (( ${#FAIL[@]} > 0 )); then
    echo "FATAL: atom-count mismatch in: ${!FAIL[*]} (ref and both concats must agree)." >&2
    exit 1
fi
echo

# ----------------------------- RMSD runs -----------------------------------
rm -f rmsd_wholepro_all.txt      rmsd_BS_all.txt \
      rmsd_wholepro_selected.txt rmsd_BS_selected.txt \
      sorted_rmsd_wholepro_all.txt      sorted_rmsd_BS_all.txt \
      sorted_rmsd_wholepro_selected.txt sorted_rmsd_BS_selected.txt

run_rmsd() {  # $1=tag (all|selected)  $2=concat_clean  $3=WHOLE_mask_con  $4=BS_mask_con
    local tag="$1" con="$2" wmask="$3" bmask="$4"
    echo "  --- $tag : $con ---"
    cat > "_cpptraj_run_${tag}.in" <<EOF
parm $con
parm ref_clean.pdb [refp]
reference ref_clean.pdb parm [refp] [R]
trajin $con
rms      WholePro_${tag} ${wmask} ${WHOLE_REF} ref [R] out rmsd_wholepro_${tag}.txt
symmrmsd BS_${tag}       ${bmask} ${BSCORE_REF} ref [R] remap out rmsd_BS_${tag}.txt
run
quit
EOF
    cpptraj -i "_cpptraj_run_${tag}.in" 2>&1 | tee -a "${OUTLOG}.run"
}

echo "=== [4/4] Running single-reference RMSD calculations ==="
: > "${OUTLOG}.run"
run_rmsd all      concat_all_clean.pdb      "$WHOLE_ALL" "$BSCORE_ALL"
run_rmsd selected concat_selected_clean.pdb "$WHOLE_SEL" "$BSCORE_SEL"

auto_reindex() {  # convert frame column to 0-index if not already
    local f="$1"
    [[ -s "$f" ]] || return 0
    if ! head -1 "$f" | grep -q '^-1[[:space:]]'; then
        awk '!/^#/{printf "%d %s\n", $1-1, $2}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    fi
}

auto_reindex rmsd_wholepro_all.txt
auto_reindex rmsd_BS_all.txt
auto_reindex rmsd_wholepro_selected.txt
auto_reindex rmsd_BS_selected.txt

ALL_OK=1
for f in rmsd_wholepro_all rmsd_BS_all rmsd_wholepro_selected rmsd_BS_selected; do
    if [[ -s "${f}.txt" ]]; then
        sort -k2,2g "${f}.txt" > "sorted_${f}.txt"
        echo "  ${f}: $(wc -l < "${f}.txt") frames -> sorted_${f}.txt"
    else
        echo "  MISSING: ${f}.txt"
        ALL_OK=0
    fi
done

[[ $ALL_OK -eq 1 ]] && echo "DONE: 8 files written." || { echo "INCOMPLETE."; exit 1; }
