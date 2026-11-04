#!/usr/bin/env bash
# Recompute the Table S2 quantities (percentage of models within 2.5 A and the
# minimum RMSD, for the binding site and for the whole protein) directly from
# the archived model ensembles, for all 22 targets and all ligand variants.
#
# It reads, for every system <n> of proteins_structural_models/gEDES/<n>/:
#   run_calc_rmsd-cpptraj*.sh   the per-system binding-site definition
#   ref_clean*.pdb.gz           the experimental reference
#   concat_selected_clean.pdb.gz  the 101 DOPE-selected models
#   concat_all_clean.pdb.gz       all 1010 models
# and calls scripts/cpptraj_rmsd_HM.sh on each.
#
# The output, results.tsv, is what tableS2_recomputed.tsv holds; compare it
# with tableS2_published.tsv, which is Table S2 as printed in the SI.
#
# Requirements
#   cpptraj on PATH (AmberTools; tested with AmberTools25)
#
# Usage
#   export ADK_ARCHIVE_ROOT=/path/to/this/archive
#   [export ADK_WORKDIR=/path/to/scratch]      # default: a temporary directory
#   bash run_all.sh
set -u

: "${ADK_ARCHIVE_ROOT:?set ADK_ARCHIVE_ROOT to the top directory of this archive}"
command -v cpptraj >/dev/null 2>&1 || {
  echo "cpptraj not found on PATH; load AmberTools first" >&2; exit 1; }

Z="$ADK_ARCHIVE_ROOT"
W="${ADK_WORKDIR:-$(mktemp -d)}"
G="$Z/proteins_structural_models/gEDES"
mkdir -p "$W"
echo "archive : $Z"
echo "work dir: $W"

RESULTS="$W/results.tsv"
printf 'system\tsuffix\tBSsel_pct\tBSsel_min\tBSall_pct\tBSall_min\tWPsel_pct\tWPsel_min\tWPall_pct\tWPall_min\n' > "$RESULTS"

# gzcat is macOS-only; gzip -dc is portable
unz() { gzip -dc "$1"; }

pct_min() {   # % of models at or below 2.5 A, then the minimum, from a cpptraj RMSD file
  awk 'NR>0{n++; if($2<=2.5)k++; if(m==""||$2<m)m=$2} END{if(n)printf "%.1f\t%.1f", 100*k/n, m; else printf "NA\tNA"}' "$1"
}

for t in $(seq 1 22); do
  for rc in "$G/$t"/run_calc_rmsd-cpptraj*.sh; do
    [ -e "$rc" ] || continue
    base=$(basename "$rc" .sh)
    sfx="${base#run_calc_rmsd-cpptraj}"; sfx="${sfx#_}"
    d="$W/T${t}${sfx:+_$sfx}"
    mkdir -p "$d" || continue
    ( cd "$d" || exit 1
      refgz="$G/$t/ref_clean${sfx:+_$sfx}.pdb.gz"
      [ -e "$refgz" ] || { echo "MISSING ref $refgz" >&2; exit 1; }
      unz "$refgz"                            > ref.pdb
      unz "$G/$t/concat_selected_clean.pdb.gz" > concat_selected.pdb
      unz "$G/$t/concat_all_clean.pdb.gz"      > concat_all.pdb
      cp "$Z/scripts/cpptraj_rmsd_HM.sh" .
      BSR=$(sed -n 's/.*BS_RESIDUES="\([^"]*\)".*/\1/p' "$rc")
      PF=$(sed -n 's/.*PROT_FIRST=\([0-9]*\).*/\1/p' "$rc")
      PL=$(sed -n 's/.*PROT_LAST=\([0-9]*\).*/\1/p' "$rc")
      BS_RESIDUES="$BSR" PROT_FIRST="$PF" PROT_LAST="$PL" \
        bash cpptraj_rmsd_HM.sh > run.log 2>&1
    ) || { echo "failed $t $sfx" >&2; continue; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$t" "${sfx:-.}" \
      "$(pct_min "$d/rmsd_BS_selected.txt")"       "$(pct_min "$d/rmsd_BS_all.txt")" \
      "$(pct_min "$d/rmsd_wholepro_selected.txt")" "$(pct_min "$d/rmsd_wholepro_all.txt")" >> "$RESULTS"
    echo "done $t $sfx"
    rm -f "$d"/concat_all.pdb "$d"/concat_selected.pdb
  done
done
echo "ALLDONE -> $RESULTS"
