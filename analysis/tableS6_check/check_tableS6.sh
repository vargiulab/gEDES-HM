#!/usr/bin/env bash
# Reproduce Table S6 from the ranking tables in docking/, for all four receptor
# conditions. Read-only; needs nothing but this archive and a POSIX shell.
#
#   bash check_tableS6.sh
#   ADK_ARCHIVE_ROOT=/path/to/archive bash check_tableS6.sh
#
# For each complex it takes the row of lowest RMSD_BM, i.e. column 6
# (pose_combined_BS-lig) of the ranking table. Columns are:
#   1 rank  2 cluster_id  3 size  4 avg_top3_score  5 pose_score
#   6 pose_combined_BS-lig  7 pose_lig_rmsd  8 pose_bs_rmsd
#   9 pose_wholepro_rmsd  10 pose_id  11 pose_fnat
#
# Note this is a minimum over clusters, not the representative pose printed in
# Table S6, which is chosen on all metrics together. Every value in Table S6
# therefore sits at or above the value printed here; the success counts, the
# percentages and the identity of the unrecovered complexes are what this
# script is meant to confirm.
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT="${ADK_ARCHIVE_ROOT:-$(cd "$HERE/../.." && pwd)}"
D="$ROOT/docking"
CUT="${CUT:-2.5}"
[ -d "$D" ] || { echo "no docking/ under $ROOT; set ADK_ARCHIVE_ROOT" >&2; exit 1; }

COMPLEXES="1 2_ADP 2_AMP 3 4 5 6 7 8 9_ADP 9_AMP 10 11 12 13 14 \
15_ADP182 15_ADP183 16 17 18 19 20 21_ADP214 21_ADP215 22"

# label : filename prefix in docking/ ("." means no prefix)
CONDS="HMgEDES:. HM4AKE_apo:APO_ HM1AKE_holo:HOLO-1ake_ HOLOexpt:redocking_HOLO-expt_"

for c in $CONDS; do
  label=${c%%:*}; pre=${c#*:}; [ "$pre" = "." ] && pre=""
  dir="$D"
  printf '\n===== %s =====\n' "$label"
  printf '%-12s %6s  %7s %7s %7s %7s  %6s\n' complex rank BM lig BS prot Fnat
  vals=""; n=0; ok=0; fails=""
  for k in $COMPLEXES; do
    s=${k%%_*}; sfx=${k#"$s"}; sfx=${sfx#_}
    f="$dir/${pre}${s}_cluster_ranks${sfx:+_$sfx}.tsv"
    [ -f "$f" ] || { printf '%-12s  MISSING %s\n' "$k" "${f#$ROOT/}"; continue; }
    set -- $(tail -n +2 "$f" | sort -t "$(printf '\t')" -k6,6n | head -1 \
             | awk -F'\t' '{printf "%s %s %s %s %s %s", $1,$6,$7,$8,$9,$11}')
    printf '%-12s %6s  %7.2f %7.2f %7.2f %7.2f  %6.2f\n' "$k" "$1" "$2" "$3" "$4" "$5" "$6"
    vals="$vals $2"; n=$((n+1))
    if awk -v a="$2" -v c="$CUT" 'BEGIN{exit !(a<=c)}'; then ok=$((ok+1)); else fails="$fails $k"; fi
  done
  med=$(printf '%s\n' $vals | sort -n | awk '{a[NR]=$1} END{if(!NR){print "NA";exit}
        if(NR%2) printf "%.2f", a[(NR+1)/2]; else printf "%.2f", (a[NR/2]+a[NR/2+1])/2}')
  printf -- '----\n'
  printf 'complexes: %s   RMSD_BM <= %s: %s (%.1f%%)   median: %s   not recovered:%s\n' \
    "$n" "$CUT" "$ok" "$(awk -v o=$ok -v n=$n 'BEGIN{print (n?100*o/n:0)}')" "$med" "${fails:- none}"
done

cat <<'EOF'

Table S6 of the SI reports, for the same four conditions in this order:
  24/26 (92.3%), 2/26 (7.7%), 24/26 (92.3%), 25/26 (96.2%)
with median RMSD_BM 2.1, 5.1, 1.1 and 0.6 A, the unrecovered complexes being
12 and 14, 12 and 14, and 12 alone for the redocking control.
EOF
