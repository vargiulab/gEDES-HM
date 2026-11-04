#!/bin/bash

# Output file to store results
output="lowest_dope_scores.txt"
> "$output"  # clear output file

# Loop over all results_* directories
for dir in results_*; do
  if [ -d "$dir" ]; then
    best_model=""
    best_score=0  # DOPE scores are negative, so any negative number is lower than 0
    
    for pdb in "$dir"/*.pdb; do
      # Skip concat.pdb files
      if [[ "$(basename "$pdb")" == "concat.pdb" ]]; then
        continue
      fi

      score=$(grep "REMARK   6 DOPE score:" "$pdb" | head -n 1 | awk '{print $NF}')
      if [ -n "$score" ]; then
        if awk "BEGIN {exit !($score < $best_score)}"; then
          best_score=$score
          best_model=$pdb
        fi
      fi
    done

    if [ -n "$best_model" ]; then
      echo "$best_model: $best_score" >> "$output"
    fi
  fi
done

# Sort the output file in natural (version) order so that results_2 comes before results_10, etc.
sort -V "$output" -o "$output"

echo "Sorted lowest DOPE scores have been written to $output"
