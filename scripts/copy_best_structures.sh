#!/bin/bash

# Create target directory if it doesn't exist
target_dir="selected_best_structures"
mkdir -p "$target_dir"

# Read each line from the sorted results file
while IFS= read -r line; do
  # Example line:
  # results_0/A0A0J9X1X4.B99990010.pdb: -21756.84961

  # Extract the full file path (everything before the colon)
  file_path=$(echo "$line" | cut -d':' -f1)

  # Extract the folder name (e.g., results_0)
  folder=$(echo "$file_path" | cut -d'/' -f1)

  # Remove the prefix "results_" to get just the numeric part
  index=${folder#results_}

  # Copy the pdb file into the target directory, renaming it as <index>.pdb
  cp "$file_path" "$target_dir/${index}.pdb"
done < lowest_dope_scores.txt

echo "Best structures have been copied into $target_dir with names matching their results folder."
