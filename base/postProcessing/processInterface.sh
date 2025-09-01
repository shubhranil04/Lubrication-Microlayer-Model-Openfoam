#!/bin/bash

# Utility to remove 'cell type', 'interfaceCentre.liquid:2' columns and duplicate rows.

# Directory containing the CSV files
input_dir="./interface"
output_dir="./processed_interface"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Process each CSV file in the directory
for file in "$input_dir"/*.csv; do
    # Extract the filename without the path
    filename=$(basename "$file")
   
    # Use awk to remove specific columns and then remove duplicate rows
    awk -F, '
        BEGIN { OFS = FS }
        NR == 1 { print $1, $3, $4 }  # Keep the header row with selected columns
        NR > 1 { print $1, $3, $4 | "sort -u" }  # Keep only unique rows for the rest
    ' "$file" > "$output_dir/$filename"

    echo "Processed: $file -> $output_dir/$filename"
done

echo "All files processed. Check the directory: $output_dir"
