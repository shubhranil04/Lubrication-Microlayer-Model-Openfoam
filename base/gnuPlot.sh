#!/bin/bash

# Check if folder name is given
if [ -z "$1" ]; then
    echo "Usage: $0 <foldername> [ymax]"
    exit 1
fi

FOLDER="$1"
XML_FILE="$FOLDER/xml.csv"
HML_FILE="$FOLDER/hml.csv"
DATA_FILE="$FOLDER/plot_data.dat"

# Optional YMAX argument
YMAX="$2"
if [ -z "$YMAX" ]; then
    YMAX="1e-5"  # Default if not provided
fi

# Check if files exist
if [ ! -f "$XML_FILE" ] || [ ! -f "$HML_FILE" ]; then
    echo "Error: xml.csv or hml.csv not found in $FOLDER"
    exit 1
fi

# Merge the files line by line into one data file (X Y)
paste -d' ' "$XML_FILE" "$HML_FILE" > "$DATA_FILE"

# Find min and max of XML (first column)
XMIN=$(awk '{print $1}' "$DATA_FILE" | sort -n | head -1)
XMAX=$(awk '{print $1}' "$DATA_FILE" | sort -n | tail -1)

# Plot using ASCII in terminal with x-range set to data bounds
gnuplot <<EOF
set terminal dumb size 100,30 enhanced
set title "HML vs XML"
set xlabel "XML"
set ylabel "HML"
set xrange [$XMIN:$XMAX]
set yrange [0:$YMAX]
plot "$DATA_FILE" using 1:2 with lines title 'HML vs XML'
EOF
