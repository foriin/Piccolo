#!/usr/bin/env bash

# Default values
CL_MIN_LEN=10000
OUTPUT_DIR="."
KEEP_ALL=false

# Get the directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

print_help() {
    echo "Usage: find_clusters.sh -i <input_bam_file> [-m <cluster_min_length>] [-o <output_directory>] [-A] [-h]"
    echo "  -i  Input BAM file."
    echo "  -m  Minimum length of clusters (default is 10000)."
    echo "  -o  Output directory (default is current directory)."
    echo "  -A  Keep all intermediate output files."
    echo "  -h  Display this help message."
}

# Parse command-line options
while getopts 'i:m:o:Ah' flag; do
    case "${flag}" in
        i) INPUT_BAM="${OPTARG}" ;;
        m) CL_MIN_LEN="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        A) KEEP_ALL=true ;;
        h) print_help
           exit 0 ;;
        *) print_help
           exit 1 ;;
    esac
done

# Check for mandatory options
if [ -z "$INPUT_BAM" ]; then
    echo "Error: Input BAM file is required."
    print_help
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Prepare filenames
PREFIX=$(basename "$INPUT_BAM" .bam)
PLUS_BG="${OUTPUT_DIR}/${PREFIX}_cov_plus.bg"
MINUS_BG="${OUTPUT_DIR}/${PREFIX}_cov_minus.bg"
COMBINED_BG="${OUTPUT_DIR}/${PREFIX}_cov.bg"
CLUSTERS_BED="${OUTPUT_DIR}/${PREFIX}_clusters.bed"
LARGE_CLUSTERS_BED="${OUTPUT_DIR}/${PREFIX}_clusters_large.bed"
FINAL_OUTPUT="${OUTPUT_DIR}/${PREFIX}_clustersClass.bed"

# Generate coverage files
bedtools genomecov -ibam "$INPUT_BAM" -bg -strand + | awk 'BEGIN {OFS = "\t"}{print $1, $2, $3, ".", $4, "+"}' > "$PLUS_BG"
bedtools genomecov -ibam "$INPUT_BAM" -bg -strand - | awk 'BEGIN {OFS = "\t"}{print $1, $2, $3, ".", $4, "-"}' > "$MINUS_BG"

# Combine coverage files
cat "$PLUS_BG" "$MINUS_BG" > "$COMBINED_BG"

# Merge and classify clusters
bedtools sort -i "$COMBINED_BG" | bedtools merge -i - -s -d 50 -c 5,6 -o mean,distinct > "$CLUSTERS_BED"
awk -v minLen="$CL_MIN_LEN" '$3 - $2 >= minLen{print}' "$CLUSTERS_BED" > "$LARGE_CLUSTERS_BED"

# Run custom C++ program (clustercall)
"$SCRIPT_DIR/bin/clustercall" "$LARGE_CLUSTERS_BED" "$FINAL_OUTPUT"

# Remove intermediate files if -A is not set
if [ "$KEEP_ALL" = false ]; then
    rm -f "$PLUS_BG" "$MINUS_BG" "$COMBINED_BG" "$CLUSTERS_BED" "$LARGE_CLUSTERS_BED"
fi

