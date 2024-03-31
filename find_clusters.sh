#!/bin/bash

# Check for input argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_bam_file> [cluster_min_length]"
    exit 1
fi

INPUT_BAM=$1
CL_MIN_LEN=${2:-10000}
PREFIX=$(basename "$INPUT_BAM" .bam)

# Generate coverage files
bedtools genomecov -ibam "$INPUT_BAM" -bg -strand + | awk 'BEGIN {OFS = "\t"}{print $1, $2, $3, ".", $4, "+"}' > "${PREFIX}_cov_plus.bg"
bedtools genomecov -ibam "$INPUT_BAM" -bg -strand - | awk 'BEGIN {OFS = "\t"}{print $1, $2, $3, ".", $4, "-"}' > "${PREFIX}_cov_minus.bg"

# Combine coverage files
cat "${PREFIX}_cov_plus.bg" "${PREFIX}_cov_minus.bg" > "${PREFIX}_cov.bg"

# Merge and classify clusters
bedtools sort -i "${PREFIX}_cov.bg" | bedtools merge -i - -s -d 50 -c 5,6 -o mean,distinct > "${PREFIX}_clusters.bed"
awk -v minLen="$CL_MIN_LEN" '$3 - $2 >= minLen{print}' "${PREFIX}_clusters.bed" > "${PREFIX}_clusters_large.bed"

# Run custom C++ program (clustercall)
./bin/clustercall "${PREFIX}_clusters_large.bed" "${PREFIX}_large_class.bed"

