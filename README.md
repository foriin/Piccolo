# Piccolo
A simple script that looks for piRNA clusters based on small RNA library coverage. It uses a simple heuristic for cluster identification and classification for uni- and double-stranded clusters. The only input it needs is the bam file resulting from mapping your small RNA reads onto your genome of choice.

<img src=https://github.com/foriin/Piccolo/blob/dev/resources/workfl.png width="600">

## Prerequisites
Bedtools

## Installation

```
git clone https://github.com/foriin/Piccolo.git
cd Piccolo
./find_clusters.sh <your_mapped_reads.bam>
```
## Running

```
Usage: find_clusters.sh -i <input_bam_file> [-m <cluster_min_length>] [-o <output_directory>] [-A] [-h]
  -i  Input BAM file (required).
  -m  Minimum length of clusters (default is 10000).
  -o  Output directory (default is current directory).
  -A  Keep all intermediate output files.
  -h  Display this help message.
```

## Tips

If you have multiple samples of your libraries, you can go two ways to merge the results: 1) be permissive and annotate overlaps of clusters in multiple replicates possibly expanding their boundaries or 2) be strict and treat only regions that intersect in all of the samples as 'proper' piRNA clusters. An example of going approach #1 using R would be:
```
library(GenomicRanges)
library(rtracklayer)
# Get a vector with paths to BED files with cluster for your samples
cldir <- dir("directory/with/Piccolo/clusters", full.names = T)
# Import BED files
clbcop <- lapply(cldir, function(x) import.bed(x))
# Unite and reduce all clusters, only clusters on the same strand would be merged
clbcop.red <- GenomicRanges::reduce(unlist(GRangesList(clbcop)))
# Save merged clusters as a BED
export.bed(clbcop.red, "/path/for/output.bed")
```

