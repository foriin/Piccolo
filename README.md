# Piccolo
A simple script that looks for piRNA clusters based on small RNA library coverage. It uses a simple heuristic for cluster identification and classification for uni- and double-stranded clusters. The only input it needs is the bam file resulting from mapping your small RNA reads onto your genome of choice.

<img src=https://github.com/foriin/Piccolo/blob/dev/resources/workflow.png width="600">

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

