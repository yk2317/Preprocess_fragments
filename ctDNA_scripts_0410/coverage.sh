#!/bin/bash

BAMFILE=$1
OUTDIR=$2

BASENAME=$(basename $BAMFILE .bam)

samtools depth -a $BAMFILE | \
awk '$1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$/ {print $1"\t"($2-1)"\t"$2"\t"$3}' > "${OUTDIR}/${BASENAME}_coverage.bed"

