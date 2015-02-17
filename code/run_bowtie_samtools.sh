#!/bin/bash
set -e
set -o pipefail

READS=$1
BASE=$(basename $READS .fna.gz)
bowtie2 preg.genes.fna.bowtieDB -f <(zcat $READS) -p 16 -S $BASE.aligned2genes.sam

samtools view -bS $BASE.aligned2genes.sam > $BASE.aligned2genes.bam
rm $BASE.aligned2genes.sam
samtools sort $BASE.aligned2genes.bam $BASE.aligned2genes.sorted
rm $BASE.aligned2genes.bam

samtools index $BASE.aligned2genes.sorted.bam
samtools idxstats $BASE.aligned2genes.sorted.bam > $BASE.mapped2genes.txt
