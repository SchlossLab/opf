#!/bin/sh
 
#script to run digital normalization for single end sequences
#USAGE: sh diginorm.sh seqs.fna
#####
#files
SEQS=$1
#paths
KHMER_SCRIPTS="/home/kiverson/khmer/scripts" #path to khmer/scripts
 
#parameters
K=20
C=20
X_PRAM='1e10'

python $KHMER_SCRIPTS/normalize-by-median.py -k $K -C $C -x $X_PRAM --savetable $SEQS.kh $SEQS
python $KHMER_SCRIPTS/filter-abund.py -V $SEQS.kh $SEQS.keep
rm $SEQS.kh

