#!/bin/bash
set -e
set -o pipefail

#datasets
DATASETS=(anterior_nares
buccal_mucosa
mid_vagina
palatine_tonsils
posterior_fornix
right_retroauricular_crease
subgingival_plaque
supragingival_plaque
throat
tongue_dorsum
vaginal_introitus)

#removed datasets (done previosuly)
#saliva
#stool
#attached_keratinized_gingiva
#left_retroauricular_crease

READS_LOC="ftp://public-ftp.hmpdacc.org/Illumina/"
GENES_LOC="ftp://public-ftp.hmpdacc.org/HMGI/"


#parallelize this
for i in ${DATASETS[*]};
do
	mkdir $i;
	mkdir $i/cluster

	wget ${READS_LOC}${i}/* -P $i/rawreads
	wget ${GENES_LOC}${i}/* -P $i/genes
done;


#everything past here is the same for all datasets


bzcat *_aa.fasta.bz2 > all.aa.genes.fasta && python ~/scripts/removeShortContigs.py all.aa.genes.fasta all.aa.genes100.fasta 100
formatdb -i all.aa.genes100.fasta -p T -n all.aa.genes100.fasta.blastDB
blastp -query all.aa.genes100.fasta -db all.aa.genes100.fasta.blastDB -out allvall.aa.genes100.out -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads 8

bzcat *_nucleotide.fasta.bz2 > all.nucleotide.genes.fasta && python ~/scripts/removeShortContigs.py all.nucleotide.genes.fasta all.nucleotide.genes300.fasta 300

bowtie-build all.nucleotide.genes300.fasta all.nucleotide.genes300.fasta.bowtieDB


READS=$1
BASE=$(basename $READS .tar.bz2)

bzcat $READS | tar -xvf -

bowtie ../genes/all.nucleotide.genes300.fasta.bowtieDB -1 $BASE/$BASE.denovo_duplicates_marked.trimmed.1.fastq -2 $BASE/$BASE.denovo_duplicates_marked.trimmed.2.fastq -p 16 -S $BASE.pe.aligned2genes.sam

samtools view -bS $BASE.pe.aligned2genes.sam > $BASE.pe.aligned2genes.bam
rm $BASE.pe.aligned2genes.sam
samtools sort $BASE.pe.aligned2genes.bam $BASE.pe.aligned2genes.sorted
rm $BASE.pe.aligned2genes.bam

samtools index $BASE.pe.aligned2genes.sorted.bam
samtools idxstats $BASE.pe.aligned2genes.sorted.bam > $BASE.pe.mapped2genes.txt
rm $BASE.pe.aligned2genes.sorted.bam


bowtie ../genes/all.nucleotide.genes300.fasta.bowtieDB -q $BASE/$BASE.denovo_duplicates_marked.trimmed.singleton.fastq -p 16 -S $BASE.se.aligned2genes.sam

samtools view -bS $BASE.se.aligned2genes.sam > $BASE.se.aligned2genes.bam
rm $BASE.se.aligned2genes.sam
samtools sort $BASE.se.aligned2genes.bam $BASE.se.aligned2genes.sorted
rm $BASE.se.aligned2genes.bam

samtools index $BASE.se.aligned2genes.sorted.bam
samtools idxstats $BASE.se.aligned2genes.sorted.bam > $BASE.se.mapped2genes.txt
rm $BASE.se.aligned2genes.sorted.bam

rm -rf $BASE


python ../../opf/code/makecountfilefrombowtie.py -f ../genes/all.aa.genes100.fasta -o saliva.genes.count -l *.txt
getMissingGenes.py kg.genes.count ../calledgenes/allvall.aa.genes100.out missingGenes.txt && cat kg.genes.count missingGenes.txt > kg.genes.ALL.count

mothur '#mgcluster(blast=/mnt/EXT/Schloss-data/kiverson/hmp_saliva/genes/allvall.aa.genes100.out, count=saliva.ALL.genes.count)'
mothur '#make.shared(count=saliva.ALL.genes.count, list=allvall.aa.genes100.an.unique_list.list)'
