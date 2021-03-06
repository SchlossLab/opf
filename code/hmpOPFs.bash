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

DATASETS=(anterior_nares
buccal_mucosa
posterior_fornix
supragingival_plaque
tongue_dorsum)

#removed datasets (done previosuly)
#saliva
#stool
#attached_keratinized_gingiva
#left_retroauricular_crease

READS_LOC="ftp://public-ftp.hmpdacc.org/Illumina/"
GENES_LOC="ftp://public-ftp.hmpdacc.org/HMGI/"

BASE_DIR=pwd

function fun_blast{
	bzcat *_aa.fasta.bz2 > all.aa.genes.fasta && python ~/scripts/removeShortContigs.py all.aa.genes.fasta all.aa.genes100.fasta 100
	formatdb -i all.aa.genes100.fasta -p T -n all.aa.genes100.fasta.blastDB
	blastp -query all.aa.genes100.fasta -db all.aa.genes100.fasta.blastDB -out allvall.aa.genes100.out -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads 8

}

function fun_bowtie_build{
	bzcat *_nucleotide.fasta.bz2 > all.nucleotide.genes.fasta && python ~/scripts/removeShortContigs.py all.nucleotide.genes.fasta all.nucleotide.genes300.fasta 300

	bowtie-build all.nucleotide.genes300.fasta all.nucleotide.genes300.fasta.bowtieDB
}

function fun_bowtie{

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
}

function fun_cluster{
	DATASET=$1
	python ../../opf/code/makecountfilefrombowtie.py -f ../genes/all.aa.genes100.fasta -o saliva.genes.count -l *.txt
	getMissingGenes.py ${DATASET}.genes.count ../calledgenes/allvall.aa.genes100.out missingGenes.txt && cat ${DATASET}.genes.count missingGenes.txt > ${DATASET}.genes.ALL.count

	mothur '#mgcluster(blast=genes/allvall.aa.genes100.out, count=${DATASET}.ALL.genes.count)'
	mothur '#make.shared(count=${DATASET}.ALL.genes.count, list=allvall.aa.genes100.an.unique_list.list)'
}

#parallelize this
for i in ${DATASETS[*]};
do
	mkdir $i;
	mkdir $i/cluster;

	wget ${READS_LOC}${i}/* -P $i/rawreads;
	wget ${GENES_LOC}${i}/* -P $i/genes;

	cd $BASE_DIR/$i/genes;
	fun_blast;
	fun_bowtie_build

	cd $BASE_DIR/$i/rawreads;
	for f in *.tar.bz2;
	do
		fun_bowtie $f;
	done;

	cd $BASE_DIR/$i/cluster;
	fun_cluster;
done;
