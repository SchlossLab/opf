MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
	.SUFFIXES:

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

READS_LOC="ftp://public-ftp.hmpdacc.org/Illumina/"
GENES_LOC="ftp://public-ftp.hmpdacc.org/HMGI/"

BASE_DIR=pwd

get_genes:
	wget ${GENES_LOC}${i}/* -P $i/genes

get_reads:
	wget ${READS_LOC}${i}/* -P $i/rawreads

bowtie_build:
	bzcat *_nucleotide.fasta.bz2 > all.nucleotide.genes.fasta && python ~/scripts/removeShortContigs.py all.nucleotide.genes.fasta all.nucleotide.genes300.fasta 300
	bowtie-build all.nucleotide.genes300.fasta all.nucleotide.genes300.fasta.bowtieDB

bowtie_run:
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

blast_build:
	bzcat *_aa.fasta.bz2 > all.aa.genes.fasta && python ~/scripts/removeShortContigs.py all.aa.genes.fasta all.aa.genes100.fasta 100
	formatdb -i all.aa.genes100.fasta -p T -n all.aa.genes100.fasta.blastDB

blast_run:
	blastp -query all.aa.genes100.fasta -db all.aa.genes100.fasta.blastDB -out allvall.aa.genes100.out -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads 8








DATASETS = hmp pregnancy twin

HMP_READS = reads
PREG_READS = reads
TWIN_READS = reads

READS = HMP_READS PREG_READS TWIN_READS

#get reads
get_data:
	wget $(READS)

trim_reads: $(FORWARD_READS) $(REV_READS)
	seqprep -f %_R1 -r %_R2 -1 %_R1_seqprep -2 %_R2_seqprep -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

shuffle_reads:
	perl shuffleSequences_fastq.pl %_R1_seqprep %_R2_seqprep combined.seqprep.fq

#diginorm parameters
K=20
C=20
X_PRAM='1e10

diginorm:
	python $KHMER_SCRIPTS/normalize-by-median.py -k $K -C $C -x $X_PRAM --savetable $COM.trimmed.kh -p $COM.trimmed.fq
	python $KHMER_SCRIPTS/filter-abund.py -V $COM.trimmed.kh $COM.trimmed.fq.keep
	rm $COM.trimmed.kh

#assembly
assemble:
	~/bin/velveth all31 31 -fasta -short all.seqprep.fq.keep.abundfilt
	~/bin/velvetg all31 -exp_cov auto -cov_cutoff 0 -scaffolding no
	~/bin/velveth all35 35 -fasta -short all.seqprep.fq.keep.abundfilt
	~/bin/velvetg all35 -exp_cov auto -cov_cutoff 0 -scaffolding no

merge:
	cat all31/contigs.fa all35/contigs.fa > all31_35.fa
	python ~/khmer2/sandbox/multi-rename.py all.fa > all.fa.renamed

cdhit:
	cd-hit-est -M 0 -i all31_35.renamed.fa -o all31_35.renamed.fa.cdhit.out -c 0.99

minimus:
	toAmos -s all31_35.renamed.fa.cdhit.out -o all31_35.renamed.cdhit.afg
	minimus2 all31_35.renamed.cdhit
	cat all31_35.renamed.cdhit.fasta all31_35.renamed.cdhit.singletons.seq > allcontigs.fa

#files:
#all31_35.renamed.cdhit.fasta -- merged contigs from minimus2
#all31_35.renamed.cdhit.singletons.seq -- singletons from minimus2

map2human:
	blastn -query allcontigs.fa -db /share/scratch/kiverson/human_g1k_v37.fasta -out allcontigsVhuman.out -evalue 1e-5 -outfmt 6 -num_threads 16
	cut -f 1 < allcontigsVhumanTopHits.out > hits2human.txt
	python removefromblast.py hits2human.txt allcontigs.fa allcontigs_humanremoved.fa

mga:
	mga allcontigs_humanremoved.fa -m > allcontigs_humanremoved.mga.out
	python ~/scripts/getgenesbypos.py allcontigs_humanremoved.mga.out allcontigs_humanremoved.fa genesProteins.fa genesnuc.fa

allvallblast:
	blastp -query outfile4000000.fa -db ../blastdb/genesProteins2.fa -out outfile4000000vall.out -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads 16
