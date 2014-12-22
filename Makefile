MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.SHELLFLAGS := -eu -o pipefail
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
	.SUFFIXES:

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

	
