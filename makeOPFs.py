'''
USAGE: ./makeOPFs genes.faa [options]
'''

import os
import sys
import argparse
import multiprocessing
from __future__ import division
import math as m
from collections import defaultdict

def which(program):
	'''from https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python'''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def get_next_fasta (fileObject):
    '''usage: for header, seq in get_next_fasta(fileObject):
    '''
    header = ''
    seq = ''

    #The following for loop gets the header of the first fasta
    #record. Skips any leading junk in the file
    for line in fileObject:
        if line.startswith('>'):
            header = line.strip()
            break

    for line in fileObject:
        if line.startswith('>'):
            yield header, seq
            header = line.strip()
            seq = ''
        else:
            seq += line.strip()
    #yield the last entry
    if header:
        yield header, seq

def make_count_file(filelist, fasta, outfile):
	countdict = defaultdict(int)
	groupdict = {}

	#list of files output from bowtie
	for file in filelist:
	    group = file.split(".")[0]
	    infile = open(file, 'r')
	    for line in infile:
	        line = line.strip()
	        line = line.split("\t")
	        gene = line[0].strip()
	        if gene != "*":
	            reads = int(line[2])
	            countdict[gene] += reads
	    groupdict[group] = countdict
	    countdict = defaultdict(int)
	    infile.close()

	groups = groupdict.keys()
	print groups

	numBases = 100

	outfile = open(outfile, 'w')

	outfile.write("Representative_Sequence\ttotal\t%s\n" % '\t'.join(groups) )

	for header, seq in get_next_fasta(fasta):
	    #header = header.split(" ")
	    #header = header[0]
	    totalreads = 0
	    groupcounts = []
	    for group in groups: #need to maintain same order from above
	        groupcounts.append(groupdict[group][header[1:]]) #list of counts for that gene (header)
	        totalreads = totalreads + groupdict[group][header[1:]] #this is the count of reads to that gene total

	    length = len(seq)
	    if totalreads == 0:
	        continue
	    totalreads = totalreads * numBases #normalize to number of bases
	    #num = int( m.ceil(totalreads/len(seq) ) )
	    groupcounts[:] = [int(m.ceil((x*numBases)/len(seq))) for x in groupcounts]
	    outfile.write("%s\t%s\t%s\n" % (header[1:], sum(groupcounts), '\t'.join(map(str, groupcounts) ) ) )


parser = argparse.ArgumentParser(description='make OPFs from a fasta file of genes')
parser.add_argument('genes', nargs='+')
parser.add_argument('--nuc_genes', default='', help='path to nucleotide genes fasta, for use with --calc_counts')
parser.add_argument('--blast', default=which('blastp'), help='path to BLAST')
parser.add_argument('--blastdb', default=which('makeblastdb'), help='path to BLAST')
parser.add_argument('--mothur', default=which('mothur'), help='path to mothur')
parser.add_argument('--nproc', default=1, help='number of processors to use')
parser.add_argument('--bowtie2', default=which('bowtie2'), help='path to bowtie2')
parser.add_argument('--bowtie2-build', default=which('bowtie2-build'), help='path to bowtie2-build')
parser.add_argument('--reads', default='', help='read file (used with bowtie2 and calc_counts option)')
parser.add_argument('--calc_counts', default=False, help='use bowtie to calculate counts')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

##make sure BLAST and mothur are installed and we can find them if they're not given to us
if args.blast == None or args.mothur == None or args.blastdb == None:
	print "Problem with path(s), can't find one or more programs:"
	print "blast: ", args.blast
	print "blastdb: ", args.blastdb
	print "mothur: ", args.mothur
	sys.exit(1)

if args.calc_counts == True and args.bowtie2 == None or args.botie2-build == None:
	print "can't find bowtie2 or bowtie2-build. Bowtie2 and bowtie2-build are required with --calc_counts"
	print "bowtie2: ", args.bowtie2
	print "bowtie2-build: ", args.bowtie2-build
	sys.exit(1)

if args.calc_counts == True and args.reads == '':
	print "--reads required with --calc_counts option"
	sys.exit(1)

##allvall blast
if args.nproc == 'all'
	args.nproc = multiprocessing.cpu_count()

else:
	args.nproc = int(args.nproc)

gene_base = args.genes.split('.')[0]



##build blast db
cmd = "{2} -in {0} -dbtype prot -out {1}.blastDB".format(args.genes, gene_base, args.blastdb)

##run blast

cmd = "{4} -query {0} -db {1}.blastDB -out {1}.blastout -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads {3}".format(args.genes, gene_base, blastout, args.nproc, args.blast)

##calculate count file
if args.calc_count:
	cmd = "{0} {1} {1}.bowtie2db".format(args.bowtie2-build, args.genes, gene_base)
	cmd = "{0} {1}.bowtie2db -f {2} -p {3} -S reads.aligned2{1}.sam".format(args.bowtie2, genes_base, args.reads, args.nproc)
	bowtieout = []
	make_count_file(bowtieout, args.nuc_genes, "{0}.count".format(gene_base))

##run mothur

cmd = "{0} '#mgcluster()'".format(args.mothur)
