'''
USAGE: ./makeOPFs genes.faa [BLAST path] [mothur path]
'''

import os
import sys
import argparse
import multiprocessing

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


parser = argparse.ArgumentParser(description='make OPFs from a fasta file of genes')
parser.add_argument('genes', nargs='+')
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

##run mothur

cmd = "{0} '#mgcluster()'".format(args.mothur)
