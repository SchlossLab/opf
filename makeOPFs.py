'''
USAGE: ./makeOPFs genes.faa [BLAST path] [mothur path]
'''

import os
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
parser.add_argument('--mothur', default=which('mothur'), help='path to mothur')
parser.add_argument('--nproc', default=1, help='number of processors to use')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

##make sure BLAST and mothur are installed and we can find them if they're not given to us
if args.blast == None or args.mothur == None:
	print "Path to blast or mothur not found"
	sys.exit(1)

##allvall blast
if args.nproc == 'all'
	nproc = multiprocessing.cpu_count()

##build blast db
cmd = "makeblastdb -in {0} -dbtype prot -out {0}.blastDB".format(args.genes)

##run blast

blastout = "output"

cmd = "blastp -query {0} -db {0}.blastDB -out {1} -evalue 1e-5 -outfmt 6 -max_target_seqs 10000 -num_threads {2}".format(args.genes, blastout, args.nproc)

##run mothur
