#!/usr/bin/python
'''Usage: python makecountfilefrombowtie.py --type p --fasta genes.faa --out mydata.count --filelist *.txt

Ths script takes in a list of genes and output from bowtie and samtools to create a count file.
Groups are named the same as the input files. For example, if given the input files group1.txt, group2.txt the resulting groups will be named group1 and group2.


'''
from __future__ import division
import sys
import os
import math as m
from collections import defaultdict
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Make countfile from bowtie/sametools output')
parser.add_argument('-t', '--type', dest ='fasta_type', default='P')
parser.add_argument('-f', '--fasta', dest = 'fasta', required = True, type=argparse.FileType('r'))
parser.add_argument('-o', '--out', dest = 'outfile', required = False, type=argparse.FileType('w'), nargs = '?', default='data.count')
parser.add_argument('-l', '--filelist', dest = 'filelist', required = True, nargs = '+')
args = parser.parse_args()

args.fasta_type = args.fasta_type.upper()

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

countdict = defaultdict(int)
groupdict = {}

for file in args.filelist:
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

if args.fasta_type == 'P':
    numBases = 100
elif args.fasta_type == 'N':
    numBases = 300
else:
    sys.exit(1)

args.outfile.write("Representative_Sequence\ttotal\t%s\n" % '\t'.join(groups) )

for header, seq in get_next_fasta(args.fasta):
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
    args.outfile.write("%s\t%s\t%s\n" % (header[1:], sum(groupcounts), '\t'.join(map(str, groupcounts) ) ) )
