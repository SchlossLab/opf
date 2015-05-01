import sys

countfile = open(sys.argv[1], 'r')
blastfile = open(sys.argv[2], 'r')
outfile = open(sys.argv[3], 'w')

counts = set()
blast = set()

for line in countfile:
    if line.startswith("R"):
        nGroups = len(line.strip().split("\t")) - 1
        continue
    line = line.strip()
    line = line.split("\t")
    gene = line[0]
    counts.add(gene)

countfile.close()

for line in blastfile:
    line = line.strip()
    line = line.split("\t")
    gene = line[0]
    blast.add(gene)
    blast.add(line[1])

blastfile.close()
zeros = [0]*nGroups
for gene in blast ^ counts:
    outfile.write("%s\t%s\n" % (gene, '\t'.join(map(str, zeros)) ) )
    #outfile.write("%s\t0\t0\t0\n" % gene )
