import sys

countfile = open(sys.argv[1], 'r')
blastfile = open(sys.argv[2], 'r')
outfile = open(sys.argv[3], 'w')

counts = []
blast = []

for line in countfile:
    if line.startswith("R"):
        continue
    line = line.strip()
    line = line.split("\t")
    gene = line[0]
    counts.append(gene)

countfile.close()

for line in blastfile:
    line = line.strip()
    line = line.split("\t")
    gene = line[0]
    blast.append(gene)
    blast.append(line[1])

blastfile.close()
zeros = [0]*21
for gene in set(blast) ^ set(counts):
    outfile.write("%s\t%s\n" % (gene, '\t'.join(map(str, zeros)) ) )
    #outfile.write("%s\t0\t0\t0\n" % gene )
