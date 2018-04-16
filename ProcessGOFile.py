from sys import argv
from collections import defaultdict

if __name__ == "__main__":
    goterms = defaultdict(list)
    pfcs = {}

    for line in open(argv[1]):
        if line.startswith('!'): continue
        if not line.strip(): continue

        data = line.split('\t')
        gene_id = data[1]
        gene_name = data[2]
        goid = data[4]
        goterms[goid].append(gene_name)
        pfcs[goid] = data[8]

    for goid in goterms:
        print(goid, pfcs[goid], *goterms[goid], sep='\t')


