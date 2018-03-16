from __future__ import print_function
from sys import argv, stdout
from collections import Counter
import re

if __name__ == "__main__":
    infile_name = argv[1]
    outfile = open(stdout if len(argv) < 3 or argv[2] == '-' else argv[2], 'w')
    exon_counts = Counter()
    fbtr_finder = re.compile('FBtr[0-9]*')
    last_fbtr = ''
    with open(infile_name) as infile:
        for line in infile:
            line = line.strip()
            fbtr = fbtr_finder.findall(line)[0]
            if fbtr != last_fbtr:
                print(line.replace('exon', 'gene', 1), file=outfile)
                print(line.replace('exon', 'transcript', 1), file=outfile)
            exon_counts[fbtr] += 1
            print(line, 'exon_id "{}_{}"; '.format(fbtr, exon_counts[fbtr]),
                  file=outfile)

    outfile.close()
