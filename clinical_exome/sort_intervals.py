#!/usr/bin/python

### Sort an intervals file (chrom:start-stop) by karyotypic order.
### This will be necessary to ensure we can quickly produce quality
### output based on GATK DepthOfCoverage.
### John Letaw 033015

import sys

handle_in = open(sys.argv[1], 'rU')
handle_out = open(sys.argv[2], 'w')

ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']

int_list = []

with handle_in as intervals:
    for line in intervals:
        line = line.rstrip('\n')
        int_list.append(line)
        
for chrom in ORDER:
    for entry in int_list:
        if entry.split(':')[0][3:] == str(chrom):
            handle_out.write(entry + '\n')

handle_out.close()
