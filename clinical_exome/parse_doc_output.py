#!/usr/bin/python

### This will take read depth data from GATK DepthOfCoverage and will calculate Q30%
### as well as depth of coverages at varying depths.  DepthOfCoverage should be run
### twice, once with the option mbq=30 and once with mbq=0.  This will allow
### an easy calculation of these metrics per interval.
### Inputs: 2 DepthOfCoverage depth outputs and an interval list in BED format
### Run match_intervals_genes.py to create interval list with associated gene names.
### Output: tsv with chrom, interval start, interval stop, q0 depth, q30 depth,
###         and 11 columns of binned coverage values.
### John Letaw 032915 

import argparse

### Set bins globally so we don't have to pass the list over and over again.

BINS = [100, 50, 25, 10, 1, 0]

def binDepths(depth):
    """
    Input: a depth value
    Output: index of which bin to add depth
    """
    for bin in BINS:
        if int(depth) >= bin:
            return BINS.index(bin)+6

def importIntervals(handle):
    """
    Input: file handle of genomic intervals
    Output: [chrom, start, stop, gene, q0, q30, d100, d50, d25, d10, d1, d0]
    """

    print("Reading intervals file.")
    interval_list = []
    with handle as intervals:
        for intval in intervals:

            intval = intval.rstrip('\n')
            chrom = intval.split('\t')[0]
            start = int(intval.split('\t')[1])
            stop = int(intval.split('\t')[2])
            gene = intval.split('\t')[3]

            interval_list.append([chrom, start, stop, gene, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return interval_list

def importDepth(handle, up_idx, int_list):
    """
    Input: file handle of depths from DepthOfCoverage, index of where to store counts, interval list
    Output: updated interval list
    """

    print("Reading DepthOfCoverage per locus output.")
    with handle as depth:
        next(depth)  ### Skip header line.
        i = 0
        for coord in depth:

            coord = coord.rstrip('\n').split('\t')
            dcoord = int(coord[0].split(':')[1])
            ddepth = int(coord[1])

            if dcoord >= int_list[i][1] and dcoord <= int_list[i][2]:
                int_list[i][up_idx] += ddepth
                if up_idx == 4:
                    int_list[i][binDepths(ddepth)] += 1
            else:
                i += 1
                if i == len(int_list):
                    return int_list
                        
    return int_list


def sumDepths(depths):
    """
    Input: Arbitrary list of depths
    Output: sum of the depths
    """
    sum = 0
    for d in depths:
        sum += d
    return sum

def formatPercent(x, y):    
    """
    Input: Two float values.
    Output: Formatted percent string.
    """
    if y != 0.0:
        return str("{0:.1f}".format(x/y*100))
    return "0.0"
    
def writeIntervalQC(handle, int_list):
    """
    Input: file handle for writing results, interval list
    Output: written file
    """

    print("Writing Interval QC.")
    handle.write("Chromosome\tStart\tStop\tGene\tQ30%\tD100\tD50\tD25\tD10\n")
    for interval in int_list:
        handle.write(interval[0] + '\t' + str(interval[1]) + '\t' + str(interval[2]) + '\t' + interval[3] + '\t')
        handle.write(formatPercent(interval[5], interval[4]))
        for i in range(7, 11):
            handle.write('\t' + formatPercent(sumDepths(interval[6:i]), sumDepths(interval[6:])))
        handle.write('\n')
    handle.close()

def createGeneQC(int_list):
     
    gene_qc = {}

    for intval in int_list:
        if intval[3] == "NONE":
            pass
        elif intval[3] not in gene_qc:
            gene_qc[intval[3]] = intval[4:]
        elif intval[3] in gene_qc:
            for i in range(8):
                gene_qc[intval[3]][i] += intval[i+4]

    return gene_qc

def writeGeneQC(handle, gene_qc):
    
    print("Writing Gene QC.")
    handle.write("Gene\tQ30%\tD100\tD50\tD25\tD10\n")
    for gene in sorted(gene_qc):
        handle.write(gene + '\t')
        handle.write(formatPercent(gene_qc[gene][1], gene_qc[gene][0]))
        for i in range(3, 7):
            handle.write('\t' + formatPercent(sumDepths(gene_qc[gene][2:i]), sumDepths(gene_qc[gene][2:])))
        handle.write('\n')
    handle.close()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("interval_list", help="List of intervals in form chr:start-stop")
    parser.add_argument("depth_0", help="DepthOfCoverage output at mbq=0")
    parser.add_argument("depth_30", help="DepthOfCoverage output at mbq=30")
    parser.add_argument("output_int", help="Output Interval TSV")
    parser.add_argument("output_gene", help="Output Gene TSV")
    args = parser.parse_args()

    interval_list = importIntervals(open(args.interval_list, 'rU'))
    interval_list = importDepth(open(args.depth_0, 'rU'), 4, interval_list)
    interval_list = importDepth(open(args.depth_30, 'rU'), 5, interval_list)
    writeIntervalQC(open(args.output_int, 'w'), interval_list)
    gene_qc = createGeneQC(interval_list)
    writeGeneQC(open(args.output_gene, 'w'), gene_qc)

if __name__ == "__main__":
    main()
