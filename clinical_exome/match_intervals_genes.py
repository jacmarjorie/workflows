#!/usr/bin/python

### Takes a GTF and an interval list, such as the Agilent Clinical Exome probe set, and
### tries to place genes names with each interval.  If it can't, the gene is assigned
### as None.  No gene level QC stats will be written when the gene is assigned as
### None.
### Usage:
### python match_intervals_genes.py <gtf> <interval file> <output bed>
### NOTE: interval_file is in format <chrom:start-stop>

import sys


def importGTF(handle):

    gtf_dict = {}
    with handle as gtf:
        for line in gtf:
            if line[0] != "#":
                if line.split('\t')[2] == "gene":
                    gene_name = line.split(';')[1].split('\"')[1]
                    chrom = line.split('\t')[0]
                    start = int(line.split('\t')[3])
                    stop = int(line.split('\t')[4])
                    gtf_dict[gene_name] = [chrom, start, stop]

    return gtf_dict


def importIntervals(handle):

    interval_list = []
    with handle as intervals:
        for intval in intervals:

            intval = intval.rstrip('\n')
            chrom = intval.split(':')[0]
            start = int(intval.split(':')[1].split('-')[0])
            stop = int(intval.split(':')[1].split('-')[1])

            interval_list.append([chrom, start, stop])
    return interval_list


def matchCoords(interval_list, gtf_dict, handle):

    i = 0

    for intval in interval_list:
        for key in gtf_dict:
            if intval[0] == gtf_dict[key][0]:
                if (intval[1] >= gtf_dict[key][1] and intval[1] <= gtf_dict[key][2]) or (intval[2] <= gtf_dict[key][2] and intval[2] >= gtf_dict[key][1]):
                    intval.append(key)
                    break

        if len(intval) == 4:
            handle.write(str(intval[0]) + '\t' + str(intval[1]) + '\t' + str(intval[2]) + '\t' + str(intval[3]) + '\n')
        elif len(intval) == 3:
            handle.write(str(intval[0]) + '\t' + str(intval[1]) + '\t' + str(intval[2]) + '\tNONE\n')
        else:
            raise Exception("Error, interval does not contain required chrom, start, and stop.")

        i += 1
        if i % 1000 == 0:
            print("Processed " + str(i) + " intervals.")


    handle.close()


def main():

    handle_gtf = open(sys.argv[1], 'rU')
    handle_ivals = open(sys.argv[2], 'rU')
    
    gtf_dict = importGTF(handle_gtf)
    interval_list = importIntervals(handle_ivals)

    handle_out = open(sys.argv[3], 'w')
    matchCoords(interval_list, gtf_dict, handle_out)


if __name__ == "__main__":
    main()
