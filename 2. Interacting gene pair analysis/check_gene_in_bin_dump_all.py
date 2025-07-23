# -*- coding: utf-8 -*- 
#This script extracts all genes within the start and end of both interacting bins. As well as genes overlapping the bins (including spanning the whole bin).
#For dumped matrices made from Juicebox .hic files.
#==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()


parser.add_argument("gff_bed", type=str,
                    help="bed/gff with gene locations")
 
parser.add_argument("int_matrix_dumped", type=str,
                    help="Sorted matrix of interaction frequency for each bin pair, dumped from juicebox. hic files")

parser.add_argument("resolution", type=int,
                    help="Resolution in bp")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#==============================================================================

def main():

    genes_bin = defaultdict(list)
    resolution = int(args.resolution)

    with open(args.gff_bed, "r") as loc_file:
        for line in loc_file:
            line = line.rstrip().split("\t")
            chromo = line[0]
            gene_start = int(line[1])                       
            gene_end = int(line[2])
            gene_bin_start = int(gene_start/resolution) #BIN ID! #This should give the anchor position of the bin in bp.
            gene_bin_end = int(gene_end/resolution)
            gene_id = line[3]
            for i in range(gene_bin_start,gene_bin_end+1): #+1 in range to include the last bin. The +1 is necessary because the range function is exclusive of the upper bound, so adding 1 ensures the last bin is included. #This will include genes in all the bins it spans.
                genes_bin[chromo+":"+str(i*resolution)].append(gene_id)

    print("Number of bins with genes in = ", len(genes_bin))
            
#If euprymna gene in bin, append to file as you loop through interaction file.
    outname = args.int_matrix_dumped.split(".txt")[0]+"_all_genes_int_freq.txt"
    with open(outname, "w") as outfile:
        outfile.write("chromosome\tbin1_start\tgenes_bin1\tbin2_start\tgenes_bin2\tint_freq\n")

        with open(args.int_matrix_dumped, "r") as interactions:
            for line in interactions:
                col = line.strip().split("\t")
                chrom = col[0]
                bin_start1 = col[1]
                bin_start2 = col[2]
                int_freq = col[3]
                bin_chr1 = chrom+":"+bin_start1
                bin_chr2 = chrom+":"+bin_start2
                if bin_chr1 in genes_bin and bin_chr2 in genes_bin: 
                    #print("both ok")
                    outfile.write(f"{chrom}\t{bin_start1}\t{','.join(genes_bin[bin_chr1])}\t{bin_start2}\t{','.join(genes_bin[bin_chr2])}\t{int_freq}\n")
        
        print("Output written to:", outname)
        
if __name__ == "__main__":
    main()
