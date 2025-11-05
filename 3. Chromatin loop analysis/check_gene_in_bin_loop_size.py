# -*- coding: utf-8 -*- 
# Get the genes in loop anchors from file in mustache format.
# This script takes all genes within the start and end of both differential loop bins. As well as genes overlapping the bins (including spanning the whole bin).
# It also adds a column of loop size.
# ==============================================================================
# Main code=====================================================================
# ==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
# ==============================================================================
# Command line options==========================================================
# ==============================================================================
parser = argparse.ArgumentParser()


parser.add_argument("gff_bed", type=str,
                    help="bed/gff with gene locations")
 
parser.add_argument("diff_loop_table", type=str,
                    help="Loop file in mustache format")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

# ==============================================================================
# ==============================================================================

def main():

    chrom_bin_locs = list()

    with open(args.diff_loop_table, "r") as diff_loops:
        for line in diff_loops:
            col = line.strip().split("\t") #split already here, so no need to split for individual columns
            chrom1 = col[0]
            bin1_start = col[1]
            bin1_end = col[2]
            chrom2 = col[3]
            bin2_start = col[4]
            bin2_end = col[5]
            fdr = col[6]
            if chrom1 == chrom2 and float(fdr) < 0.5: 
                chrom_bin_locs.append(chrom1+":"+bin1_start + "-" + bin1_end + "_" +  bin2_start + "-" + bin2_end) #Add instead of append to make it a set. This assigns bin id to the chromosomal "base" location (so chr1 40000 or chr5 80000 etc).
    
    print("Number of significant loops = ", len(chrom_bin_locs))

    genes_bin1 = defaultdict(list)
    genes_bin2 = defaultdict(list)
    with open(args.gff_bed, "r") as loc_file:
        for line in loc_file:
            line = line.rstrip().split("\t")
            chromo = line[0]               
            gene_bin_start = int(line[1]) #These are in actual bp in mustache, so no need to convert to bin id.
            gene_bin_end = int(line[2])
            gene_id = line[3]
            for i in chrom_bin_locs:
                chrom1 = i.split(":")[0]
                bin_start1 = int(i.split(":")[1].split("_")[0].split("-")[0])
                bin_end1 = int(i.split(":")[1].split("_")[0].split("-")[1])
                bin_start2 = int(i.split(":")[1].split("_")[1].split("-")[0])
                bin_end2 = int(i.split(":")[1].split("_")[1].split("-")[1])
                if chromo == chrom1 and ((gene_bin_start >= bin_start1 and gene_bin_start <= bin_end1) or 
                                        (gene_bin_end >= bin_start1 and gene_bin_end <= bin_end1) or
                                        (gene_bin_start <= bin_start1 and gene_bin_end >= bin_end1)):
                    genes_bin1[i].append(gene_id) #Add that gene to that bin if it falls within bin start and end or if the gene overlaps with the bin.
                if chromo == chrom1 and ((gene_bin_start >= bin_start2 and gene_bin_start <= bin_end2) or 
                                        (gene_bin_end >= bin_start2 and gene_bin_end <= bin_end2) or
                                        (gene_bin_start <= bin_start2 and gene_bin_end >= bin_end2)):
                    genes_bin2[i].append(gene_id) #Add that gene to that bin if it falls within bin start and end or if the gene overlaps with the bin.
                  
                    #print("added: "+i +" -> "+gene_id)

    print("Number of loop bins with genes in (loop start) = ", len(genes_bin1))
    print("Number of loop bins with genes in (loop end) = ", len(genes_bin2))
    print("Number of loop bins with genes in both bins = ", len(set(genes_bin1.keys()).intersection(genes_bin2.keys())))
            
   #Write to the output file
    outname = args.diff_loop_table + ".genes"
    with open(outname, "w") as outfile:
        # Write the header
        outfile.write("chromosome\tbin1_start (bp)\tgenes_bin1\tbin2_start (bp)\tgenes_bin2\tfdr\tloop_size\n")
        
        #Read the differential loop table again, could be done in loop above, but this is more readable
        with open(args.diff_loop_table, "r") as loops_diff:
            for line in loops_diff:
                col = line.rstrip().split("\t")
                chrom1 = col[0]
                bin1_start = col[1]
                bin1_end = col[2]
                chrom2 = col[3]
                bin2_start = col[4]
                bin2_end = col[5]
                fdr = col[6]
                loop_size = col[8]
                bin_chr1 = chrom1 + ":" + bin1_start + "-" + bin1_end + "_" + bin2_start + "-" + bin2_end
                bin_chr2 = chrom2 + ":" + bin1_start + "-" + bin1_end + "_" + bin2_start + "-" + bin2_end
                if bin_chr1 in genes_bin1 and bin_chr2 in genes_bin2:
                    # Convert lists to comma-separated strings
                    genes_bin1_str = ", ".join(genes_bin1[bin_chr1])
                    genes_bin2_str = ", ".join(genes_bin2[bin_chr2])
                    # Write the data to the output file
                    outfile.write(f"{chrom1}\t{bin1_start}\t{genes_bin1_str}\t{bin2_start}\t{genes_bin2_str}\t{fdr}\t{loop_size}\n")
    print("Output written to: ", outname)

if __name__ == "__main__":
    main()
