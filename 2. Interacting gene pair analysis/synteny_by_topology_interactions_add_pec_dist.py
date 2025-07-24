# -*- coding: utf-8 -*- 
#Get distances between gene pairs in P. maximus.
#Code (like all the geom dist scripts) ignores pairs of genes that are genes within genes or are the same genes. Distance is only written to outfile when they overlap or are in completely different places.
#=============================================================================
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

parser.add_argument("pmax2_orthos", type=str,
                    help="Input file of ceph species and scallop orthologs")

parser.add_argument("pmax2_bed", type=str,
                    help="pmax2 bed file of chromosome | start | end | gene")

parser.add_argument("merged_int_dist", type=str,
                    help="Merged file of eup vs obi interactions and their genomic distances and pec chromosome status. Output of plot_geom_dist_2spp_for_categories R script")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#==============================================================================

def main():

    with open(args.pmax2_orthos, "r") as orthos:
        ortho_dict = {}
        for line in orthos:
            line = line.rstrip()
            ceph_gene = line.split("\t")[0]
            pmax_gene = line.split("\t")[1]
            ortho_dict[ceph_gene] = pmax_gene

    gene_loc_pmax2 = {} #Dictionary of genes and their genomic locations
    with open(args.pmax2_bed, "r") as pmax2_bed:
        for line in pmax2_bed:
            line = line.rstrip()
            gene = line.split("\t")[3]
            chrom = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            loc =  chrom + ":" + start + ":" + end
            gene_loc_pmax2[gene] = loc

    counter = 0
    counter_all = 0
    outname = args.merged_int_dist.split(".txt")[0] + "_PEC_dist.txt"
    with open(outname, "w") as outfile:
        with open(args.merged_int_dist, "r") as interaction_file:
            for line in interaction_file:
                line = line.rstrip()
                if line.startswith("Orth"):
                    outfile.write(f"{line}\tPecten_dist\n")
                elif line.split("\t")[2] == "same_pec_chrs": #If the pecten genes are on the same chromosome, calculate distance between them. Don't need to check they're on the same chromosomes on exclude scaffolds because previous script has already cehcked this.
                    outfile.write(f"{line}\t")
                    genes_bin1 = line.split("\t")[0].split(",")[0].split(";")[0]
                    genes_bin2 = line.split("\t")[0].split(",")[1].split(";")[0]
                    if genes_bin1 in ortho_dict and genes_bin2 in ortho_dict:
                        #Look up the genes in ortho_dict and get their locations
                        pec_gene1 = ortho_dict[genes_bin1]
                        pec_gene2 = ortho_dict[genes_bin2]
                        if pec_gene1 in gene_loc_pmax2 and pec_gene2 in gene_loc_pmax2: 
                            loc1 = gene_loc_pmax2[pec_gene1]
                            start1 = loc1.split(":")[1]
                            end1 = loc1.split(":")[2]
                            loc2 = gene_loc_pmax2[pec_gene2]
                            start2 = loc2.split(":")[1]
                            end2 = loc2.split(":")[2]
                            if int(start2) > int(start1) and int(end2) > int(end1): #If first gene comes first (and is not the same gene)
                                counter +=1
                                counter_all +=1
                                genom_dist = int(start2) - int(end1) 
                                if genom_dist > 0:
                                    outfile.write(f"{genom_dist}\n")
                                if genom_dist <= 0: #If the octopus/pecten genes are overlapping, give a value of zero.
                                    outfile.write("0\n")
                            elif int(start2) < int(start1) and int(end2) < int(end1): #If second gene comes first (and is not the same gene)
                                counter +=1
                                counter_all +=1
                                genom_dist = int(start1) - int(end2)
                                if genom_dist > 0: #If the octopus/pecten genes are overlapping, give a value of zero.
                                    outfile.write(f"{genom_dist}\n")
                                if genom_dist <= 0:
                                    outfile.write("0\n")
                            else: #If the gene pair is the same gene, or the gene is within another gene, you don't want this in your final dataset, so add NA.
                                counter_all +=1
                                outfile.write("NA\n")
    
    print("Number of P. maximus distances added to file = ", counter)
    print("Number of P. maximus distances added to file including NAs (the same gene or 'within' one another) = ", counter_all)
    print(f"Output written to: {outname}")
if __name__ == "__main__":
    main()
