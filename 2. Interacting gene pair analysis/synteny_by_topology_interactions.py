# -*- coding: utf-8 -*- 
#Categorise orthologolous P. maximus genes as being on same/diff, same/same chromosomes for pairs of cephalopod interactions.
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

parser.add_argument("orth_interactions", type=str,
                    help="Output of eup_vs_obi_int_freq.py script, orthologous interactions and their freqs between eup and obi")


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

    gene_chrom_pmax2 = {} #Dictionary of genes and their genomic locations
    with open(args.pmax2_bed, "r") as pmax2_bed:
        for line in pmax2_bed:
            line = line.rstrip()
            gene = line.split("\t")[3]
            chrom = line.split("\t")[0]
            gene_chrom_pmax2[gene] = chrom

    counter_same = 0
    counter_diff = 0
    count_gene_not_in_bed = []
    outname = args.orth_interactions.split(".txt")[0] + "_PEC_synteny_by_topology.txt"
    with open(outname, "w") as outfile:
        outfile.write("Orth_interacting_gene_pair") #This is basaed on the euprymna interaction matrix?
        outfile.write("\t")
        outfile.write("Interaction_frequency_spp1")
        outfile.write("\t")
        outfile.write("Interaction_frequency_spp2")
        outfile.write("\t")
        outfile.write("Pecten_chromosome_status")
        outfile.write("\n")
        with open(args.orth_interactions, "r") as interaction_file:
            for line in interaction_file:
                if not line.startswith("Orth"):
                    line = line.rstrip()
                    g_bin1 = line.split("\t")[0].split(",")[0].split(";")[0]
                    g_bin2 = line.split("\t")[0].split(",")[1].split(";")[0]
                    if not g_bin1 == g_bin2:
                        if g_bin1 in ortho_dict and g_bin2 in ortho_dict:
                            pec_orth_bin1 = ortho_dict[g_bin1]
                            pec_orth_bin2 = ortho_dict[g_bin2]
                            if pec_orth_bin1 in gene_chrom_pmax2 and pec_orth_bin2 in gene_chrom_pmax2:
                                pec_chr_bin1 = gene_chrom_pmax2[pec_orth_bin1]
                                pec_chr_bin2 = gene_chrom_pmax2[pec_orth_bin2]
                                if pec_chr_bin1 == pec_chr_bin2: #No scaffolds in Pecten bed file so don't need to check this
                                    counter_same += 1
                                    outfile.write(line)
                                    outfile.write("\t")
                                    outfile.write("same_pec_chrs")
                                    outfile.write("\n")
                                if not pec_chr_bin1 == pec_chr_bin2:
                                    counter_diff += 1
                                    outfile.write(line)
                                    outfile.write("\t")
                                    outfile.write("diff_pec_chrs")
                                    outfile.write("\n")
                            if pec_orth_bin1 not in gene_chrom_pmax2:
                                count_gene_not_in_bed.append(pec_orth_bin1)
                            if pec_orth_bin2 not in gene_chrom_pmax2:
                                count_gene_not_in_bed.append(pec_orth_bin2)

    print("Number of ceph interactions on same Pecten chromosomes =", counter_same)
    print("Number of ceph interactions on different Pecten chromosomes =", counter_diff)
    print("Output written to:", outname)
    print("Note this output should not contain 'swapped bin' duplicates")

if __name__ == "__main__":
    main()
