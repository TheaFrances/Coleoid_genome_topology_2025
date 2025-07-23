#Add S. officinalis distances to gene pair pair distance file for E. scolopes and O. bimaculoides.
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

parser.add_argument("sof_bed", type=str,
                    help="Sepia officinalis bed file of chromosome | start | end | gene")

parser.add_argument("orth_distances", type=str,
                    help="File of: orthologous interacting gene pair | eupsc distance | octbi distance")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#==============================================================================

def main():
      
    sof_genes = {}
    for line in open(args.sof_bed, "r"):
        line = line.rstrip()
        line = line.split("\t")
        chrom = line[0]
        gene = line[3]
        gene_start = line[1]
        gene_end = line[2]
        sof_genes[gene] = [chrom, gene_start, gene_end]

    input_interactions = 0
    sof_orth_same_chrom = 0
    sof_orths_diff_chrom = 0
    counter_no_ortho_gene = 0

    with open(args.orth_distances, "r") as orth_distances:
        next(orth_distances) #Skip header
        outname = args.orth_distances.split("_octbi")[0] + "_sepof_dists_no_dups.txt"
        with open(outname, "w") as outfile:
            outfile.write("Orth_pair_based_on_eupsc_interaction_matrix\tGenomic_distance_eupsc_bp\tGenomic_distance_octbi_bp\tGenomic_distance_sepof_bp\n")
            for line in orth_distances:
                input_interactions += 1
                line = line.rstrip()
                orth_gene_pair = line.split("\t")[0]
                eupsc_gene1 = orth_gene_pair.split(",")[0].split(";")[0]
                eupsc_gene2 = orth_gene_pair.split(",")[1].split(";")[0]
                if eupsc_gene1 in sof_genes and eupsc_gene2 in sof_genes:
                    sof_chrom1 = sof_genes[eupsc_gene1][0]
                    sof_chrom2 = sof_genes[eupsc_gene2][0]
                    sof_start1 = sof_genes[eupsc_gene1][1]
                    sof_start2 = sof_genes[eupsc_gene2][1]
                    sof_end1 = sof_genes[eupsc_gene1][2]
                    sof_end2 = sof_genes[eupsc_gene2][2]
                    if sof_chrom1 == sof_chrom2:
                        sof_orth_same_chrom += 1
                        if int(sof_start1) < int(sof_start2): #If first gene comes before second gene in the DNA sequence
                            sof_dist = int(sof_start2) - int(sof_end1)
                            if int(sof_dist) < 0: #If genes overlap and sof_dist is negative, set sof_dist to 0.
                                outfile.write(f"{line}\t0\n")
                            else:
                                outfile.write(f"{line}\t{sof_dist}\n")
                        elif int(sof_start2) < int(sof_start1): #If second gene comes before first gene in the DNA sequence
                            sof_dist = int(sof_start1) - int(sof_end2)
                            if int(sof_dist) < 0: #If genes overlap and sof_dist is negative, set sof_dist to 0.
                                outfile.write(f"{line}\t0\n")
                            else:
                                outfile.write(f"{line}\t{sof_dist}\n")
                        else:
                            outfile.write(f"{line}\t0\n")
                            #print("Genes start at the same position in the Sepia officinalis DNA sequence.")
                            #print("Gene1:", eupsc_gene1, "Gene2:", eupsc_gene2)
                            #print("Gene1 start:",sof_start1, "Gene1 end:",sof_end1, "Gene2 start:",sof_start2, "Gene2 end:",sof_end2)
                    elif sof_chrom1 != sof_chrom2:
                        sof_orths_diff_chrom += 1
                elif eupsc_gene1 not in sof_genes or eupsc_gene2 not in sof_genes:
                    counter_no_ortho_gene += 1

    print("Note this script is only compatible with 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged_rm_dups.txt file. Should be no duplicate interacting gene pairs in this file. Which means there are no duplictes in the numbers below.")
    print("Number of input interactions = ", input_interactions)
    print("Number of Sepia officinalis orthologous gene pairs on same chromosome with distances included in outfile = ", sof_orth_same_chrom)
    print("Number of Sepia officinalis orthologous gene pairs on different chromosomes not written to outfile = ", sof_orths_diff_chrom)
    print("Number of orthologous gene pairs with at least one ortholog not in Sepia officinalis bed file = ", counter_no_ortho_gene)
    print("Output file written to:", outname)

if __name__ == "__main__":
    main()

