#Add Sepia officinalis interactions to output of synteny_by_topology_interactions.py script.
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

parser.add_argument("sof_interactions", type=str,
                    help="Input file of: chromosome	bin1_start	genes_bin1	bin2_start	genes_bin2	int_freqs in Sepia officinalis, output of check_gene_in_bin_dump_all.py script ran of sof dumped matrix.")


parser.add_argument("orth_interactions", type=str,
                    help="Output of MACIs_by_topology_interactions.py script, orthologous interactions and their freqs between eup and obi and P. max chromosome status")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#==============================================================================

def main():


    sof_dict = defaultdict(list)  
    with open(args.sof_interactions, "r") as sof_file:
        next(sof_file) #Skip header
        for line in sof_file:
            line = line.rstrip().split("\t")
            genes_bin1 = line[2].split(",")
            genes_bin2 = line[4].split(",") 
            int_freq = line[5]
            for gene1 in genes_bin1:
                for gene2 in genes_bin2:
                    sof_dict[gene1 + "_" + gene2].append(int_freq)
                    #print("Added", gene1, gene2, "to sof_dict with int_freq", int_freq)
    
    print("Number of Sepia officinalis interactions (note contains 'swapped bin' duplicates) = ", len(sof_dict))

    counter_orth_gene_pairs = 0
    counter_interactions = 0
    with open(args.orth_interactions, "r") as orth_file:
        outname = args.orth_interactions.split("PEC_MACIs_by_topology.txt")[0] + "PEC_chr_status_with_sof.txt"
        with open(outname, "w") as outfile:
            outfile.write("Orth_interacting_gene_pair\tInteraction_frequency_spp1\tInteraction_frequency_spp2\tPecten_chromosome_status\tInteraction_frequency_Sepia\n")
            next(orth_file) #Skip header
            for line in orth_file:
                line = line.rstrip()
                eupsc_gene1 = line.split("\t")[0].split(",")[0].split(";")[0]
                eupsc_gene2 = line.split("\t")[0].split(",")[1].split(";")[0]
                eupsc_genes = eupsc_gene1 + "_" + eupsc_gene2
                if eupsc_genes in sof_dict:
                   # print("Found", eupsc_genes, "in sof_dict")
                    counter_orth_gene_pairs += 1
                    sof_int_freq = sof_dict[eupsc_genes]
                    for i in sof_int_freq:
                        outfile.write(f"{line}\t{i}\n")
                        counter_interactions += 1

    print("Number of orthologous gene pairs in input file that also have Sepia officinalis orthologs  = ", counter_orth_gene_pairs)
    print("Number of interactions appended to output file = ", counter_interactions)
    print("NOTE both above numbers have duplicate interactions with different interaction frequencies as in the two input files, but swapped bin duplicates were already removed previously from the second input file, and so will not be in the output file")
    print("Output file written to:", outname)

if __name__ == "__main__":
    main()
