# -*- coding: utf-8 -*- 
# Check how many loops with genes are categorised as interacting gene pairs.
# ==============================================================================
# Imports=======================================================================
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

parser.add_argument("int_file", type=str, help="Interaction pairs file")
parser.add_argument("loop_file", type=str, help="File of loops with genes in")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# ==============================================================================
# Main code=====================================================================
# ==============================================================================

def main():

    loop_species = args.loop_file.split("/")[-1].split("_")[0]
    print(f"Loop species: {loop_species}")
    interacting_gene_pair_set = set()
    non_interacting_gene_pair_set = set()
    # Process each input 
    with open(args.int_file, "r") as loops:
        for line in loops:
            if line.startswith("Orth"):
                continue
            else:
                line = line.rstrip().split("\t")
                interaction_status = line[1]
                gene_pair = line[0]
                if loop_species == "eupsc" or loop_species == "sepof":
                    if interaction_status == "interacting_deca_only" or interaction_status == "interacting_all_species":
                        gene_1 = gene_pair.split(",")[0].split(";")[0]
                        gene_2 = gene_pair.split(",")[1].split(";")[0]
                        gene_pair_one_spp = f"{gene_1}_{gene_2}"
                        interacting_gene_pair_set.add(gene_pair_one_spp)
                    elif interaction_status == "not_interacting_any_species":
                        gene_1 = gene_pair.split(",")[0].split(";")[0]
                        gene_2 = gene_pair.split(",")[1].split(";")[0]
                        gene_pair_one_spp = f"{gene_1}_{gene_2}"
                        non_interacting_gene_pair_set.add(gene_pair_one_spp)
                elif loop_species == "octbi":
                    if interaction_status == "interacting_octbi_only" or interaction_status == "interacting_all_species":
                        gene_1 = gene_pair.split(",")[0].split(";")[1]
                        gene_2 = gene_pair.split(",")[1].split(";")[1]
                        gene_pair_one_spp = f"{gene_1}_{gene_2}"
                        interacting_gene_pair_set.add(gene_pair_one_spp)
                    elif interaction_status == "not_interacting_any_species":
                        gene_1 = gene_pair.split(",")[0].split(";")[1]
                        gene_2 = gene_pair.split(",")[1].split(";")[1]
                        gene_pair_one_spp = f"{gene_1}_{gene_2}"
                        non_interacting_gene_pair_set.add(gene_pair_one_spp)
                else:
                    print("Species not recognised, script needs modifying to account for this")
                    sys.exit(1)
                
    
    print(f"Number of unique interacting gene pairs in interaction file: {len(interacting_gene_pair_set)}")
    print(f"Number of unique non-interacting gene pairs in interaction file: {len(non_interacting_gene_pair_set)}")

    count_loops_with_genes = 0
    count_int_gps_in_loops = 0
    count_non_int_gps_in_loops = 0
    gene_pairs_to_check = set()
    with open(args.loop_file, "r") as loops:
        for line in loops:
            if line.startswith("chromosome"):
                continue
            else:
                line = line.rstrip().split("\t")
                count_loops_with_genes += 1
                genes_1 = line[2]
                genes_2 = line[4]
                for gene_1 in genes_1.split(","):
                    for gene_2 in genes_2.split(","):
                        gene_pair = f"{gene_1}_{gene_2}"
                        gene_pair_rev = f"{gene_2}_{gene_1}" # Reverse the gene pair to account for the fact that the gene pair could be in the opposite order in interacting gene pairs file
                        gene_pairs_to_check.add(gene_pair)
                        gene_pairs_to_check.add(gene_pair_rev)

    for gene_pair in gene_pairs_to_check:
        if gene_pair in interacting_gene_pair_set: #Swapped bin duplicates are removed from interacting gene pairs file so you don't need to remove them when checking against loop gene pairs
          count_int_gps_in_loops += 1
    for gene_pair in non_interacting_gene_pair_set:
        if gene_pair in gene_pairs_to_check:
            count_non_int_gps_in_loops += 1                   

    print(f"Number of loops with genes = {count_loops_with_genes}")
    print(f"Number of interacting gene pairs in loop file = {count_int_gps_in_loops}")
    print(f"Percentage of interacting gene pairs in loop file = {(count_int_gps_in_loops/len(interacting_gene_pair_set))*100}")
    print(f"Number of non-interacting gene pairs in loop file = {count_non_int_gps_in_loops}")
    print(f"Percentage of non-interacting gene pairs in loop file = {(count_non_int_gps_in_loops/len(non_interacting_gene_pair_set))*100}")
    print(f"Percentage of loops in interacting gene pairs  = {(count_int_gps_in_loops/count_loops_with_genes)*100}")
    print(f"Percentage of loops in non-interacting gene pairs  = {(count_non_int_gps_in_loops/count_loops_with_genes)*100}")

    
if __name__ == "__main__":
    main()
