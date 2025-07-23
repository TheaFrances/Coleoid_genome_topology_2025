# -*- coding: utf-8 -*- 
#Format output of last python script: check_orthos_dump.py
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


parser.add_argument("spp2_orthos_int", type=str,
                    help="Output file from print_genes_int_freq_same_diff.py script. Formatted like: Ortholog_pair_genes_bin1 | Ortholog_pair_genes_bin2 | Interaction_frequency")

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
#==============================================================================

def main():

    int_ortho_dict = defaultdict(list)
    with open(args.spp2_orthos_int, "r") as input_file:
        for line in input_file:
            if not line.startswith("Ortholog"):
                line = line.rstrip()
                int_freq = line.split("\t")[2]
                if not line.split("\t")[0] == "NA":
                    if not line.split("\t")[1] == "NA":
                        orthos_bin1 = set(line.split("\t")[0].split(",")) #Get only unqiue distance values because otherwise you will just end up plotting the same interaction freq + distance twice in R.
                        orthos_bin2 = set(line.split("\t")[1].split(","))
                        for i in orthos_bin1:
                            if not i == "":
                                for x in orthos_bin2:
                                    if not x == "": #This gets rid of the bit after the last comma
                                        if not i == x: #This gets rid of self interacting genes
                                            ortho_pair = i + "," + x
                                            int_ortho_dict[int_freq].append(ortho_pair)
    
    #Test dict for testing code for removing 'swapped bin duplicates'
    #     int_ortho_dict = {
    # 1: ["geneA,geneB", "geneC,geneD"],
    # 2: ["geneX,geneY", "geneZ,geneW", "geneW,geneZ"],
    # 3: ["geneX,geneY", "geneZ,geneW", "geneW,geneZ"]}

    new_int_ortho_dict = defaultdict(set)  # New defaultdict to store non-swapped pairs. Set because once you swap the bins, there will be duplicates and you don't want these.
    for int_freq, gene_pairs_list in int_ortho_dict.items():
        for gene_pairs in gene_pairs_list:
            pair1, pair2 = gene_pairs.split(',')
            if pair1 > pair2:  #check if the pairs are in alphabetical or numerical order. If pair1 is greater than pair2, it implies that the pairs are in the reverse order, and they need to be swapped to maintain a consistent order. You need this bit so you're not swapping both bins round.
                pair1, pair2 = pair2, pair1  # Swap pairs if necessary
                new_pairs_value = f"{pair1},{pair2}"
                new_int_ortho_dict[int_freq].add(new_pairs_value)
            else:
                new_int_ortho_dict[int_freq].add(gene_pairs)
    
    #Check if the swapped bin duplicates are removed            
    # print(new_int_ortho_dict)
    
    counter = 0
    outname = args.spp2_orthos_int.split(".txt")[0] +"_form.txt"
    with open(outname, "w") as outfile:
        outfile.write("Interacting_orthologous_gene_pair")
        outfile.write("\t")
        outfile.write("Interaction_frequency")
        outfile.write("\n")
        for interaction in new_int_ortho_dict: 
            orthos = new_int_ortho_dict[interaction]
            for ortho in orthos:
                counter +=1
                outfile.write(ortho)
                outfile.write("\t")
                outfile.write(interaction)
                outfile.write("\n")
    
    print("Number of interacting gene pairs in outfile (should have 'swapped bin' duplicates removed) = ", counter)
    print("Output written to:", outname)

if __name__ == "__main__":
    main()
