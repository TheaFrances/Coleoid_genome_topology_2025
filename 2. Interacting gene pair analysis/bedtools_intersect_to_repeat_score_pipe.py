# -*- coding: utf-8 -*- 
#As bedtools intersect command is running, pipe the output into this script to get the number of each repeat type for each gene pair.     
#=============================================================================
#=============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections

def main():

  #Open nested defaultdict
    nested_dict = lambda: defaultdict(int)
    gene_repeat_counts = defaultdict(nested_dict)

    total_repeats = 0

    #Open the sys.stdin line by line.  sys.stdin.read() reads the entire file at once which is not what we want.
    for line in sys.stdin:
        line = line.rstrip().split("\t")
        gene_pair = line[3]
        repeat_type = line[7].split("Target=")[1]
        #Increment the count for the gene pair and repeat type
        gene_repeat_counts[gene_pair][repeat_type] += 1
        #Increment the total repeats counter
        total_repeats += 1


    #Print the nested dictionary so each count is on one line 
    for gene_pair, repeat_types in gene_repeat_counts.items():
        for repeat_type, count in repeat_types.items():
            print(f"{gene_pair}\t{repeat_type}\t{count}")

    #Save to file/pipe to file here
            
    #Print the total number of repeats to bottom of file
    #print(f"Total number of repeats = {total_repeats}")

if __name__ == "__main__":
    main()
