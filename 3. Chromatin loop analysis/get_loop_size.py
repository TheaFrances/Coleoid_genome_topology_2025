# -*- coding: utf-8 -*- 
# Add a column of loop size (end of end bin minus start of first bin) to mustache output file.
# ==============================================================================
# Imports=====================================================================
# ==============================================================================
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

parser.add_argument("loop_file", type=str,
                    help="File of loops formatted like the output of mustache")
 
if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()
# ==============================================================================
# Main code=====================================================================
# ==============================================================================

def main():

    if args.loop_file.endswith(".tsv"):
        output_file = args.loop_file.split(".tsv")[0] + ".loopsize.tsv"
        loop_list = []
        with open(args.loop_file, 'r') as infile:
            with open(output_file, 'w') as outfile:
                for line in infile:
                    line = line.rstrip()
                    fdr = line.split('\t')[6]
                    if line.startswith('BIN'):
                        outfile.write(f"{line}\tLOOP_SIZE\n")
                    elif float(fdr) < 0.5:
                        bin1_start = int(line.split('\t')[1])
                        bin2_end = int(line.split('\t')[5])
                        #Calculate the loop size
                        loop_size = bin2_end - bin1_start
                        outfile.write(f"{line.rstrip()}\t{loop_size}\n")
                        loop_list.append(loop_size)

    
    print("Total number of significant, differential loops = ", len(loop_list))
    print("Average loop size = ", sum(loop_list)/len(loop_list), "Mb")
    print("Output written to: ", output_file)

if __name__ == "__main__":
    main()
