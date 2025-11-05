# -*- coding: utf-8 -*- 
# Make bed file of loop start and end coordinates (anchor points only) out of loop tsv mustache outfile.
#==============================================================================
# Imports======================================================================
#==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
#==============================================================================
# Command line options=========================================================
#==============================================================================
parser = argparse.ArgumentParser()

 
parser.add_argument("loop_table", type=str,
                    help="Loop file in mustache format")


if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

#==============================================================================
# Main code +++++++++==========================================================
#==============================================================================

def main():

    outname = args.loop_table.split(".tsv")[0] + ".bed"
    with open(args.loop_table, "r") as loops:
        with open(outname, "w") as outfile:
            for line in loops:
                if not line.startswith("BIN"):
                    col = line.strip().split("\t")
                    bin1_start = col[1]
                    bin1_end = col[2]
                    chrom2 = col[3]
                    bin2_start = col[4]
                    bin2_end = col[5]
                    loop_name = chrom1 + ":" + bin1_start + "-" + bin1_end + "_" + bin2_start + "-" + bin2_end
                    outfile.write(f"{chrom1}\t{bin1_start}\t{bin1_end}\t{loop_name}_start\n")
                    outfile.write(f"{chrom2}\t{bin2_start}\t{bin2_end}\t{loop_name}_end\n")

    print(f"Output written to: {outname}")

if __name__ == "__main__":
    main()
