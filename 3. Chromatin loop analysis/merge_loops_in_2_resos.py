# -*- coding: utf-8 -*- 
#Merge loops from two resolutions, keeping the loops in the first file if they appear in both files, so put higher resolution loops first.
#Loops found in both resolutions are considered redundant if their positions fall within a defined distance (tolerance).

#The `--tolerance` parameter sets the maximum allowed difference (in base pairs) between the start and end coordinates of loops at 50 kb and 100 kb resolution for them to be considered overlapping and thus merged. 
#If a 100 kb loop is within this window of a 50 kb loop, it is discarded in favor of the higher-resolution (50 kb) call.

#This script is specifically for the output of mustache, where the columns are in a specific order.
#==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
import pandas as pd
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser(description="Merge loops from two resolutions, keeping the loops in the first file if they appear in both files.")
parser.add_argument("file_reso1", type=str, help="Path to the first and higher resolution loops file. This file will be the 'priority' list.")
parser.add_argument("file_reso2", type=str, help="Path to the second resolution loops file.")
parser.add_argument("output_file", type=str, help="Path to save the merged loops file.")
parser.add_argument("--tolerance", type=int, default=50000, help="Maximum allowed difference between loop positions to be considered overlapping.")
args = parser.parse_args()
#==============================================================================
#==============================================================================

#Function to check if two loops overlap based on a given tolerance
def loops_overlap(loop1, loop2, tolerance=args.tolerance):
    """
    This function checks if two loops overlap by comparing their coordinates 
    within a specified tolerance. The loops overlap if they have the same 
    chromosome for both bins and if the start and end positions are within 
    the tolerance for both bins.
    
    Args:
    - loop1: A loop represented as a dictionary (row from DataFrame).
    - loop2: A loop represented as a dictionary (row from DataFrame).
    - tolerance: Maximum allowed difference between loop positions to be considered overlapping.
    
    Returns:
    - True if the loops overlap, otherwise False.
    """
    return (loop1['BIN1_CHR'] == loop2['BIN1_CHR'] and
            loop1['BIN2_CHROMOSOME'] == loop2['BIN2_CHROMOSOME'] and
            abs(loop1['BIN1_START'] - loop2['BIN1_START']) <= tolerance and
            abs(loop1['BIN1_END'] - loop2['BIN1_END']) <= tolerance and
            abs(loop1['BIN2_START'] - loop2['BIN2_START']) <= tolerance and
            abs(loop1['BIN2_END'] - loop2['BIN2_END']) <= tolerance)

#Main function to process the files and merge loops
def main():
    # Load the two resolution loops into pandas DataFrames
    loops_reso1 = pd.read_csv(args.file_reso1, sep='\t')
    loops_reso2 = pd.read_csv(args.file_reso2, sep='\t')

    # Filter out loops with FDR >= 0.5 from both files (FDR is in the 7th column, index 6)
    loops_reso1 = loops_reso1[loops_reso1.iloc[:, 6] < 0.5]
    loops_reso2 = loops_reso2[loops_reso2.iloc[:, 6] < 0.5]

    merged_loops = []
    overlapping_loops_reso1 = []

    #Add all loops from the first resolution file to the merged list
    merged_loops.extend(loops_reso1.to_dict(orient='records'))

    #Check each loop in the second resolution file
    for _, loop_reso2 in loops_reso2.iterrows():
        overlapping = False
        for loop_reso1 in loops_reso1.to_dict(orient='records'):
            if loops_overlap(loop_reso2, loop_reso1, args.tolerance):
                overlapping = True
                overlapping_loops_reso1.append(loop_reso1)  #Collect overlapping loops from reso1
                break
        
        #If no overlap, add the loop with resolution 2 to the merged list
        if not overlapping:
            merged_loops.append(loop_reso2.to_dict())

    #Convert the merged loops back to a DataFrame
    merged_df = pd.DataFrame(merged_loops)

    #Save the merged loops to the output file as tsv
    merged_df.to_csv(args.output_file, sep='\t', index=False)

    #Print out the counts of loops
    print(f"Number of loops in first resolution file = {len(loops_reso1)}")
    print(f"Number of loops in second resolution file = {len(loops_reso2)}")
    print(f"Number of loops in merged file = {len(merged_df)}")

    #Print loops from reso1 that were found in both files
    #print("\nLoops in both files (from first file) with tolerance of", args.tolerance,":")
    #for loop in overlapping_loops_reso1:
    #    print("\t".join(map(str, loop.values())))

if __name__ == "__main__":
    main()
