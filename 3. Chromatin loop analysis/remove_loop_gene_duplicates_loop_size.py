# -*- coding: utf-8 -*- 
# Check loops in merged files that have the same sets of genes interacting. If they do, consider this a duplicate loop.
# Then merge these loops into one loop with the outermost coordinates of the two loops as the new coordinates.
# For loop size files.
# ==============================================================================
# Imports=======================================================================
# ==============================================================================
import argparse
import sys
from collections import defaultdict
# ==============================================================================
# Command line options==========================================================
# ==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument("loop_genes_file", type=str, help="Loop file with genes and loop sizes with suffix .genes")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
# ==============================================================================
# Main code=====================================================================
# ==============================================================================

def main():
    counter_loop_input = 0
    loop_dict = defaultdict(list)
    loop_size_list = []
    species_name = args.loop_genes_file.split("/")[-1].split("_")[0]
    
    # Read the input file and populate the loop_dict
    with open(args.loop_genes_file, "r") as loop_file:
        for line in loop_file:
            if not line.startswith("chromosome"):
                counter_loop_input += 1
                line = line.strip().split("\t")
                chrom, start, genes1, end, genes2, fdr, loop_size = line
                # Convert start and end to integers
                start, end = int(start), int(end)
                # Store all information including FDR
                loop_dict[(genes1, genes2)].append((chrom, start, end, fdr, int(loop_size)))

    counter_loop_output = 0
    output = args.loop_genes_file.split(".genes")[0] + ".genes_rm_dups"
    
    # Write the merged loops to the output file
    with open(output, "w") as outfile:
        # Write the header
        outfile.write("chromosome\tbin1_start (bp)\tgenes_bin1\tbin2_start (bp)\tgenes_bin2\tfdr\tloop_size\n")
        
        # Iterate through the loop_dict to process each unique gene pair
        for (genes_bin1, genes_bin2), values in loop_dict.items():
            if len(values) > 1:  # If there are multiple values, merge them
                chrom = values[0][0]  # All entries should have the same chromosome
                start = min(v[1] for v in values)  # Take the minimum start position
                end = max(v[2] for v in values)  # Take the maximum end position
                fdr = min(v[3] for v in values)  # Take the minimum FDR (most significant)
                loop_size = max(v[4] for v in values) # Take the maximum loop size
                loop_size_list.append(loop_size)

            else:
                # If there's only one entry, use its values directly
                chrom, start, end, fdr, loop_size = values[0]
                loop_size_list.append(loop_size)
            
            # Write the merged or single loop to the output file
            outfile.write(f"{chrom}\t{start}\t{genes_bin1}\t{end}\t{genes_bin2}\t{fdr}\t{loop_size}\n")
            counter_loop_output += 1

    print("Species = ", species_name)
    print("Number of conserved loops in input file = ", counter_loop_input)
    print("Number of loops in output with duplicates removed = ", counter_loop_output)
    print("Average loop size (with genes in) = ", sum(loop_size_list)/len(loop_size_list), "Mb")
    if counter_loop_input == counter_loop_output:
        print("All loops in input file were unique and no duplicates were removed.")
    
    print("Output written to: ", output)

if __name__ == "__main__":
    main()




                    



