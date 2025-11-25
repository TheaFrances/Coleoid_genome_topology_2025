# -*- coding: utf-8 -*- 
# Count intervening genes in loops and write intervening genes to the outfile, normalising forloop size.
# The script identifies genes that are located **between** the end of the first bin and the start of the second bin.
# Note this script outputs the genes in the first half of the loop to the bin1 column, and the genes in the second half of the TAD/loop to the bin2 column.
#==============================================================================

import argparse
import sys
from collections import defaultdict
import os
import statistics
import bisect

parser = argparse.ArgumentParser()
parser.add_argument("gff_bed", type=str, help="bed/gff with gene locations")
parser.add_argument("loop_table", type=str, help="tsv loop file from mustache")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def main():
    species = args.gff_bed.split("/")[-1].split(".")[0]
    print(f"Species: {species}")

    # Read loop table
    loop_table = []
    with open(args.loop_table, 'r') as infile:
        header = infile.readline().strip().split("\t")
        for line in infile:
            parts = line.strip().split("\t")
            loop_table.append({header[i]: parts[i] for i in range(len(header))})

    # Read and index gene locations
    gene_dict = defaultdict(list)
    with open(args.gff_bed, 'r') as gff:
        for line in gff:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chromosome, start, end, gene_id = parts[0], int(parts[1]), int(parts[2]), parts[3]
                gene_dict[chromosome].append((min(start, end), max(start, end), gene_id))  # Ensure start < end
    
    # Sort genes for binary search
    for chrom in gene_dict:
        gene_dict[chrom].sort()

    significant_loops_count = 0
    intervening_gene_counts = []
    loop_sizes = []
    gene_coverage_fractions = []

    outname = args.loop_table.replace(".tsv", ".tsv.intervening_genes")
    with open(outname, 'w') as outfile:
        outfile.write("chromosome\tbin1_start (bp)\tbin1_end (bp)\tbin2_start (bp)\tbin2_end (bp)\tfdr\tinner_loop_size (bp)\tintervening_gene_ids\tintervening_gene_count\ttotal_gene_bp\tgene_coverage_percentage\n")
        
        for loop in loop_table:
            try:
                fdr = float(loop['FDR'])
            except ValueError:
                continue

            if fdr < 0.05:
                significant_loops_count += 1
                chrom1 = loop['BIN1_CHR']
                bin1_start = int(loop['BIN1_START'])
                bin1_end = int(loop['BIN1_END'])
                bin2_start = int(loop['BIN2_START'])
                bin2_end = int(loop['BIN2_END'])

                # Ensure bins are in correct order (bin1 before bin2). This should always be the case for mustache output.
                if bin1_end < bin2_start:
                    region_start, region_end = bin1_end, bin2_start  # Normal case
                else:
                    region_start, region_end = bin2_start, bin1_end  # Swap if reversed

                # Ensure valid intervening region
                loop_size = region_end - region_start
                if loop_size <= 0:
                    print(f"Skipping overlapping or adjacent bins: ({bin1_start}-{bin1_end}) and ({bin2_start}-{bin2_end})")
                    continue  # Skip invalid regions

                genes = gene_dict.get(chrom1, [])
                idx1 = bisect.bisect_left(genes, (region_start,)) 
                idx2 = bisect.bisect_right(genes, (region_end,))
                # Get genes in the region and filter out genes >1.5Mb (misannotations)
                gene_intervals = [(gene[0], gene[1], gene[2]) for gene in genes[idx1:idx2] if (gene[1] - gene[0]) <= 1_500_000]    

                #Unhash the the two lines below and hash the line above if you don't want to filter out genes >1.5Mb
                #intervening_genes = [gene[2] for gene in genes[idx1:idx2] if region_start < gene[0] < region_end and region_start < gene[1] < region_end]
                #total_gene_bp = sum(abs(gene[1] - gene[0]) for gene in genes[idx1:idx2] if region_start < gene[0] < region_end and region_start < gene[1] < region_end)    

                intervening_genes = [gene[2] for gene in gene_intervals]
                # Merge overlapping genes and calculate total unique coverage
                merged_intervals = []
                for start, end, _ in sorted(gene_intervals):
                    if merged_intervals and start <= merged_intervals[-1][1]:  
                        merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], end))  
                    else:
                        merged_intervals.append((start, end))

                total_gene_bp = sum(end - start for start, end in merged_intervals)

                # Calculate gene coverage percentage
                # Prevent division by zero and cap at 100, avoiding issues with rounding up
                gene_coverage_percentage = min((total_gene_bp / loop_size) * 100, 100) if loop_size > 0 else 0

                intervening_gene_counts.append(len(intervening_genes))
                loop_sizes.append(loop_size)
                gene_coverage_fractions.append(gene_coverage_percentage)

                outfile.write(f"{chrom1}\t{loop['BIN1_START']}\t{bin1_end}\t{bin2_start}\t{loop['BIN2_END']}\t{fdr}\t{loop_size}\t" +
                              f"{','.join(intervening_genes) if intervening_genes else '0'}\t{len(intervening_genes)}\t{total_gene_bp}\t{gene_coverage_percentage:.6f}%\n")

    if significant_loops_count > 0:
        print(f"Total significant loops = {significant_loops_count}")
        print(f"Mean number of intervening genes = {statistics.mean(intervening_gene_counts):.6f}")
        print(f"Median number of intervening genes = {statistics.median(intervening_gene_counts)}")
        print(f"Mean gene coverage (percentage of loop occupied by genes) = {statistics.mean(gene_coverage_fractions):.6f}%")
        print(f"Median gene coverage (percentage of loop occupied by genes) = {statistics.median(gene_coverage_fractions):.6f}%")
    else:
        print("\nNo significant loops found.")

    print(f"Number of significant loops = {significant_loops_count}")
    print(f"Output written to {outname}:")

if __name__ == "__main__":
    main()
