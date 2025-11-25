# -*- coding: utf-8 -*- 
# Count intervening genes in gene pairs and write intervening genes to the outfile, this time accounting for region size.
# The script identifies genes that are located **between** the end of the first gene and the start of the second gene in the gene pair.
# This outputs all four categories of gene pairs: interacting, not_interacting, interacting_octbi_only, interacting_deca_only.
#==============================================================================
import argparse
import pandas as pd
import statistics
from intervaltree import IntervalTree

# Command line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("gff_bed", type=str, help="BED/GFF file with gene locations")
parser.add_argument("gene_pairs_file", type=str, help="TSV file with gene pairs")
parser.add_argument("interaction_status", type=str, help="Interaction status: interacting, not_interacting, or none")
args = parser.parse_args()

def preprocess_gene_locations(gff_bed):
    """Load gene locations into an interval tree for fast lookups."""
    gene_locations = {}
    chrom_intervals = {}
    
    with open(gff_bed, 'r') as gff:
        for line in gff:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom, start, end, gene_id = parts[0], int(parts[1]), int(parts[2]), parts[3]
                gene_locations[gene_id] = (chrom, start, end)
                
                if chrom not in chrom_intervals:
                    chrom_intervals[chrom] = IntervalTree()
                chrom_intervals[chrom].addi(start, end, gene_id)
    
    return gene_locations, chrom_intervals

def merge_intervals(intervals, region_start, region_end):
    """Merge overlapping intervals while ensuring they remain within (region_start, region_end)."""
    if not intervals:
        return 0
    
    # Remove duplicates
    intervals = list(dict.fromkeys(intervals))  # Removes duplicates while keeping order

    # Sort by start position
    intervals.sort()

    merged = []
    current_start, current_end = intervals[0]

    for s, e in intervals[1:]:
        if s <= current_end:  # Overlapping interval
            current_end = max(current_end, e)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = s, e

    # Append last interval
    merged.append((current_start, current_end))

    # Adjust intervals to be within the given region
    final_intervals = [(max(s, region_start), min(e, region_end)) for s, e in merged if e > region_start and s < region_end]

    # Compute total base pairs
    return sum(e - s for s, e in final_intervals)


def find_intervening_genes(gene_pairs, gene_locations, chrom_intervals, species, outname):
    """Finds genes located between gene pairs."""
    intervening_counts = []
    coverage_percentages = []
    results = []
    
    with open(outname, 'w') as outfile:
        outfile.write("gene1\tgene2\tchromosome\tgene1_start\tgene1_end\tgene2_start\tgene2_end\tintervening_genes\tintervening_count\ttotal_gene_bp\tgene_coverage_percentage\n")
        
        for pair in gene_pairs:
            genes = pair.split(",")
            if len(genes) != 2:
                continue
            
            if species == "octbi":
                gene1, gene2 = genes[0].split(";")[1], genes[1].split(";")[1]  # Extract XP IDs
            else:
                gene1, gene2 = genes[0].split(";")[0], genes[1].split(";")[0]
            
            if gene1 not in gene_locations or gene2 not in gene_locations:
                continue
            
            chrom1, start1, end1 = gene_locations[gene1]
            chrom2, start2, end2 = gene_locations[gene2]
            
            if chrom1 != chrom2:
                continue
            
            if end1 <= start2:
                start, end = end1, start2  # Normal case
            elif end2 <= start1:
                start, end = end2, start1  # Swap case
            else:
                print(f"Skipping overlapping gene pair: {gene1} ({start1}-{end1}) and {gene2} ({start2}-{end2})")
                continue

            intervening_genes = [iv.data for iv in chrom_intervals[chrom1].overlap(start, end)] # Genes in range
            #gene_intervals = [(iv.begin, iv.end) for iv in chrom_intervals[chrom1].overlap(start, end)] # Gene intervals without filtering
            # Filter out genes longer than 1.5 Mb (likely misannotations)
            gene_intervals = [(iv.begin, iv.end) for iv in chrom_intervals[chrom1].overlap(start, end) if (iv.end - iv.begin) <= 1_500_000]

            
            total_gene_bp = merge_intervals(gene_intervals, start, end) if gene_intervals else 0
            intervening_count = len(intervening_genes)
            intervening_counts.append(intervening_count)
            gene_coverage_percentage = (total_gene_bp / max(1, (end - start))) * 100  # Prevent division by zero
            gene_coverage_percentage = min(gene_coverage_percentage, 100)  # Ensure it does not exceed 100%
            coverage_percentages.append(gene_coverage_percentage)

            #Debugging
            if gene_coverage_percentage > 100:  
                print(f"Error, gene coverage percentage exceeds 100%: {gene_coverage_percentage:.2f}%")
                continue
            
            intervening_gene_list = ",".join(intervening_genes) if intervening_genes else "0"
            outfile.write(f"{gene1}\t{gene2}\t{chrom1}\t{start1}\t{end1}\t{start2}\t{end2}\t{intervening_gene_list}\t{intervening_count}\t{total_gene_bp}\t{gene_coverage_percentage:.6f}%\n")
    
    if intervening_counts:
        print(f"Total gene pairs = {len(gene_pairs)}")
        print(f"Mean intervening genes = {statistics.mean(intervening_counts):.6f}")
        print(f"Median intervening genes = {statistics.median(intervening_counts)}")
        print(f"Mean gene coverage = {statistics.mean(coverage_percentages):.6f}%")
        print(f"Median gene coverage = {statistics.median(coverage_percentages):.6f}%")
    else:
        print("No valid gene pairs found with intervening genes.")
    
    print(f"Output written to {outname}:")

def main():
    species = args.gff_bed.split("/")[-1].split(".")[0]
    print("Species =", species)
    print("Interaction status =", args.interaction_status)
    
    gene_pairs = []
    with open(args.gene_pairs_file, 'r') as infile:
        for line in infile:
            if not line.startswith("Orth"):
                parts = line.strip().split("\t")
                gene_pair = parts[0]
                int_status = parts[1]
                if args.interaction_status == "interacting_all_species":
                    if int_status == "interacting_all_species":
                        gene_pairs.append(gene_pair)
                if args.interaction_status == "interacting_octbi_only":
                    if int_status == "interacting_octbi_only":
                        gene_pairs.append(gene_pair)
                if args.interaction_status == "interacting_deca_only":
                    if int_status == "interacting_deca_only":
                        gene_pairs.append(gene_pair)
                if args.interaction_status == "not_interacting_any_species":
                    if int_status == "not_interacting_any_species":
                        gene_pairs.append(gene_pair)
    
    gene_locations, chrom_intervals = preprocess_gene_locations(args.gff_bed)
    outname = args.gene_pairs_file.replace(".txt", f".{species}_{args.interaction_status if args.interaction_status != 'none' else ''}_intervening_genes.txt")
    find_intervening_genes(gene_pairs, gene_locations, chrom_intervals, species, outname)

if __name__ == "__main__":
    main()
