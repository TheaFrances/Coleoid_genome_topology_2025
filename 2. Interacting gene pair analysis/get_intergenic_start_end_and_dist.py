# -*- coding: utf-8 -*- 
#Check all pairs of genes in input interaction file and output intergenic coordinates and genomic distance.
#Input files with orthologous gene pairs must be formatted for one gene pair per row and in the first column.
#=============================================================================
#=============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
#==============================================================================
# Command-line options
#==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument("spp_bed_files", type=str, nargs=3, help="Bed files of chromosome | start | end | gene")
parser.add_argument("int_ortho_file", type=str, help="File of orthologous genes and interaction frequencies")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

#==============================================================================
# Function to parse BED files
#==============================================================================
def parse_bed_file(file_path):
    species = file_path.split("/")[-1].split(".")[0]  # Extract species name from file name
    gene_locations = {}  # Store gene locations
    with open(file_path, "r") as bed_file:
        for line in bed_file:
            chrom, start, end, gene = line.strip().split("\t")[:4]
            gene_locations[gene] = (chrom, int(start), int(end))
    return species, gene_locations

#==============================================================================
# Function to compute intergenic distances (Handles species2 separately)
#==============================================================================
def get_intergenic_distances(species, gene_locations, ortho_file, is_species2=False):
    counter = 0
    outname = ortho_file.split(".txt")[0] + f"_{species}_intergenic_genom_dist.txt"
    
    with open(ortho_file, "r") as ortho_ints, open(outname, "w") as outfile:
        outfile.write(f"{species}_chromosome\tOrth_pair_based_on_eupsc_interaction_matrix\tInteraction_status\tGenomic_distance_{species}_bp\tIntergenic_start_{species}_bp\tIntergenic_end_{species}_bp\n")
        
        for line in ortho_ints:
            if line.startswith("Orth"):
                continue  # Skip header
            
            gene_pair = line.strip().split("\t")[0]
            int_status = line.strip().split("\t")[1]
            
            if is_species2:
                # Extract the second part (after the semicolon) for species2
                gene1, gene2 = gene_pair.split(",")[0].split(";")[1], gene_pair.split(",")[1].split(";")[1]
            else:
                # Extract the first part (before the semicolon) for species1 and species3
                gene1, gene2 = gene_pair.split(",")[0].split(";")[0], gene_pair.split(",")[1].split(";")[0]
            
            if gene1 in gene_locations and gene2 in gene_locations:
                chrom1, start1, end1 = gene_locations[gene1]
                chrom2, start2, end2 = gene_locations[gene2]

                if chrom1 == chrom2 and not chrom1.startswith(("scaffold", "Sc", "CAY")):
                    counter += 1
                    outfile.write(f"{chrom1}\t{gene_pair}\t{int_status}\t")
                    if start2 > start1:
                        genom_dist = start2 - end1
                        outfile.write(f"{max(genom_dist, 0)}\t{end1 if genom_dist > 0 else 'NA'}\t{start2 if genom_dist > 0 else 'NA'}\n") 
                    
                    else:
                        genom_dist = start1 - end2
                        outfile.write(f"{max(genom_dist, 0)}\t{end2 if genom_dist > 0 else 'NA'}\t{start1 if genom_dist > 0 else 'NA'}\n")

    return counter, outname

#==============================================================================
# Main function
#==============================================================================
def main():
    species_list = []
    gene_loc_dict = {}

    # Read BED files
    for file in args.spp_bed_files:
        species, gene_locations = parse_bed_file(file)
        species_list.append(species)
        gene_loc_dict[species] = gene_locations

    # Check species order
    if species_list != ["eupsc", "octbi", "sepof"]:
        print("ERROR: Order of BED files should be: eupsc, octbi, sepof")
        sys.exit(1)

    # Get genomic distances for each species
    for i, species in enumerate(species_list):
        is_species2 = (i == 1)  # Set to True only for species2
        count, output_file = get_intergenic_distances(species, gene_loc_dict[species], args.int_ortho_file, is_species2)
        print(f"Number of intrachromosomal interactions with orthologs in {species}: {count}")
        print(f"Output file: {output_file}")

if __name__ == "__main__":
    main()
