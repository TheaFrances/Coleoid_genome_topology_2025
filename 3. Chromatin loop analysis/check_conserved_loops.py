# -*- coding: utf-8 -*- 
# This script defines a loop as conserved if there is at least one orthologous gene pair in bin 1 and one in bin 2 (allowing reversed bin orientation across species).
# ==============================================================================
# Imports=====================================================================
# ==============================================================================
import sys
import argparse
# ==============================================================================
# Main code=====================================================================
# ==============================================================================
def main():
    # Check for correct number of arguments
    if len(sys.argv) != 5:
        print("Usage: python3  check_conserved_loops_outfile.py  <ortholog_file> <species1_loops_file> <species2_loops_file> <output_file>")
        sys.exit(1)

    # Input files
    ortholog_file = sys.argv[1]
    species1_loops_file = sys.argv[2]
    species2_loops_file = sys.argv[3]
    output_file = sys.argv[4]  # The output file

    print(f"Ortholog file: {ortholog_file}")
    print(f"Species 1 loops file: {species1_loops_file}")
    print(f"Species 2 loops file: {species2_loops_file}")
    print(f"Output file: {output_file}")
    species_1 = species1_loops_file.split("/")[-1].split("_")[0]
    species_2 = species2_loops_file.split("/")[-1].split("_")[0]
    print(f"Species 1: {species_1}")
    print(f"Species 2: {species_2}")

    # Load orthologs into a dictionary
    ortho_dict = {}
    with open(ortholog_file, "r") as orthos:
        for line in orthos:
            col = line.strip().split("\t")
            if len(col) < 2:
                continue  # Skip malformed lines
            spp1_g = col[0]
            spp2_g = col[1]
            ortho_dict[spp1_g] = spp2_g

    # Initialize a set to track unique conserved loops
    conserved_loops = set()

    # Open loop files for species 1 and 2
    with open(species1_loops_file, "r") as loops1, open(species2_loops_file, "r") as loops2:
        # Skip the headers
        next(loops1)
        next(loops2)

        # Open the output file for writing
        with open(output_file, "w") as outfile:
            # Write header to output file
            outfile.write(f"{species_1}_loop\t{species_1}_genes_bin1\t{species_1}_genes_bin2\t{species_2}_loop\t{species_2}_genes_bin1\t{species_2}_genes_bin2\n")
            
            # Process loops for species 1
            for line1 in loops1:
                col1 = line1.strip().split("\t")
                if len(col1) < 5:
                    continue  # Skip malformed lines
                chrom1 = col1[0]
                bin1_start1 = col1[1]
                bin1_genes1 = col1[2].split(", ")
                bin2_start1 = col1[3]
                bin2_genes1 = col1[4].split(", ")

                # Get orthologous genes for species 1 bins
                ortho_bin1 = {ortho_dict[gene] for gene in bin1_genes1 if gene in ortho_dict}
                ortho_bin2 = {ortho_dict[gene] for gene in bin2_genes1 if gene in ortho_dict}

                # Reset loops2 to the start for comparing every loop in species 1 with all loops in species 2
                loops2.seek(0)
                next(loops2)  # Skip the header again

                # Process loops for species 2
                for line2 in loops2:
                    col2 = line2.strip().split("\t")
                    if len(col2) < 5:
                        continue  # Skip malformed lines
                    chrom2 = col2[0]
                    bin1_start2 = col2[1]
                    bin1_genes2 = col2[2].split(", ")
                    bin2_start2 = col2[3]
                    bin2_genes2 = col2[4].split(", ")

                    # Check for conserved loops
                    if ((ortho_bin1.intersection(bin1_genes2) and ortho_bin2.intersection(bin2_genes2)) or 
                        (ortho_bin1.intersection(bin2_genes2) and ortho_bin2.intersection(bin1_genes2))):
                        
                        # Create a unique identifier for the loop
                        loop_id = f"{chrom1}:{bin1_start1}-{bin2_start1}::{chrom2}:{bin1_start2}-{bin2_start2}"
                        if loop_id not in conserved_loops:
                            conserved_loops.add(loop_id)

                            # Print to terminal
                            print(f"Conserved loop found between:")
                            print(f"  Species 1: {chrom1}:{bin1_start1}-{bin2_start1}")
                            print(f"    Genes in bin 1: {', '.join(bin1_genes1)}")
                            print(f"    Genes in bin 2: {', '.join(bin2_genes1)}")
                            print(f"  Species 2: {chrom2}:{bin1_start2}-{bin2_start2}")
                            print(f"    Genes in bin 1: {', '.join(bin1_genes2)}")
                            print(f"    Genes in bin 2: {', '.join(bin2_genes2)}")
                            
                            # Write to output file
                            outfile.write(f"{chrom1}:{bin1_start1}-{bin2_start1}\t{', '.join(bin1_genes1)}\t{', '.join(bin2_genes1)}\t{chrom2}:{bin1_start2}-{bin2_start2}\t{', '.join(bin1_genes2)}\t{', '.join(bin2_genes2)}\n")
    
    print(f"Total number of conserved loops = {len(conserved_loops)}")
    print(f"Conserved loops written to: {output_file}")
  
if __name__ == "__main__":
    main()
