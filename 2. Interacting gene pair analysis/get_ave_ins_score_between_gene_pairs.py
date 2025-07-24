# -*- coding: utf-8 -*- 
# Get the mean insulation score between the gene pairs (between gene 1 end and gene 2 start). This doesn't include insulation scores in bins overlapping the gene pairs.
#=============================================================================
import argparse
import sys
import pandas as pd
#==============================================================================
# Command line options
#==============================================================================
parser = argparse.ArgumentParser(description="Calculate average insulation scores for gene pairs.")

parser.add_argument("insulation_file", type=str,
                    help="File of insulation scores in BED format")

parser.add_argument("gene_locations_file", type=str,
                    help="File of gene locations in BED format")
                    
parser.add_argument("gene_pairs_file", type=str,
                    help="File of gene pairs with interaction information")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#==============================================================================
def main():

    print("This script must be used with the gene pairs file based on the eupsc interaction matrix, with eupsc gene names before the semi colons and octbi gene names after the semi colons. Must be modified if using it another way.")
    # Determine species based on the input BED file name
    species = args.gene_locations_file.split('/')[-1].split('.bed')[0]
    window_size = args.insulation_file.split('/')[-1].split('_')[-1].split('.bed')[0]
    print(f"Species = {species}")
    print(f"Window size = {window_size}")
    if species not in ["octbi", "eupsc", "sepof"]:
        print("Species not recognized. Please provide a BED file with 'octbi' or 'eupsc' or 'sepof' in the name or modify code to add more species.")
        sys.exit(1)
    
    # Load the gene locations into a DataFrame
    gene_locations = pd.read_csv(args.gene_locations_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])

    # Load the insulation scores into a DataFrame
    insulation_scores = pd.read_csv(args.insulation_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'dot', 'score', 'dot2'])

    # Filter out scaffolds starting with "Sc"
    insulation_scores = insulation_scores[~insulation_scores['chrom'].str.startswith("Sc")]

    # Load the gene pairs into a DataFrame
    gene_pairs = pd.read_csv(args.gene_pairs_file, sep='\t')

    # Function to calculate average insulation score.
    # Note the base number of the average insulation score is caluclated as the average of the insulation scores in the bins between the gene pairs, not the distance between the two genes.
    def calculate_average_insulation_score(row, gene_locations, insulation_scores):
        genes = row['Orth_pair_based_on_eupsc_interaction_matrix'].split(',')
        
        if species == "octbi":
            gene1 = genes[0].split(';')[1]
            gene2 = genes[1].split(';')[1]
        elif species == "eupsc":
            gene1 = genes[0].split(';')[0]
            gene2 = genes[1].split(';')[0]
        elif species == "sepof":
            gene1 = genes[0].split(';')[0]
            gene2 = genes[1].split(';')[0]
        else:
            return 'nan'

        # Get the start and end positions of the gene pair
        gene1_location = gene_locations[gene_locations['gene'] == gene1]
        gene2_location = gene_locations[gene_locations['gene'] == gene2]

        if gene1_location.empty or gene2_location.empty:
            return 'nan'

        start = min(gene1_location['start'].values[0], gene2_location['start'].values[0])
        end = max(gene1_location['end'].values[0], gene2_location['end'].values[0])

        # Filter insulation scores for the region
        region_scores = insulation_scores[
            (insulation_scores['chrom'] == gene1_location['chrom'].values[0]) &
            (insulation_scores['start'] >= start) &
            (insulation_scores['end'] <= end)
        ]

        if region_scores.empty:
            return 'nan'

        # Calculate the average score, ignoring 'nan' values
        scores = pd.to_numeric(region_scores['score'], errors='coerce')
        avg_score = scores[scores.notna()].mean()

        return avg_score if pd.notna(avg_score) else 'nan' #Pandas will automatically ignore 'nan' values when calculating the mean

    # Apply the function to each row in the gene pairs DataFrame
    gene_pairs[f'Average_insulation_score_{species}'] = gene_pairs.apply(calculate_average_insulation_score, axis=1, gene_locations=gene_locations, insulation_scores=insulation_scores)

    # Save the results to a new file
    output_file = args.gene_pairs_file.split(".txt")[0] + "_" + species + "_window" + window_size + "_ins_score.txt"
    gene_pairs.to_csv(output_file, sep='\t', index=False)

    print(f"Outpu file saved to:{output_file}")

if __name__ == "__main__":
    main()
