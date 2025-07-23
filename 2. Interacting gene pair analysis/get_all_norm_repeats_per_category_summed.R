#This script outputs repeat files for each interaction category with columns representing repeat type, total repeat count, total genomic distance between gene pairs, normalised repeat count.

# Load necessary libraries
library(readr)
library(dplyr)
library(optparse)

# Command-line argument parsing
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input file path", metavar="character"),
  make_option(c("-d", "--dists"), type="character", default=NULL, help="E. scolopes dists file path (must be the file that has NA removed and is already averaged for multiple interaction categories per gene pair)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="top_repeats_output", help="Output file name prefix [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$dists)) {
  print_help(opt_parser)
  stop("Both repeat overlap and and 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.txt file for Euprymna distances must be supplied", call.=FALSE)
}

# Load cleaned repeat data
repeats_eupsc <- read_delim(opt$input, delim = "\t", col_names = FALSE)
colnames(repeats_eupsc) <- c("GenePair", "RepeatType", "RepeatCount")

# Load Euprymna dists data
eupsc_dists <- read_delim(opt$dists, delim = "\t")

# Merge datasets
repeats_allinfo <- merge(eupsc_dists, repeats_eupsc, by.x="Orth_pair_based_on_eupsc_interaction_matrix", by.y="GenePair")

# Remove rows with zero distances
repeats_allinfo <- repeats_allinfo %>% filter(Genomic_distance_eupsc_bp > 0)

# Define function to get summed repeats
get_summed_repeats <- function(data) {
  data %>%
    group_by(RepeatType) %>%
    summarise(TotalRepeatCount = sum(RepeatCount),
              TotalGenomicDistance = sum(Genomic_distance_eupsc_bp)) %>%
    mutate(NormCount = TotalRepeatCount / TotalGenomicDistance)
}

# Filter data by Interaction_status and save results to separate files
interaction_statuses <- unique(repeats_allinfo$Interaction_status)
print(interaction_statuses)

for (status in interaction_statuses) {
  filtered_data <- repeats_allinfo %>% filter(Interaction_status == status)
  summed_repeats <- get_summed_repeats(filtered_data)
  write.table(summed_repeats, paste0(opt$output, "_", status, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
