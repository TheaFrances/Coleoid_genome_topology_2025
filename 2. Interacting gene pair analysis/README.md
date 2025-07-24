# Interacting gene pair analysis

## Contents

- [Prepare files of gene pairs, orthologs, and extract their ancestral chromosome status](#prepare-files-of-gene-pairs-orthologs-and-extract-their-ancestral-chromosome-status)
  - [Make orthology annotation file for *S. officinalis*](#make-orthology-annotation-file-for-s-officinalis)
  - [Match gene pairs to bins](#match-gene-pairs-to-bins)
  - [Check for orthologous genes](#check-for-orthologous-genes)
  - [Format interaction matrix files](#format-interaction-matrix-files)
  - [Merge interaction frequency values across species](#merge-interaction-frequency-values-across-species)
  - [Classify interactions by their *P. maximus* chromosome status](#classify-interactions-by-their-p-maximus-chromosome-status)
  - [Add interaction frequency from *S. officinalis*](#add-interaction-frequency-from-s-officinalis)
- [Get genomic distance between gene pairs for each species](#get-genomic-distance-between-gene-pairs-for-each-species)
  - [Get genomic distances for *E. scolopes* and *O. bimaculoides* orthologous interactions](#get-genomic-distances-for-e-scolopes-and-o-bimaculoides-orthologous-interactions)
  - [Sort *E. scolopes* and *O. bimaculoides* genomic distance output files and add headers](#sort-e-scolopes-and-o-bimaculoides-genomic-distance-output-files-and-add-headers)
  - [Merge sorted genomic distance files](#merge-sorted-genomic-distance-files)
  - [Remove duplicate gene pairs](#remove-duplicate-gene-pairs)
  - [Add *S. officinalis* genomic distances to the merged file](#add-s-officinalis-genomic-distances-to-the-merged-file)
- [Plot scatterplots to define the interaction frequency threshold for interacting gene pairs](#plot-scatterplots-to-define-the-interaction-frequency-threshold-for-interacting-gene-pairs)
- [Plot boxplots of genomic distances and barplots summarising interaction and ancestral chromosomal origin categories for gene pairs](#plot-boxplots-of-genomic-distances-and-barplots-summarising-interaction-and-ancestral-chromosomal-origin-categories-for-gene-pairs)
- [Repeat association analyses of interacting gene pairs](#repeat-association-analyses-of-interacting-gene-pairs)
  - [Run RepeatModeler and RepeatMasker \[TBC\]](#run-repeatmodeler-and-repeatmasker-tbc)
  - [Get intergenic start and end positions for gene pairs](#get-intergenic-start-and-end-positions-for-gene-pairs)
  - [Convert output files to BED format and sort](#convert-output-files-to-bed-format-and-sort)
  - [Convert RepeatMasker output to GFF and then BED format](#convert-repeatmasker-output-to-gff-and-then-bed-format)
  - [Run bedtools intersect with a processing script](#run-bedtools-intersect-with-a-processing-script)
  - [Get repeat content summaries per interaction category](#get-repeat-content-summaries-per-interaction-category)
  - [Get barplots of simple repeats and RND elements \[TBC\]](#Get-barplots-of-simple-repeats-and-RND-elements)
- [Calculate co-expression of gene pairs across different interaction categories \[TBC\]](#Calculate-co-expression-of-gene-pairs-across-different-interaction-categories)


This folder documents the interacting gene pair analyses. Initital steps are demonstrated using only the *E. scolopes* (stage 29) sample 403493 at 100 kb resolution, which is later merged with the *O. bimaculoides* interaction matrix at 50 kb resolution and the *S. officinalis* interaction matrix at 100 kb resolution. Boxplots of genomic distances are also only demonstrated using *E. scolopes* distance, but based on this merged interaction matrix. For the species *S. officinalis*, no gene annotation was available for the *S. officinalis* reference genome at the time of writing this paper. Therefore, additional commands are provided at the start to classify orthologous genes as well as for some downstream analyses.

## Prepare files of gene pairs, orthologs, and extract their ancestral chromosome status

**Notes on input files of orthologous genes:**

BLASTP (v2.16.0) was run with the parameters: -evalue 1E-2 -max_target_seqs 1 -outfmt 6 to classify reciprocal best hit orthologs between *E. scolopes*, O. bimaculoides* and *P. maximus*. *E. scolopes* was used as the reference (database) in all cases, except for identifying *P. maximus* orthologs of *O. bimaculoides* genes, where *O. bimaculoides* was used as the database. Text files of 1:1 orthologs in both directions were generated: **EUPgeneOBI.txt**, **EUPgenePEC.txt**, **OBIgeneEUP.txt**, and **OBIgenePEC.txt**, where: EUP = *E. scolopes*, OBI = *O. bimaculoides*, PEC = *P. maximus*. If the species abbreviation appears first in the filename, this indicates that the species' genes are in the first column of the file. For example, in EUPgeneOBI.txt, *E. scolopes* genes are in the first column, and their best reciprocal orthologs in *O. bimaculoides* are in the second. This directionality matters in the downstream analyses.

**Note on species bed files:**

Species BED files should contain four columns: chromosome, gene start, gene end, and gene name (in that order), with no additional features. 
Here, **eupsc.bed** is used as an example, corresponding to *E. scolopes*. 

### Make orthology annotation file for *S. officinalis*

*E. scolopes* proteins were mapped to the *S. officinalis* genome using miniprot. The command used creates a miniprot index at the same time. 
This step is documented in the [`index_map_miniprot.sh`](index_map_miniprot.sh) script.

The result of the miniprot is a mixture of gff and paf entries, so we filtered to keep just gff entries:

```bash
grep -v "^##PAF" sanger_sepof_eup_prot.gff > sanger_sepof_eup_prot_clean.gff
```
Extract only highest ranks hits:

```bash
grep -E 'mRNA' sanger_sepof_eup_prot_clean.gff | grep "Rank=1;" > sanger_sepof_eup_prot_best_hits.gff

```
Remove scaffolds:

```bash
grep OZ sanger_sepof_eup_prot_best_hits.gff > sanger_sepof_eup_prot_best_hits_rm_scaff.gff
```

Count number of genes left in file:

```bash
wc -l  sanger_sepof_eup_prot_best_hits_rm_scaff.gff
24431  sanger_sepof_eup_prot_best_hits_rm_scaff.gff # There are 24431 E. scolopes proteins that mapped to the S. officinalis genome that we considered orthologs.
```
The resulting **sanger_sepof_eup_prot_best_hits_rm_scaff.gff** file was then converted into bed format into a file called **sepof.bed**, with chromosome names being from the *S. officinalis* genome file, and gene names being identical to the *E. scolopes* gene names that were the best hits mapped to the *S. officinalis* genome. 

Orthology files were also made for downstream analyses: **EUPgeneSOF.txt**, **SOFgeneEUP.txt**, although gene names were identical in both columns, the files were formatted this way to be compatible with downstream scripts.


### Match gene pairs to bins

We used a custom Python script  [`check_gene_in_bin_dump_all.py`](check_gene_in_bin_dump_all.py) to identify gene pairs in bins pairs in dumped matrices. The input file is **409493_intrachrom_allchrs_KR_100000.dumped.hic.txt** is from the end of step 1 which creates a dumped matrix of all intrachromosomal interactions. This script includes any genes overlapping with the bins in the output files:

```bash
python3 check_gene_in_bin_dump_all.py eupsc.bed 409493_intrachrom_allchrs_KR_100000.dumped.hic.txt 100000
```
Where the species BED file matches the species whose interaction matrix it is. This produced an output like:

```bash
Number of bins with genes in = 26197
Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
```

This command was repeated for *S. officinalis* at 100kb resolution and *O. bimaculoides* at 50 kb resolution.

### Check for orthologous genes

The script [`check_orthos_dump.py`](check_orthos_dump.py) was then used to assess whether gene pairs in bin interactions have orthologs in other species. The species whose orthologs appear in the first column of the ortholog file (and whose abbreviation appears first in the filename) must correspond to the species from which the interaction matrix was generated.

```bash
python3 check_orthos_dump.py EUPgeneOBI.txt 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
# Example output:
# Number of reciprocal best hit orthologs for EUP and OBI = 12002 # Where EUP is E. scolopes and OBI is O. bimaculoides
# Number of EUP interactions with at least one ortholog in OBI = 1569457 # Note this file likely contains duplicate entries that have not yet been removed, as well as multiple interacting gene pairs per line. As such, the reported number is not biologically meaningful at this stage
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos.txt
```

Commands were repeated for the following species comparisons:

*E. scolopes* ↔ *O. bimaculoides*

*E. scolopes* ↔ *P. maximus*

*O. bimaculoides* ↔ *E. scolopes*

*O. bimaculoides* ↔ *P. maximus*

**Note:** Since *S. officinalis*–*E. scolopes* orthologs were named identically to *E. scolopes* genes, it was not necessary to extract ortholog names separately in this case. *S. officinalis* interactions were added downstream.


### Format interaction matrix files

Before comparing interaction frequencies across species, run the [`format_interactions.py`](format_interactions.py) script to:
- Remove NA values
- Standardise pairwise gene interactions so that each row contains exactly one pair of interacting genes (i.e. some bins contain many genes, so all pairwise combinations are extracted)
- Remove swapped-bin duplicates (i.e. if bin1-bin2 and bin2-bin1 are both present, keep only one)

```bash
# Format E. scolopes (100 kb) with O. bimaculoides orthologs
python3 format_interactions.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos.txt
# Example output: 
# 2,565,954 interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos_form.txt

# Format E. scolopes (100 kb) with P. maximus orthologs
python3 format_interactions.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_PECorthos.txt
# Example output: 
# 2,426,852 interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_PECorthos_form.txt
```

### Merge interaction frequency values across species

Use [`eup_vs_obi_int_freq_form.py`](eup_vs_obi_int_freq_form.py) to merge gene-pair interaction frequencies between species, based on *E. scolopes* bin pairs at 100 kb resolution, handling cases where gene order differs between files.

```bash
python3 eup_vs_obi_int_freq_form.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos_form.txt \
                                    212489_intrachrom_allchrs_KR_50000.dumped.hic_all_genes_int_freq_EUPorthos_form.txt \
                                    EUPgeneOBI.txt
# Example output: 
# 3,185,945 matched gene-pair interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq.txt
```
**Note:** the resulting file still contains duplicate gene pairs with differing interaction frequencies, along with all possible gene pair combinations across the two species. This contributes to the large number of interactions and file size. These redundancies will be resolved in later steps.


### Classify interactions by their ***P. maximus*** chromosome status

This step uses [`syteny_by_topology_interactions.py`](syteny_by_topology_interactions.py) to classify gene-pair interactions as:
- Gene pair on the same *P. maximus* chromosome
- Gene pair on different *P. maximus* chromosomes

```bash
python3 synteny_by_topology_interactions.py EUPgenePEC.txt pmax2.bed \
       409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq.txt
# Example output: 
# Interactions on same P. maximus chromosomes: 917,352
# Interactions on different P. maximus chromosomes: 1,493,463
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_MACIs_by_topology.txt
```

### Add interaction frequency from *S. officinalis*

This step can be done without an ortholog file, since *S. officinalis* orthologs were named the same as *E. scolopes* genes. The script [`add_sof_interactions.py`](add_sof_interactions.py) was used.

```bash
python3 add_sof_interactions.py 992270_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt \
      409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_synteny_by_topology.txt
# Example output: 
# S. officinalis interactions: 4,355,779
# Gene pairs with S. officinalis orthologs: 1,056,029
# Appended interactions: 4,772,919
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_chr_status_with_sof.txt
```
**Note**:  Swapped-bin duplicates were already removed from the second input file and so will not appear in the output of this command. However, the resulting file still contains duplicate gene pairs with differing interaction frequencies, along with all possible gene pair combinations across the three species. This contributes to the large file size. These redundancies will be resolved in later steps.

## Get genomic distance between gene pairs for each species


### Get genomic distances for ***E. scolopes*** and ***O. bimaculoides*** orthologous interactions

The script [`eup_vs_obi_genom_dist_form.py `](eup_vs_obi_genom_dist_form.py) was used to generate one output file per species, each containing the genomic distances between gene pairs in that species.

```bash
python3 eup_vs_obi_genom_dist_form.py eupsc.bed octbi.bed 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos_form.txt
# Example output:
# Number of E. scolopes interactions with genomic distances: 2,507828
# Number of O. bimaculoides interactions with genomic distances: 1,307932
# Output written to: 409493_intrachrom_allchrs_KR_100000_eupsc_genom_dist_orthos_sorted.txt and 409493_intrachrom_allchrs_KR_100000_octbi_genom_dist_orthos_sorted.txt
```
### Sort ***E. scolopes*** and ***O. bimaculoides*** genomic distance output files and add headers

```bash
(echo -e "Orth_pair_based_on_eupsc_interaction_matrix\tGenomic_distance_eupsc_bp"; tail -n +2 eupsc_genom_dist_orthos.txt | sort) > 409493_intrachrom_allchrs_KR_100000_eupsc_genom_dist_orthos_sorted.txt
(echo -e "Orth_pair_based_on_eupsc_interaction_matrix\tGenomic_distance_octbi_bp"; tail -n +2 octbi_genom_dist_orthos.txt | sort) > 409493_intrachrom_allchrs_KR_100000_octbi_genom_dist_orthos_sorted.txt
```

### Merge sorted genomic distance files

```bash
join -t $'\t' -1 1 -2 1 eupsc_genom_dist_orthos_sorted.txt octbi_genom_dist_orthos_sorted.txt > 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged.txt
```

### Remove duplicate gene pairs

Even though duplicate gene pairs may have different multiple different interaction frequencies, they all have the same genomic distance, so these can be removed from genomic distance files.

```bash
awk '!seen[$0]++' 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged.txt > 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged_rm_dups.txt
```

### Add ***S. officinalis*** genomic distances to the merged file

Again, this could be done seperately and with a simpler script than the one used for *O. bimculoides* above, because the *S. officinalis* orthologous gene names are the same as *E. scolopes*. So, the script [`add_sof_dists.py`](add_sof_dists.py) was used to add *S. officinalis* distances. The input file for this script already has duplicate gene pairs removed from the previous command so this is not necessary to do this step again. 

```bash
python3 add_sof_dists.py sepof.bed 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_genom_dist_sorted_merged_rm_dups.txt
# Example output:
# Number of S. officinalis distances: 410,738
# Number of S. officinalis gene pairs on different chromosomes not written to outfile: 244,02
# Output file written to: 409493_intrachrom_allchrs_KR_100000_eupsc_octbi_sepof_dists_no_dups.txt
```
## Plot scatterplots to define the interaction frequency threshold for interacting gene pairs

The interaction frequency threshold for interacting gene pairs was determined by examining scatterplots generated by R scripts. The script [`get_interaction_threshold_logged_eupsc_sepof.R`](get_interaction_threshold_logged_eupsc_sepof.R) plots *E. scolopes* vs. *S. officinalis* interaction frequency for every gene pair, and dots are coloured by genomic distance. Similar code was used for comparing *E. scolopes* and *O. bimaculoides*.

## Plot boxplots of genomic distances and barplots summarising interaction and ancestral chromosomal origin categories for gene pairs

Next we used the R script [`classify_int_status_plot_boxplots_barplots_by_pec_status.R`](classify_int_status_plot_boxplots_barplots_by_pec_status.R) to:
- Classify different interacting gene pair categories, where gene pairs with an interaction frequency of 10 or more (normalised mapped reads) were classified as interacting
- Average interaction frequencies for duplicate gene pairs within each interaction category. This should unbias results caused by longer genes, but doing it per interaction category keeps cases of genes that may be e.g. on TAD borders (partially interacting and non interacting)
-  Save the resulting dataframe with interaction statuses and averaged interaction frequencies for future analyses (this file is named **409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt**)
- With the new dataframe, plot boxplots of genomic distances between gene pairs across interaction categories, coloured by *P. maximus* chromosomal origin
- Then, test for significant differences between distances for different *P. maximus* chromosome categories per interaction status using a pairwise Wilcoxon test with BH correction
- Finally, plot barplots and stacked percentage barplots for the number of gene pairs each interaction category and *P. maximus* chromosome status group

## Repeat association analyses of interacting gene pairs

This section outlines the steps taken to identify associations between repetitive elements and regions between gene pairs in each species.

### Run RepeatModeler and RepeaMasker [TBC]

Repeats were identified RepeatModeler v.2.0.674 and RepeatMasker v.4.1.875 using default parameters. Repeatmodeler identifies highly repeated regions and constructs consensus sequences and repeatmasker searches these consensus sequences in the genome an identified their positions.

### Get intergenic start and end positions for gene pairs
  
We outputted intergenic start and end positions for gene pairs using the sctipt ['output_intergenic_start_end_and_dist.py'](output_intergenic_start_end_and_dist.py). This script also recalculates and outputs genomic distance and can be run on any file with the format eupsc_gene1;octbi_gene1, eupsc_gene2;octbi_gene2 in the first column. This script outputs one file per species bedfile inputed.

```bash
python3 output_intergenic_start_end_and_dist.py eupsc.bed octbi.bed sepof.bed 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt
# Example output:
# Number of interactions in eupsc: 59,761  
# Output file written to: 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.txt  
# Number of interactions in octbi: 59,761  
# Output file written to: 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_octbi_intergenic_genom_dist.txt  
# Number of in sepof: 59,761  
# Output file written to: 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_sepof_intergenic_genom_dist.txt

```
Where 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt is the file outputted from the script the R script [`classify_int_status_plot_boxplots_barplots_by_pec_status.R`](classify_int_status_plot_boxplots_barplots_by_pec_status.R) used previosuly.


### Convert output files to BED format and sort

Headers and lines with `NA` values are removed, and the gene pair ID is placed in the last column.

```bash
tail -n +2 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.txt | awk '$0 !~ /NA/ {print $1 "\t" $4 "\t" $5 "\t" $2}' | sort | uniq > 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.bed
```
Then, the output file was sorted for downstream bedtools analysis:

```bash
sort -k1,1 -k2,2n 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.bed > 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.bed
```

### Convert RepeatMasker output to GFF and then BED format

Repeat GFF files were created from RepeatMasker output `.out` files using `rmOutToGFF3.pl` which is documented in the script repeats_to_gff.sh. The `rmOutToGFF3.pl` script was obtained from utils folder of RepeatMasker and can be found on the [RepeatMasker GitHub repo](https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/rmOutToGFF3.pl).

Next, the output gff was converted into bed format:

```bash
awk 'BEGIN{OFS="\t"} { if ($3 == "dispersed_repeat") print $1, $4 - 1, $5, $9 }' repeats_eupsc.gff3 > repeats_eupsc.bed
```

Then, the output bed file was sorted for downstream bedtools analyses:

```bash
sort -k1,1 -k2,2n repeats_eupsc.bed > repeats_eupsc_sorted.bed
```

### Run `bedtools intersect` with a processing script

Intersect intergenic gene pair regions with repeats and pipe the output to the Python script [`bedtools_intersect_to_repeat_score_pipe.py `](bedtools_intersect_to_repeat_score_pipe.py) to reduce issues with output file size. This step is documented in the script [`intersect_interactions_repeats.sh`](intersect_interactions_repeats.sh).

### Get repeat content summaries per interaction category

Use [`get_all_norm_repeats_per_category_summed.R`](get_all_norm_repeats_per_category_summed.R) script to summarise repeat content by repeat type and interaction category. This step is documented in the script [`get_all_repeats_sum_norm.sh`](get_all_repeats_sum_norm.sh). This outputs files for each interaction category with columns representing repeat type, total repeat count, total genomic distance between gene pairs, normalised repeat count. Normalised repeat count is calculated as the number of repeats divided by the number of basepairs between gene pairs.

### Get barplots of simple repeats and RND elements

This was done using the R script [`barplot_simple_vs_rnd_repeats.R`](barplot_simple_vs_rnd_repeats.R). Barplots and Wilcoxon significance tests with BH correction were done to compare the proportion of simple repeats RND elements per interaction status. [TBC, need to check wilcox is correct, should be quoted in paper and tile if so].

## Calculate co-expression of gene pairs across different interaction categories [TBC]

[TBC, check why you used paired = false for correlation expression calculation wicox and change if necess. also add all pairwise comparisons sig diff in the legend for fig 3. is this default??]

The R script [`coexpression_analysis_of_interacting_gene_pairs_and_categories_expression_logged.R`](coexpression_analysis_of_interacting_gene_pairs_and_categories_expression_logged.R) was used to:
- Log and TPM normalise *E. scolopes* expression data across tissues
- Calculate  Pearson’s correlation coefficients for co-expression per gene pair across *E. scolopes* tissues for each interaction category and plot it as a density plot
- Calculate significant differences (Wilcoxon test with BH correction), means, and medians in co-expression coefficients for gene pairs across different interaction categories.




