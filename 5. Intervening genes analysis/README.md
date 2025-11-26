# Quantify intervening gene coverage in loops and interacting gene pairs

This folder contains scripts and output summaries for quantifying intervening gene content and gene coverage within chromatin loops and interacting gene pairs, using *E. scolopes* as an example.

## Contents

- [Intervening genes in loops](#intervening-genes-in-loops)
- [Intervening genes across four interaction categories](#intervening-genes-across-four-interaction-categories)
- [Compare gene coverage in gene pairs and loops](#compare-gene-coverage-in-gene-pairs-and-loops)

### Intervening genes in loops

The script [`count_intervening_genes_loops_size_norm.py`](count_intervening_genes_loops_size_norm.py) was used to count intervening genes and calculate the percentage of each loop covered by genes. Input files include a BED file with species-specific gene coordinates (see example [here](../Test%20input%20files/eupsc.bed)) and TSV Mustache output file generated from the [`3. Chromatin loop analysis`](../3.%20Chromatin%20loop%20analysis) step with significant loops. This outputs a file with with the suffix 'intervening_genes' and **inner_loop_size (bp)**, **intervening_gene_ids**, **intervening_gene_count**, **total_gene_bp**, **gene_coverage_percentage** columns added to the input TSV.

**Example command** (*E. scolopes* sample 409493 example):
```bash
python3 count_intervening_genes_loops_size_norm.py eupsc.bed eupsc_loops_50k+100k.tsv
```

**Example output:**

* Total significant loops: 329
* Mean intervening genes: 13.56
* Median intervening genes: 4
* Mean gene coverage: 21.65%
* Median gene coverage: 16.58%
* Output written to: eupsc_loops_50k+100k.tsv.intervening_genes

### Intervening genes across four interaction categories

The script [`count_intervening_genes_gene_pair_size_norm_4cats.py`](count_intervening_genes_gene_pair_size_norm_4cats.py) was used to calculate the number of intervening genes and the proportion of the gene pair span covered by genes within gene pairs for the four interaction categories per species:  *Interacting across the coleoids*, *Decapodiform-only interacting*, *O. bimaculoides-only interacting*, and *Not in conserved coleoid interacting gene pair*. Input files include a BED file with species-specific gene coordinates (see example [here](../Test%20input%20files/eupsc.bed)) and the merged interaction threshold file containing gene pair information and interaction categories generated in the [`2. Interacting gene pair analysis`](../2.%20Interacting%20gene%20pair%20analysis) step. This outputs a file with the species name and interaction category and 'intervening_genes' as the suffix e.g. eupsc_interacting_all_species_intervening_genes.txt, and with the columns **intervening_genes**, **intervening_count**, **total_gene_bp**, **gene_coverage_percentage** added to the input file.


**Example command** (*E. scolopes*):

```bash
python3 count_intervening_genes_gene_pair_size_norm_4cats.py eupsc.bed 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.txt  interacting_all_species
```

**Example output (*E. scolopes*, interacting_all_species)**

* Total gene pairs: 1592
* Mean intervening genes: 3.75
* Median intervening genes: 1
* Mean gene coverage: 26.60%
* Median gene coverage: 8.76%
* Output written to: 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof.eupsc_interacting_all_species_intervening_genes.txt

> **Note:** Final plots and statistics use the *R-derived* values after filtering gene pairs with more than one gene per bin. The Python-derived medians were not reported in the paper.


### Compare gene coverage in gene pairs and loops

The R script [`boxplots_intervening_genes_gene_pairs_vs_loops.R`](boxplots_intervening_genes_gene_pairs_vs_loops.R) compares gene coverage percentages across different gene pair interaction categories and chromatin loops, using outputs from previous intervening gene analyses.

**Summary of script functionality**:

- Loads gene pair coverage files for four interaction categories per species:  
  *Interacting across the coleoids*, *Decapodiform-only interacting*, *O. bimaculoides-only interacting*, and *Not in conserved coleoid interacting gene pair*.  
- Filters gene pairs to include only those with at least one intervening gene.
- Adds `region_type` labels to distinguish loop data from gene pair categories.
- Combines loop and gene pair data into a single dataset for comparison.
- Generates:
  - A boxplot comparing gene coverage (%) across interaction types and species.
  - A combined plot comparing gene coverage in loops vs. gene pairs.
- Performs Wilcoxon rank-sum tests with BH correction to compare:
  - Interaction categories within species.
  - Loops vs. all gene pair categories within species.
- Outputs boxplots, significance tests, and tables of median gene coverage per category

