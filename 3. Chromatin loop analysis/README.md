# Chromatin loop analyses

This folder documents the chromatin loop analyses demonstrated using the *E. scolopes* (stage 29) sample 212493 at 50 and 100 kb resolution, and the differential chromatin loop analyses samples demonstrated using the *E. scolopes* (stage 25) sample 212492 and 212493 (stage 29) at 50 and 100 kb resolution.


## Contents
- [Get loops and genes in loop anchors](#Get-loops-and-genes-in-loop-anchors)
- [Loop anchor gene expression analyses](#Loop-anchor-gene-expression-analyses)

## Get loops and genes in loop anchors

### Running Mustache for loop calling

[Mustache](https://github.com/ay-lab/mustache) was run in differential mode using `.hic` files, to produce files of differential loops across developmental stages, as well as the total number of loops in each stage.  Loop calling was restricted to annotated chromosomes, and datasets containing unplaced scaffolds were excluded from analysis. Below are examples of the commands used:

**Example: Stage 25 vs stage 29 (100 kb resolution)**
```bash
python3 mustache/diff_mustache.py \
  -f1 212492_intrachrom.allValidPairs.hic \
  -f2 212493_intrachrom.allValidPairs.hic \
  -pt 0.05 -pt2 0.1 -st 0.8 -r 100000 \
  -o eupsc_25vs29_100k
```
**Additional comparisons**

These comparisons were also run at various resolutions (50 kb, 100 kb):

- Stage 20 vs Stage 29
- Stage 20 vs Stage 25

Mustache outputs chromatin loops in a simple **tab-separated format** with the following columns:

| Column | Description                          |
|--------|--------------------------------------|
| 1      | Chromosome of start anchor           |
| 2      | Start position of start anchor (bin) |
| 3      | End position of start anchor (bin)   |
| 4      | Chromosome of end anchor             |
| 5      | Start position of end anchor (bin)   |
| 6      | End position of end anchor (bin)     |
| 7      | q-value of the loop                  |
| 8      | log10(p-value)                       |
| 9      | Normalized interaction count         |


When comparing loop calls across stages using Mustache’s differential mode, the following output files are generated:

| File name           | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `output.loop1`      | Loops found in the first dataset (e.g. `data1.hic`)                         |
| `output.loop2`      | Loops found in the second dataset (e.g. `data2.hic`)                        |
| `output.diffloop1`  | Loops **present in `data1.hic`** but **weakened or absent** in `data2.hic` |
| `output.diffloop2`  | Loops **present in `data2.hic`** but **weakened or absent** in `data1.hic` |


### Merge loops across resolutions
Finally, reproducible loops at 50 kb and 100 kb were merged using the [`merge_loops_in_2_resos.py`](merge_loops_in_2_resos.py) script. Loops in the 100 kb file considered duplicates and removed if they fell within a 50 kb window of those in the 50 kb file. This 50 kb window is specified by the --tolerance parameter.

**Example: Stage 25 vs stage 29 50 kb + 100 kb resolution**
```bash
python merge_loops_in_2_resos.py \
  eupsc_25vs29_50k.diff_loop1 \
  eupsc_25vs29_100k.diff_loop1 \
  eupsc_25vs29_50k+100k.diff_loop1 \
  --tolerance 50000

  python merge_loops_in_2_resos.py \
  eupsc_25vs29_50k.diff_loop2 \
  eupsc_25vs29_100k.diff_loop2 \
  eupsc_25vs29_50k+100k.diff_loop2 \
  --tolerance 50000
```

This command generates an output file with the name specified in the command, e.g., eupsc_25vs29_50k+100k.diff_loop1. 

**Note** the same procedure was applied to files with the suffixes loop1 and loop2, which contain all detected loops for each sample, not just the differential ones, for example:

```bash
python merge_loops_in_2_resos.py \
eupsc_25vs29.loop2 \
eupsc_25vs29.loop2 \
tot_loops_eupsc_29_50k+100k.tsv \
  --tolerance 50000
  ```

**Example output:**

```
Number of loops in first resolution file = 7
Number of loops in second resolution file = 153
Number of loops in merged file = 159
```

### Extract genes from differential loop and total loop files

To identify genes located within differential chromatin loop anchors, the script [`check_gene_in_bin_loops.py`](check_gene_in_bin_loops.py) was used on both differential loop files and  the files of total loops in each stage/ This script compares loop anchor coordinates to gene locations.

**Examples:**

```bash
python3 check_gene_in_bin_loops.py eupsc.bed eupsc_25vs29_50k+100k.diff_loop1
python3 check_gene_in_bin_loops.py eupsc.bed tot_loops_eupsc_29_50k+100k.tsv
```

- Where the first input file is a BED file with gene coordinates (e.g. eupsc.bed for *E. scolopes*)
- And the second input file is a differential loop file in Mustache output format with anchor coordinates

This script prints:
- Number of differential loops
- Number of loop anchors containing genes in the start bin
- Number of loop anchors containing genes in the end bin
- Number of loops with genes in both bins

**Example output:**

```
Number of significant, differential loops =  37
Number of loop bins with genes in (loop start) =  22
Number of loop bins with genes in (loop end) =  23
Number of loop bins with genes in both bins =  15
```
Each processed file produces a new output file with the `.genes` suffix, indicating that genes in loop anchors have been added to the file.


### Extract gene lists from annotated loop files

Once loop files are annotated with genes (e.g. .genes files), gene lists were extracted, cleaned, and deduplicated using the following command:

**Examples:**
```bash
awk -F '\t' '{print $3 "\n" $5}' eupsc_25vs29_50k+100k.diff_loop1.genes | tr ', ' '\n' | grep -v '^$' | sort | uniq > eupsc_25vs29_50k+100k.genes_list.txt
awk -F '\t' '{print $3 "\n" $5}' tot_loops_eupsc_29_50k+100k.tsv.genes | tr ', ' '\n' | grep -v '^$' | sort | uniq > tot_loops_eupsc_29_50k+100k.tsv.genes_list.txt
```
This:
- Extracts gene ID columns (assumed to be columns 3 and 5)
- Splits comma- or space-separated values into individual lines
- Removes empty lines
- Sorts and deduplicates entries

Repeat for each .genes file to generate clean gene lists per condition or comparison.

### Remove duplicate loops from .genes files based on their interactions

To remove duplicate chromatin loops based on overlapping gene interactions, the script [`remove_loop_gene_replicates.py`](remove_loop_gene_replicates.py) was ran on the .genes files.  
This script identifies and removes redundant loops where the same sets of genes appear multiple times across loop anchors.

**Examples:**

```bash
    python3 remove_loop_gene_replicates.py eupsc_25vs29_50k+100k.diff_loop1.genes
    python3 remove_loop_gene_replicates.py tot_loops_eupsc_29_50k+100k.tsv.genes
```
**Example output:**
```
Processing loop file: 25_29_50k+100k.diffloop1.genes
Number of conserved loops in input file =  69
Number of loops in output with duplicates removed =  68
Output written to: 25_29_50k+100k.diffloop1.genes_rm_dups
```

Each processed file produces a new output file with the `_rm_dups` suffix, indicating that duplicate loops have been filtered out.

## Loop anchor gene expression analyses

### Clustered heatmaps of loop-associated gene expression across developmental stages

This analysis visualizes expression patterns of genes located within *E. scolopes* chromatin loops for each developmental stage (20, 25, and 29).  
The script [`heatmaps_tau_coexp.R`](heatmaps_tau_coexp.R) was run **separately for each stage**, using the corresponding loop–gene mapping file.

**Summary of script functionality:**

This script:

- **Loads and processes loop-gene mappings**  
  Loops from all three stages (20, 25, 29) are read and mapped to their associated genes using unique loop IDs derived from genomic bin coordinates.

- **Identifies stage-specific loops**  
  For the focal stage (e.g. stage 29), only loops with unique gene content not shared with other stages are retained for visualisation.

- **Merges expression data**  
  TPM-normalised expression data are log-transformed and merged with the loop–gene mappings to generate per-loop gene expression matrices.

- **Generates a clustered heatmap**  
  A heatmap is created using `pheatmap`, displaying expression of genes grouped by loop.  
  Genes appearing in multiple loops are shown separately for each loop occurrence, with `gaps_row` separating different loops.

- **Calculates Tau tissue-specificity**  
  The tissue-specificity index (Tau) is computed for each gene based on log-transformed expression values, then summarised per loop (mean and median Tau).

- **Calculates intra-loop co-expression**  
  For each loop, mean pairwise Pearson correlation between genes is calculated to assess co-regulation within loop domains.

- **Outputs summary statistics and visualizations**  
  Final outputs include:
  - A heatmap for each stage showing loop-level gene expression patterns  
  - Summary statistics for Tau and co-expression across loops  

Each stage (20, 25, and 29) was analysed independently, producing separate heatmaps and summary outputs.

### Gene expression comparison across developmental stage-specific loops

The script [`loop_exp_boxplots.R`](loop_exp_boxplots.R) compares gene expression profiles of loop-associated genes that are **unique to stage-specific loops** (Stages 20, 25, and 29) in *E. scolopes*. Expression values are averaged per developmental stage, log-transformed, and statistically compared using Wilcoxon rank-sum tests.

**Overview of analysis steps**:

1. **Input preparation:**
   - Gene lists were extracted from stage-specific loop–gene mappings (`*.genes_list.txt`).
   - TPM-normalized expression data were read and averaged across replicates for each developmental stage (Stages 14–28).

2. **Stage-specific filtering:**
   - For each stage (20, 25, 29), expression data were subset to include only genes uniquely present in that stage’s loop set (i.e., genes not found in loops from other stages).

3. **Expression matrix preparation:**
   - Data were reshaped into long format (gene × tissue × stage).
   - Expression values were log-transformed using a pseudocount of 0.01.

4. **Statistical comparisons:**
   - All pairwise comparisons between stages were tested using **Wilcoxon rank-sum tests** within each tissue.
   - Significance levels were annotated using the `ggsignif` package.

5. **Visualization:**
   - A faceted boxplot was created showing gene expression (logTPM) across tissues, grouped by developmental stage.
   - Custom color schemes and angled axis labels were used for clarity.

6. **Output:**
   - A multi-panel boxplot (per tissue) of gene expression for genes from stage-specific loops.
   - A significance summary table including:
     - Wilcoxon test statistic and p-value
     - Mean and median expression values per stage

### Notes:
- Genes appearing in multiple stages were **excluded** to focus on genes unique to a single stage’s loops.
- Optional sections allow toggling between keeping all genes or filtering for stage-specific ones.








