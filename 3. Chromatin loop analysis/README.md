# Chromatin loop analyses

This folder documents the chromatin loop analyses demonstrated using the *E. scolopes* (stage 29) sample 212493 at 50 and 100 kb resolution, and the differential chromatin loop analyses samples demonstrated using the *E. scolopes* (stage 25) sample 212492 and 212493 (stage 29) at 50 and 100 kb resolution as examples. The analyses of chromatin loops across species were done with higher coverage samples and an example shown is for the sample 409493 (samples 200409 and 212493 merged).

## Contents
- [Get loops and genes in loop anchors](#get-loops-and-genes-in-loop-anchors)
  - [Running Mustache for loop calling](#running-mustache-for-loop-calling)
  - [Merge loops across resolutions](#merge-loops-across-resolutions)
  - [Extract genes from differential loop and total loop files](#extract-genes-from-differential-loop-and-total-loop-files)
  - [Extract gene lists from annotated loop files](#extract-gene-lists-from-annotated-loop-files)
  - [Remove duplicate loops from .genes files based on their interactions](#remove-duplicate-loops-from-genes-files-based-on-their-interactions)
- [Loop anchor gene expression analyses](#loop-anchor-gene-expression-analyses)
  - [Clustered heatmaps of loop-associated gene expression across developmental stages](#clustered-heatmaps-of-loop-associated-gene-expression-across-developmental-stages)
  - [Gene expression comparison across developmental stage-specific loops](#gene-expression-comparison-across-developmental-stage-specific-loops)
- [Loop anchor GO analysis](#loop-anchor-go-analysis)
- [ATAC-seq signal normalisation, peak filtering and file conversion and formatting](#atac-seq-signal-normalisation-peak-filtering-and-file-conversion-and-formatting)
  - [ATAC-seq signal normalisation](#atac-seq-signal-normalisation)
  - [File formatting and peak filtering](#file-formatting-and-peak-filtering)
  - [Converting and subsetting ATAC-seq BigWig files](#converting-and-subsetting-atac-seq-bigwig-files)
- [Plot triangle loop figures with annotation tracks](#plot-triangle-loop-figures-with-annotation-tracks)
  - [Plot loops using Plotgardner](#plot-loops-using-plotgardner)
- [Plotting differential insulation score](#plotting-differential-insulation-score)
- [Prepare files of loops for cross-species comparisons](#prepare-files-of-loops-for-cross-species-comparisons)
- [Check whether loop genes in *E. scolopes* fall on the same or different chromosomes in *O. bimaculoides*](#check-whether-loop-genes-in-e-scolopes-fall-on-the-same-or-different-chromosomes-in-o-bimaculoides)
- [Conserved chromatin loop analyses across species](#conserved-chromatin-loop-analyses-across-species)
  - [Loop calling with Mustache](#loop-calling-with-mustache)
  - [Merge loop calls across resolutions](#merge-loop-calls-across-resolutions)
  - [Annotate loops with genes](#annotate-loops-with-genes)
  - [Remove duplicate gene-pair annotations](#remove-duplicate-gene-pair-annotations)
  - [Identify conserved loops between species](#identify-conserved-loops-between-species)
  - [Identify loops conserved across all three species](#identify-loops-conserved-across-all-three-species)
  - [Check chromosomal status and intergenic distances of conserved loops](#check-chromosomal-status-and-intergenic-distances-of-conserved-loops)
  - [Expression heatmap of conserved loop anchor genes](#expression-heatmap-of-conserved-loop-anchor-genes)
  - [Correlation of conserved loop size with genome size](#correlation-of-conserved-loop-size-with-genome-size)

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
| 9      | Normalised interaction count         |


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

This analysis visualises expression patterns of genes located within *E. scolopes* chromatin loops for each developmental stage (20, 25, and 29).  
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
   - TPM-normalised expression data were read and averaged across replicates for each developmental stage (Stages 14–28).

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

**Notes:**
- Genes appearing in multiple stages were **excluded** to focus on genes unique to a single stage’s loops.
- Optional sections allow toggling between keeping all genes or filtering for stage-specific ones.

## Loop anchor GO analysis

The script [`GO_analysis_diff_loops.R`](GO_analysis_diff_loops.R) was used to perform gene ontology analyses on loop anchor genes, using files containing gene lists from all loop anchors in each stage, e.g. tot_loops_eupsc_29_50k+100k.tsv.genes_rm_dups_list.txt as the gene list to test for enrichment. **Note** GO analyses were performed the same way as in [GO_analyses_interacting_all_octbi.R](../2.%20Interacting%20gene%20pair%20analysis/GO_analyses_interacting_all_octbi.R)) from the interacting gene pair analysis. 

## ATAC-seq signal normalisation, peak filtering and file conversion and formatting

### ATAC-seq signal normalisation

To allow direct comparison of ATAC-seq signal across developmental stages, we normalised the `.bw` signal tracks to a common mean depth based on `bigWigInfo` statistics. Specifically, we scaled Stage 25 and Stage 29 signals to match the depth of Stage 20, which had the highest mean coverage.

**Mean signal values from `bigWigInfo`:**

| Stage | Mean signal | Relative depth          | Interpretation   |
|-------|-------------|-------------------------|------------------|
| 20    | **3.43**    | highest                 | deepest coverage |
| 25    | **3.35**    | ≈ 0.98×                 | nearly equal     |
| 29    | **2.81**    | ≈ 0.82×                 | shallower        |

**Scaling factors (relative to Stage 20):**

| Stage | Mean  | Scaling factor |
|-------|-------|----------------|
| 20    | 3.428 | **1.00**       |
| 25    | 3.352 | **1.02**       |
| 29    | 2.815 | **1.22**       |

Scaling was performed using [`bigwigCompare`](https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html) from deepTools:
```bash
# Stage 25 → scale to Stage 20
bigwigCompare -b1 97308_Stage25_bbduk_aln.bw -b2 97308_Stage25_bbduk_aln.bw \
  --scaleFactors 1.02288469:1 --operation first \
  -o 97308_Stage25_bbduk_aln_normalised.bw

# Stage 29 → scale to Stage 20
bigwigCompare -b1 97309_Stage29_bbduk_aln.bw -b2 97309_Stage29_bbduk_aln.bw \
  --scaleFactors 1.21804881:1 --operation first \
  -o 97309_Stage29_bbduk_aln_normalised.bw
  ```

### File formatting and peak filtering

Normalised files were then converted to bedGraph:
```bash
bigWigToBedGraph 97309_Stage29_bbduk_aln_normalised.bw 97309_Stage29_bbduk_aln_normalised.bedgraph
bigWigToBedGraph 97308_Stage25_bbduk_aln_normalised.bw 97308_Stage25_bbduk_aln_normalised.bedgraph
bigWigToBedGraph 97307_Stage20_bbduk_aln.bw 97307_Stage20_bbduk_aln_normalised.bedgraph  # For consistency
```

Next, regions were filtered for signal > 1, i.e. regions of open chromatin:
```bash
awk '$4 > 1' 97309_Stage29_bbduk_aln_normalised.bedgraph > 97309_Stage29_bbduk_aln_normalised_peaks.bedgraph
```
The same command was then applied to the other stages as well (same with subsequent commands).

Next, files were converted to 3-column BED format:
```bash
cut -f1-3 97309_Stage29_bbduk_aln_normalised_peaks.bedgraph > 97309_Stage29_bbduk_aln_normalised_peaks_3col.bed
```

Then, unplaced scaffolds were removed:
```bash
grep Lachesis 97309_Stage29_bbduk_aln_normalised_peaks_3col.bed > 97309_Stage29_bbduk_aln_normalised_peaks_3col_rm_scaf.bed
```

BED files were sorted using a sorted chromosome size file (with unplaced scaffolds removed):
```bash
bedtools sort -i 97309_Stage29_bbduk_aln_normalised_peaks_3col_rm_scaf.bed \
  -g ../Lachesis_assembly_chrsize_rm_scaf.sorted.txt \
  > 97309_Stage29_bbduk_aln_normalised_peaks_3col_sorted_rm_scaf.bed
  ```

### Converting and subsetting ATAC-seq BigWig files

To handle large `.bw` files for visualisation or downstream analysis, we subset individual chromosomes and from whole-genome normalised bigWig files to reduce file size for when uploading to Plotgardener. The script used was [`subsetBigWig.py`](https://gist.github.com/dpryan79/740f2d00ce6b509ab9644fc43418996c), which requires the `pyBigWig` and `numpy` packages.

Example using the chromosome Lachesis_group4__56_contigs__length_174030884 (chromosome 5 in *E. scolopes*):
```bash
python subsetBigWig.py /path/to/97309_Stage29_bbduk_aln_normalised.bw /output/97309_st29_chr5.bw Lachesis_group4__56_contigs__length_174030884
```

## Plot triangle loop figures with annotation tracks

### Plot loops using Plotgardner

To generate triangle heatmaps of chromatin loops alongside gene models and regulatory tracks (e.g., ATAC-seq), we used the R package [`plotgardner`](https://phanstiellab.github.io/plotgardner/). The R script [`plot_diff_loop_plotgardener_with_atac.R`](plot_diff_loop_plotgardener_with_atac.R)  was used to create Figure 4D of the manuscript, highlighting a developmentally dynamic loop in *E. scolopes* chromosome 5.

The script performs the following steps:

1. Loads a `.hic` file containing KR-normalised Hi-C data for the relevant developmental stage.
2. Builds a genome assembly using:
   - A `.fasta` genome file  
   - A `TxDb` object created from an exon-corrected GFF  
   - A `BSgenome` object generated from a `.2bit` file  
   - A custom `OrgDb` package (e.g. `org.Eesc.eg.db`)
3. Plots a triangular Hi-C map over a specified genomic region.
4. Overlays gene models using `plotGenes()`, with specific genes of interest highlighted.
5. Loads and plots ATAC-seq signal from normalised `.bw` files for three developmental stages, scaling the signal tracks using the 99.5th percentile to match the y-axis range.
6. Saves the figure.

**Required input files:**

- `.hic` file (KR-normalised Hi-C matrix)  
- GFF file with corrected exon annotations  
- `.2bit` genome sequence file  
- BigWig (`.bw`) files containing normalised ATAC-seq signal  
- Optional: `OrgDb` object for gene metadata (used to display gene names)

This method enables high-quality visualisation of long-range chromatin interactions and their associated regulatory context across developmental stages in *E. scolopes*.

## Plotting differential insulation score

The R script [`genova_diff_insulation_score.R`](genova_diff_insulation_score.R) uses the [GENOVA R package](https://github.com/robinweide/GENOVA) to calculate insulation scores from *E. scolopes* Micro-C iced matrices and visualise differences in insulation scores across the three developmental stages for the loop region shown in Figure 4D. The script then annotates potential differential loop regions to produce Figure 4E. Input files for this script are taken directly from the HiC-Pro mapping output (detailed in [1. Preprocessing and mapping of Micro-C reads](../1.%20Preprocessing%20and%20mapping%20of%20Micro-C%20reads/)).

## Prepare files of loops for cross-species comparisons

### Loop calling with Mustache

Chromatin loops were identified at 50 kb and 100 resolution using the [Mustache](https://github.com/ay-lab/mustache) peak caller. Loops were called on whole-tissue `.hic` files normalised with KR as follows:

```bash
python3 mustache.py -f 409493_intrachrom.allValidPairs.hic -r 50kb -norm KR -pt 0.01 -o eupsc_loops_50k.tsv
```

Where eupsc_50k_loops.tsv is the specified outfile. This was repeated for each species (samples 212489 - *O. bimaculoides* and 992270 - *S. officinalis*) at both 50 and 100 kb resolution.

### Merge loop calls across resolutions

As in the previous loop calling steps, output from resolutions 50 kb and 100 kb resolutions were merged using the script `merge_loops_in_2_resos.py`, allowing a 50 kb tolerance between loop anchors:

```bash
python3 merge_loops_in_2_resos.py eupsc_loops_50k.tsv eupsc_loops_100k.tsv eupsc_loops_50k+100k.tsv --tolerance 50000
```

### Annotate loops with genes

Also as in the previous loop calling steps, genes located in the start and end bins of each loop were identified using `check_gene_in_bin_diff_loops.py`. The input is a BED file of gene coordinates and the loop file:

```bash
python3 check_gene_in_bin_loops.py eupsc.bed eupsc_loops_50k+100k.tsv
```

This reports how many loop bins contain genes at the start and end, and outputs a `.genes` file.

### Remove duplicate gene-pair annotations

Also as in the previous loop calling steps, to avoid double-counting of orthologous gene pairs that appear in multiple loop calls, duplicate gene-pair entries were removed using:

```bash
python3 remove_loop_gene_replicates.py eupsc_loops_50k+100k.tsv.genes
```

This generates a `.genes_rm_dups` file for downstream analysis.

## Check whether loop genes in *E. scolopes* fall on the same or different chromosomes in *O. bimaculoides*

To assess whether interacting gene pairs in *E. scolopes* loops are located on the same chromosome in *O. bimaculoides*, the script [`loop_chrom_status.py`](loop_chrom_status.py) was used. This script takes in a filtered loop-gene file, a 2-column ortholog mapping file (with *E. scolopes* gene IDs in the first column and *O. bimaculoides* gene IDs in the second), and a BED file of gene coordinates in *O. bimaculoides*.

**Example command:**

```bash
python3 loop_chrom_status.py \
  EUPgeneOBI.txt \
  octbi.bed eupsc_50k+100k.tsv.genes_rm_dups
  ```

  #YOU ARE HERE 

## Conserved chromatin loop analyses across species

To identify chromatin loops that are conserved across coleoid cephalopods (*Euprymna scolopes*, *Octopus bimaculoides*, and *Sepia officinalis*), the following multi-step pipeline was used:

### Identify conserved loops between species

To identify conserved loops between species, ortholog files were used (e.g. `EUPgeneSOF.txt`, `EUPgeneOBI.txt` with *E. scolopes* genes in the first column and *S officinalis* or *O. bimaculoides* genes in the second column. The script [`check_conserved_loops.py`](check_conserved_loops.py) defines a loop as conserved if there is at least one orthologous gene pair in bin 1 and one in bin 2 (allowing reversed bin orientation across species):

```bash
python3 check_conserved_loops.py EUPgeneSOF.txt eupsc_loops_50k+100k.genes_rm_dups  sepof_loops_50k+100k.genes_rm_dups eupsc_sepof_consv_loops_50k+100k.tsv
```

Where eupsc_sepof_consv_loops_50k+100k.tsv is the specified output file.

> **Note:** Even after duplicate removal, some loops may share some of the same the same orthologous gene pairs. These should be manually reviewed and considered as a single conserved loop where appropriate.

### Identify loops conserved across all three species

To identify loops conserved across *all three* species, orthologous gene IDs found in conserved loops between species pairs were cross-referenced. For a loop to be considered fully conserved, orthologous genes must appear in both bins across all three species.

### Check chromosomal status and intergenic distances of conserved loops

To assess whether genes involved in conserved loop anchors remain on the same chromosome in other species, the script [`check_chrom_of_conserved_loops_and_dist.py`](check_chrom_of_conserved_loops_and_dist.py) was used. This script checks whether the orthologous genes in each loop anchor are found on the same or different chromosomes in the compared species.

- If **all genes** in a loop are located on the **same chromosome**, the script calculates the **shortest genomic distance** between the genes at the start and end of the loop.
- If the genes are found on **different chromosomes**, the loop is considered to have been disrupted by a chromosomal rearrangement.

**Example command:**
```bash
python3 check_chrom_of_conserved_loops_and_dist.py EUPgeneOBI.txt octbi.bed eupsc_sepof_consv_loops_50k+100k.tsv
```
In this example:
- `eupsc_sepof_consv_loops_50k+100k.tsv` contains loops conserved between *E. scolopes* and *S. officinalis*.
- The script checks the chromosomal locations of the orthologous genes in *O. bimaculoides*, using `EUPgeneOBI.txt` (ortholog mapping) and `octbi.bed` (gene coordinates).

**Example output:**
```bash
Checking chromosome status of conserved loops in: octbi
Number of conserved loops in input file =  28
Number of conserved loops with genes on the same chromosome =  20
Number of conserved loops with genes on different chromosomes =  4
Output written to: eupsc_sepof_consv_loops_50k+100k.tsv.octbi_chrom_status_dist.tsv
```

The output file provides a summary of:
- How many conserved loops are intact (genes still on the same chromosome).
- How many loops have likely been disrupted (genes on different chromosomes).
- The minimum genomic distances between start and end loop anchor genes for loops where all genes remain on the same chromosome.


### Expression heatmap of conserved loop anchor genes

The R script [`conserved_loop_exp_heatmap.R`](conserved_loop_exp_heatmap.R) visualises the expression of orthologous genes found in chromatin loops conserved across species across a range of *E. scolopes* tissues. The script shows an example for loops conserved across the decapodiform species and performs the following steps:

- Loads TPM-normalised gene expression data across multiple tissues.  
- Log-transforms the TPM values (adding a pseudocount of 0.01).  
- Averages expression values across biological replicates for each tissue.  
- Extracts the *E. scolopes* gene IDs (beginning with `cluster_`) found in conserved loops using the file `eupsc_sepof_consv_loops_50k+100k.txt`.  
- Filters the expression matrix to include only genes found in conserved loops.  
- Generates a heatmap of log-transformed, row-scaled expression values using the [`pheatmap`](https://cran.r-project.org/web/packages/pheatmap/index.html) package.  
- Saves the plot.  

**Note:** Tissue names in the expression matrix are manually mapped to more readable names, and replicates are averaged per tissue group before plotting.

This visualisation is used to explore whether conserved loop anchor genes show tissue-specific expression patterns, particularly in brain and neural tissues.

##Get conserved loop sizes

### Calculate loop sizes for conserved loops

The R script [`conserved_loop_sizes.R`](conserved_loop_sizes.R) calculates loop sizes for chromatin loops conserved between species by subtracting the start coordinate of the first loop bin from the end coordinate of the second bin. The script integrates loop‑size data across *E. scolopes*, *S. officinalis*, and *O. bimaculoides*, matching them to orthologous loop pairs identified in the conserved loop analyses.

**Summary of script functionality:**

1. **Input data processing**
   - Reads per‑species loop‑size tables from previously processed `.genes_rm_dups` files for *E. scolopes*, *S. officinalis*, and *O. bimaculoides*.
   - Standardises column names and ensures numeric coordinate formats for the start position of the first loop bin.
   - Uses a helper function to parse loop IDs of the form `<chromosome>:<start>-<end>` into separate chromosome and start‑position columns.
   - These parsed identifiers enable consistent joining of loop‑size information across species.

2. **Getting conserved loop sizes**
   - For each species pair (*E. scolopes–S. officinalis*, *E. scolopes–O. bimaculoides*, and *S. officinalis–O. bimaculoides*), the script joins loop‑size data by matching chromosome and bin‑start coordinates.
   - This produces new files such as:
     - `eupsc_sepof_consv_loops_with_sizes.tsv`
     - `eupsc_octbi_consv_loops_with_sizes.tsv`
     - `sepof_octbi_consv_loops_with_sizes.tsv`
   - A combined file containing loops conserved across all three species is read in and joined with loop‑size tables for each species.
   - The resulting table (`eupsc_octbi_sepof_consv_loops_with_sizes.tsv`) includes loop‑size information from all three species side by side.

3. **Diagnostics**
   - A simple summary table reports:
     - The number of loops processed per comparison.
     - How many loops are missing size information in each species (useful for detecting coordinate mismatches).

**Example output table (diagnostic summary):**
```
# A tibble: 4 × 5
  dataset       n_rows n_missing_eupsc n_missing_sepof n_missing_octbi
  <chr>           <int>           <int>           <int>           <int>
1 eupsc-sepof        65               0               3              NA
2 eupsc-octbi        72               1              NA               2
3 sepof-octbi        68              NA               0               1
4 all-three          45               0               2               1
```
**Note:**
- Minor coordinate mismatches between files can lead to missing loop‑size values, which were manually corrected when due to small shifts in overlapping bins.

### Correlation of conserved loop size with genome size

The R script [`conserved_loop_size_correlation.R`](conserved_loop_size_correlation.R) assesses whether the physical size of conserved chromatin loops scales with overall genome size across species.

**Input:**  
The input file is `loop_sizes_start_end.txt`, which is a simplified version of the full three‑species conserved loop table (`eupsc_octbi_sepof_consv_loops_with_sizes.tsv`). It contains the following columns, with one row per conserved loop, where missing values indicate the loop is not conserved in that species:

```
eupsc_size    octbi_size    sepof_size
130000        600000        1450000
6900000       2300000        5950000
...
```

**Summary of analysis steps:**

1. **Input data processing**
   - Adds a loop ID to each row to track conserved loops.
   - Removes known outlier loops (e.g., a single loop in *E. scolopes* which is likely due to a genome missassembly).
   - Pivots the wide format into long format (one row per species per loop).
   - Assigns genome size values to each species.

2. **Per-loop correlation analysis**
  - **Linear regression slope**, estimating bp increase per 1 Gb genome size.
   - Summarises the number and proportion of loops showing a positive correlation with genome size.
   - Performs Wilcoxon signed-rank tests on the distribution of slopes (vs. 0).
   - This tests whether conserved loops tend to scale with genome size across species.

3. **Scatter plot of individual loops**
   - Generates a clean scatter/line plot showing each loop's size across species. Species genome sizes are shown on the x-axis; loop sizes on the y-axis.
   - Each loop is assigned a unique colour using a distinct HCL-based palette.
   - The final figure is saved.

