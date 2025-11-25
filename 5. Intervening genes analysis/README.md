### Intervening gene analysis in loops and interacting gene pairs

This folder contains scripts and output summaries for quantifying intervening gene content and gene coverage within chromatin loops and interacting gene pairs, using *Euprymna scolopes* as the primary example.

---

### 1. Intervening genes in significant loops

The script `count_intervening_genes_loops_size_norm.py` was used to count intervening genes and calculate the percentage of each loop covered by genes. Input files include:

* A BED file with species-specific gene coordinates
* A TSV file with significant loops

**Summary of script behavior (eupsc example):**

* Total significant loops: 329
* Mean intervening genes: 13.56
* Median intervening genes: 4
* Mean gene coverage: 21.65%
* Median gene coverage: 16.58%

**Output:**
A file with intervening gene and gene coverage statistics for each loop:
`eupsc_29cat_50k+100k.tsv.intervening_genes`

---

### 2. Intervening genes in interacting gene pairs

The script `count_intervening_genes_gene_pair_size_norm.py` was used to quantify intervening genes and loop coverage for gene pairs grouped by interaction status: `interacting`, `not_interacting`, or `none`.

**Example input:**

* A BED file with gene coordinates (e.g., `eupsc.bed`)
* A gene pair list with interaction annotations

**Summary output for eupsc:**

* Interacting:

  * Total gene pairs: 2544
  * Median intervening genes: 1
  * Median gene coverage: 21.43%

* Not interacting:

  * Total gene pairs: 57217
  * Median intervening genes: 163
  * Median gene coverage: 32.97%

* None:

  * Total gene pairs: 59761
  * Median intervening genes: 154
  * Median gene coverage: 32.91%

Filtered out gene pairs >1.5Mb apart.

---

### 3. Intervening genes across four interaction categories

The script `count_intervening_genes_gene_pair_size_norm_4cats.py` was used to summarize intervening gene content across four gene pair categories:

* `interacting_all_species`
* `interacting_deca_only`
* `interacting_octbi_only`
* `not_interacting_any_species`

**Summary output for eupsc:**

* `interacting_all_species`: Median = 1 intervening gene; 8.76% coverage
* `interacting_deca_only`: Median = 2 intervening genes; 35.46% coverage
* `interacting_octbi_only`: Median = 118.5 intervening genes; 33.00% coverage
* `not_interacting_any_species`: Median = 165 intervening genes; 32.97% coverage

**Note:** Median values from R were used in final plots, as Python script did not exclude gene pairs with >1 intervening gene.

---

**Scripts used:**

* `count_intervening_genes_loops_size_norm.py`
* `count_intervening_genes_gene_pair_size_norm.py`
* `count_intervening_genes_gene_pair_size_norm_4cats.py`

All scripts accept a BED file and a loop or gene pair file as input and output a summary file with intervening gene and coverage values.

