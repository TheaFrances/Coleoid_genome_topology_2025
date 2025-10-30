# Chromatin loop analysis

## Contents

This folder documents the chromatin loop analyses demonstrated using the *E. scolopes* (stage 29) sample 212493 at 50 and 100 kb resolution.

## Differential loop calling

## Running Mustache for differential loop calling

[Mustache](https://github.com/ay-lab/mustache) was run in differential mode using `.hic` files and absolute BED files for each stage. Below are examples of the commands used:

### Example: Stage 25 vs stage 29 (100 kb resolution)
```bash
python3 mustache/diff_mustache.py \
  -f1 212492_intrachrom.allValidPairs.hic \
  -f2 212493_intrachrom.allValidPairs.hic \
  -pt 0.05 -pt2 0.1 -st 0.8 -r 100000 \
  -o eupsc_25vs29_100k
```

### Additional comparisons

These comparisons were also run at various resolutions (50 kb, 100 kb):

- Stage 20 vs Stage 29
- Stage 20 vs Stage 25

### Merge loops across resolutions
Finally, reproducible loops at 50 kb and 100 kb were merged using the [`merge_loops_in_2_resos.py`](merge_loops_in_2_resos.py) script. Loops in the 100 kb file considered duplicates and removed if they fell within a 50 kb window of those in the 50 kb file. This 50 kb window is specified by the --tolerance parameter.


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

**Example output:**

```
Number of loops in first resolution file = 7
Number of loops in second resolution file = 153
Number of loops in merged file = 159
```


## Extract genes from differential loops

To identify genes located within differential chromatin loop anchors, we used the script [`check_gene_in_bin_diff_loops.py`](check_gene_in_bin_diff_loops.py), which compares loop anchor coordinates to gene locations.

Usage:

python3 check_gene_in_bin_diff_loops.py eupsc.bed eupsc_25vs29_50k+100k.diff_loop1

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


## Extract gene lists from annotated loop files

Once loop files are annotated with genes (e.g. .genes files), gene lists were extracted, cleaned, and deduplicated using the following command:

```bash
awk -F '\t' '{print $3 "\n" $5}' eupsc_25vs29_50k+100k.diff_loop1.genes | tr ', ' '\n' | grep -v '^$' | sort | uniq > eupsc_25vs29_50k+100k.genes_list.txt
```
This:
- Extracts gene ID columns (assumed to be columns 3 and 5)
- Splits comma- or space-separated values into individual lines
- Removes empty lines
- Sorts and deduplicates entries

Repeat for each .genes file to generate clean gene lists per condition or comparison.

## Remove duplicate loops from .genes files based on their interactions

To remove duplicate chromatin loops based on overlapping gene interactions, the script [`remove_loop_gene_replicates.py`](remove_loop_gene_replicates.py) was ran on the .genes files.  
This script identifies and removes redundant loops where the same sets of genes appear multiple times across loop anchors.

### Example: *E. scolopes* 25 vs 29 (50 kb + 100 kb)


```bash
    python3 remove_loop_gene_replicates.py eupsc_25vs29_50k+100k.diff_loop1.genes
done
```

**Example output:**

```
Processing loop file: 25_29_50k+100k.diffloop1.genes
Number of conserved loops in input file =  69
Number of loops in output with duplicates removed =  68
Output written to: 25_29_50k+100k.diffloop1.genes_rm_dups
```

Each processed file produces a new output file with the `_rm_dups` suffix, indicating that duplicate loops have been filtered out.




