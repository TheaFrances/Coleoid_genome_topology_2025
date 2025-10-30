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

## Extract genes from differential loops

To identify genes located within differential chromatin loop anchors, we used the script [`check_gene_in_bin_diff_loops.py`](check_gene_in_bin_diff_loops.py), which compares loop anchor coordinates to gene locations.

Usage:

python3 check_gene_in_bin_diff_loops.py eupsc.bed eupsc_25vs29_50k+100k.diff_loop1

- Where the first input file is a BED file with gene coordinates (e.g. eupsc.bed for *E. scolopes*)
- And the seciond input file is a differential loop file in Mustache output format with anchor coordinates

This script prints:
- Number of differential loops
- Number of loop anchors containing genes in the start bin
- Number of loop anchors containing genes in the end bin
- Number of loops with genes in both bins

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




