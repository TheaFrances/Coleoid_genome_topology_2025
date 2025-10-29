# Interacting gene pair analysis

## Contents

This folder documents the chromatin loop analyses demonstrated using the *E. scolopes* (stage 29) sample 212493 at 50 and 100kb resolution.

## Differential loop calling

## Running Mustache for Differential Loop Calling

Mustache was run in differential mode using `.hic` files and absolute BED files for each stage. Below are examples of the commands used:

### Example: Stage 25 vs Stage 29 (100 kb resolution)
```bash
python3 mustache/diff_mustache.py \
  -f1 212492_intrachrom.allValidPairs.hic \
  -f2 212493_intrachrom.allValidPairs.hic \
  -pt 0.05 -pt2 0.1 -st 0.8 -r 100000 \
  -o eupsc_25vs29
```

### Additional Comparisons

These comparisons were also run at various resolutions (50 kb, 100 kb):

- Stage 20 vs Stage 29
- Stage 20 vs Stage 25



