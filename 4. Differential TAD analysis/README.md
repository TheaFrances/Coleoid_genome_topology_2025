# Differential TAD Analysis

This README documents the pipeline used for differential TAD analysis, using the data for *E. scolopes* developmental stages as an example.

## Prepare files

### Step 1: Extract KR-normalised intrachromosomal matrices from .hic files

To perform differential TAD calling with TADCompare, we first extracted KR-normalized contact matrices for each chromosome from `.hic` files at 100 kb resolution using `juicer_tools dump observed KR`. This was done separately for each developmental stage: stage 20 (sample 320995), stage 25 (sample 212492), stage 29 (sample 212493). Each script loops over all chromosomes, generates a `.dumped.hic` file per chromosome, and writes output to a dedicated folder.

### Dump individual chromosomes

Example command for *E. scolopes* sample 29 (212492):

```bash
 #!/bin/bash

# Define the command and its components
juicer_tools="juicer_tools_1.22.01.jar"
input_file="212492_intrachrom.allValidPairs.hic"
output_directory="/dumped_100k/"
bin_size="100000"

# Define an array of chromosome names
chromosomes=(
    "Lachesis_group0__68_contigs__length_203855924"
    "Lachesis_group1__63_contigs__length_187967259"
    "Lachesis_group10__55_contigs__length_138218214"
    "Lachesis_group11__73_contigs__length_131769366"
    "Lachesis_group12__46_contigs__length_133395132"
    "Lachesis_group13__58_contigs__length_132030582"
    "Lachesis_group14__44_contigs__length_129294686"
    "Lachesis_group15__55_contigs__length_125787758"
    "Lachesis_group16__53_contigs__length_120303996"
    "Lachesis_group17__52_contigs__length_118315660"
    "Lachesis_group18__53_contigs__length_106174747"
    "Lachesis_group19__43_contigs__length_108131866"
    "Lachesis_group2__55_contigs__length_182973478"
    "Lachesis_group20__56_contigs__length_105542877"
    "Lachesis_group21__43_contigs__length_103963499"
    "Lachesis_group22__46_contigs__length_104782409"
    "Lachesis_group23__35_contigs__length_97660015"
    "Lachesis_group24__36_contigs__length_98752639"
    "Lachesis_group25__42_contigs__length_94613958"
    "Lachesis_group26__46_contigs__length_93059273"
    "Lachesis_group27__44_contigs__length_93323795"
    "Lachesis_group28__49_contigs__length_91470296"
    "Lachesis_group29__46_contigs__length_89259793"
    "Lachesis_group3__52_contigs__length_175295898"
    "Lachesis_group30__51_contigs__length_87052432"
    "Lachesis_group31__46_contigs__length_81539927"
    "Lachesis_group32__40_contigs__length_78634074"
    "Lachesis_group33__48_contigs__length_77384337"
    "Lachesis_group34__47_contigs__length_76047342"
    "Lachesis_group35__42_contigs__length_70637049"
    "Lachesis_group36__35_contigs__length_73496112"
    "Lachesis_group37__50_contigs__length_69233352"
    "Lachesis_group38__38_contigs__length_61810156"
    "Lachesis_group39__38_contigs__length_50967332"
    "Lachesis_group4__56_contigs__length_174030884"
    "Lachesis_group40__41_contigs__length_53238745"
    "Lachesis_group41__39_contigs__length_54644242"
    "Lachesis_group42__34_contigs__length_44581599"
    "Lachesis_group43__26_contigs__length_23511657"
    "Lachesis_group44__24_contigs__length_15433969"
    "Lachesis_group45__11_contigs__length_9560946"
    "Lachesis_group46__3_contigs__length_1604139"
    "Lachesis_group47__3_contigs__length_1393269"
    "Lachesis_group5__70_contigs__length_158228327"
    "Lachesis_group6__62_contigs__length_156955521"
    "Lachesis_group7__52_contigs__length_145002158"
    "Lachesis_group8__53_contigs__length_148332072"
    "Lachesis_group9__67_contigs__length_144638433"
)

# Loop through each chromosome and execute the command
for chromosome in "${chromosomes[@]}"
do
    output_file="${output_directory}212492_${chromosome}_KR_${bin_size}.dumped.hic"
    java -jar "$juicer_tools" dump observed KR "$input_file" "$chromosome" "$chromosome" BP "$bin_size" > "$output_file"
    echo "Processed chromosome $chromosome"
done

```

### Clean dumped matrix files

The dumped matrix files sometimes contain Juicer warning lines (e.g. `WARN: ...`). These lines were removed using:

```bash
for file in *.dumped.hic; do
    grep -v "^WARN" "$file" > "${file%.dumped.hic}_cleaned.dumped.hic"
done
```

Uncleaned files were deleted to avoid redundancy:

```bash
rm *100000.dumped.hic
```
## Differential TAD Analysis in TADCompare

The script [`TAD_compare.R`](TAD_compare.R) was used for differential TAD analysis using the [`TADcompare R package`](https://github.com/dozmorovlab/TADCompare).

**Summary of script functionality**:

- loads the KR-normalised contact matrix dumps for all chromosomes and tissues (100 kb resolution).
- Performs pairwise comparisons between all tissues using the TADCompare R package, with a 1 Mb sliding window.
- Converts sparse matrices to dense format and analyses TAD structure differences for each chromosome and comparison.
- Categorises change types (e.g. Split, Merge, Strength Change) and normalises them as percentages per chromosome and across the genome.
- Outputs a summary table with differential TAD counts and percentages and a stacked bar plot showing differential TAD distributions across comparisons.

## Upset plot of shared differential TADs across samples

The script [`diff_TAD_overlap_upset_plot.R`](diff_TAD_overlap_upset_plot.R) visualises the overlap of differential TAD boundaries across all stage comparisons in *E. scolopes* using an UpSet plot.

**Summary of script functionality**:

- Loads all differential TAD tables from differential_tad files generated in the previous [`TAD_compare.R`](TAD_compare.R) R script.
- Keeps only TADs that are marked as differential in at least one comparison and removes any with NA.
- Assigns unique TAD IDs from chromosome and boundary coordinates. 
- Converts data to a binary matrix (1 = present, 0 = absent) per comparison.
- Generates an UpSet plot showing intersections of differential TADs across comparisons and outputs it. 


