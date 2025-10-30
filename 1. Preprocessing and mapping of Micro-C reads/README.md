# Preprocessing and mapping of Micro-C reads

This folder contains the commands used to process raw Micro-C reads. These steps include trimming, concatenating libraries for higher coverage, and mapping with HiC-Pro, and generating a Juicebox `.hic` file. *E. scolopes* stage 29 samples are used as examples throughout this folder. Where applicable, *S. officinalis* is used instead, for example, in the downsampling step, which was not necessary for *E. scolopes*.


## Contents

- [Downsampling](#downsampling)
- [Trimming](#trimming)
- [Read concatenation](#read-concatenation)
- [Create Bowtie2 index of reference genome for HiC-Pro mapping](#create-bowtie2-index-of-reference-genome-for-hic-pro-mapping)
- [Make file of chromosome sizes for HiC-Pro mapping](#make-file-of-chromosome-sizes-for-hic-pro-mapping)
- [Mapping with HiC-Pro](#mapping-with-hic-pro)
- [Remove interchromosomal interactions from .allValidPairs file](#remove-interchromosomal-interactions-from-allvalidpairs-file)
- [Remove unplaced scaffolds from .allValidPairs file](#remove-unplaced-scaffolds-from-allvalidpairs-file)
- [Generate the .hic file and apply KR normalisation](#generate-the-hic-file-and-apply-kr-normalisation)
- [Dump .hic matrix at 100 kb resolution with KR normalisation](#dump-hic-matrix-at-100-kb-resolution-with-kr-normalisation)
- [Remove lines from dumped matrices and concatenate into one file](#remove-lines-from-dumped-matrices-and-concatenate-into-one-file)

### Downsampling
We downsampled relevant tissue-specific Micro-C samples to approximately 350 million reads to normalise for read depth across tissues.  
For example, for the *S. officinalis* eye tissue, this was applied to read 1 using the command below, where `0.8559` was the downsampling fraction used.

```bash
module load seqtk
seqtk sample -s100 327271_R1.fastq.gz 0.8559 | gzip > /lisc/scratch/molevo/thea/sepof_microc/microc3/trim50_data/ds_EY/327271/327271_R1_ds.fastq.gz
```

### Trimming

Reads were trimmed to 50 bp using **Trimmomatic** with the `CROP:50` setting.  
This step is documented for one stage 29 *E. scolopes* sample in the [`trimming.sh`](trimming.sh) script.


### Read concatenation

To increase read depth, we concatenated two trimmed samples: **200409** and **212493**. These were concatenated into a single sample referred to as **409493**, which was used for downstream mapping, using the following command:

```bash
zcat 200409_R1_trimmed.fastq.gz 212493_R1_trimmed.fastq.gz | gzip -c > 409493_R1.fastq.gz
zcat 200409_R2_trimmed.fastq.gz 212493_R2_trimmed.fastq.gz | gzip -c > 409493_R2.fastq.gz
```

Note: The same approach was also applied to *S. officinalis*, where samples 320992 and 327270 were concatenated into sample 992270. However, only the *E. scolopes* processing is documented here.

### Create Bowtie2 index of reference genome for HiC-Pro mapping

```bash
module load bowtie
bowtie2-build Lachesis_assembly.fasta Lachesis_assembly
```

### Make file of chromosome sizes for HiC-Pro mapping

This was run using the script [`getSize.pl`](getSize.pl).

```bash
perl getSize.pl Lachesis_assembly.fasta > Lachesis_assembly_chrsize.txt
```

### Mapping with HiC-Pro

The concatenated sample **409493** was mapped using **HiC-Pro**.

This step was run using the script [`hicpro_mapping.sh`](hicpro_mapping.sh).  
The HiC-Pro configuration used is provided in [`config-hicpro_eupsc.txt`](config-hicpro_eupsc.txt).

### Remove interchromosomal interactions from .allValidPairs file (from the HiC-Pro output)

```bash
awk '$2 == $5' 409493.allValidPairs > 409493_intrachrom.allValidPairs
```

### Remove unplaced scaffolds from .allValidPairs file

```bash
awk '$2 == $5 && $2 ~ /Lachesis/' 409493_intrachrom.allValidPairs > 409493_intrachrom_rm_scaffolds.allValidPairs
```

### Generate the .hic file and apply KR normalisation

The filtered data was formatted into a Juicebox .hic file and KR normalised using the script [`run_hicpro2juicebox.sh`](run_hicpro2juicebox.sh).
The `hicpro2juicebox.sh` script was obtained from the [HiC-Pro GitHub repo](https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh).

### Dump .hic matrix at 100 kb resolution with KR normalisation for every chromosome for downstream analyses

```bash
#!/bin/bash

# Define the command and its components
juicer_tools="/scratch/molevo/thea/juicer_tools_1.22.01.jar"
input_file="409493_intrachrom.allValidPairs.hic"
output_directory="/path/to/output/directory/"
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
    output_file="${output_directory}409493_${chromosome}_KR_${bin_size}.dumped.hic"
    java -jar "$juicer_tools" dump observed KR "$input_file" "$chromosome" "$chromosome" BP "$bin_size" > "$output_file"
    echo "Processed chromosome $chromosome"
done
```

### Remove lines from dumped matrices that start with "WARN" and concatenate all matrices of single chromosomes into one file

A chromosome identifier is added to the first column of each line to track the source.

```bash

#!/bin/bash

# Directory containing your dumped files
input_dir="409493"
output_file="409493_intrachrom_allchrs_KR_100000.dumped.hic.txt"

# Remove the output file if it already exists
rm -f $output_file

# Loop through each file in the directory
for file in "$input_dir"*_KR_50000.dumped.hic; do
    # Extract chromosome name from the filename
    chromosome=$(basename "$file" | cut -d'_' -f2)
    
    # Process the file: remove WARN lines and add chromosome name as the first column
    grep -v '^WARN' "$file" | awk -v chr="$chromosome" '{print chr, $0}' >> "$output_file"
done

echo "Concatenated file created: $output_file"


```






