# Preprocessing and mapping of Micro-C reads

This folder contains the commands used to process raw Micro-C reads. These steps include trimming, merging libraries for higher coverage, and mapping with HiC-Pro, and generating a Juicebox `.hic` file. *E. scolopes* stage 29 samples are used as examples throughout this folder.  
Where applicable, *S. officinalis* is used instead, for example, in the downsampling step, which was not necessary for *E. scolopes*.


## Downsampling
We downsampled relevant tissue-specific Micro-C samples to approximately 350 million reads to normalise for read depth across tissues.  
For example, for the *S. officinalis* eye tissue, this was applied to read 1 using the command below, where `0.8559` was the downsampling fraction used.

```bash
module load seqtk
seqtk sample -s100 327271_R1.fastq.gz 0.8559 | gzip > /lisc/scratch/molevo/thea/sepof_microc/microc3/trim50_data/ds_EY/327271/327271_R1_ds.fastq.gz
```

## Trimming

Reads were trimmed to 50 bp using **Trimmomatic** with the `CROP:50` setting.  
This step is documented for one stage 29 *E. scolopes* sample in the [`trimming.sh`](trimming.sh) script.


## Read concatenation – *E. scolopes*

To increase read depth, we concatenated two trimmed samples:

- **200409**
- **212493**

These were merged into a single sample referred to as **409493**, which was used for downstream mapping.

The following command was run directly in the terminal:

```bash
zcat 200409_R1_trimmed.fastq.gz 212493_R1_trimmed.fastq.gz | gzip -c > 409493_R1.fastq.gz
zcat 200409_R2_trimmed.fastq.gz 212493_R2_trimmed.fastq.gz | gzip -c > 409493_R2.fastq.gz
```

Note: The same approach was also applied to Sepia officinalis, where samples 320992 and 327270 were concatenated into sample 992270. However, only the E. scolopes processing is documented here.

## Create Bowtie2 index of reference genome for HiC-Pro step

```bash
module load bowtie
bowtie2-build Lachesis_assembly.fasta Lachesis_assembly
```

## Make file of chromosome sizes for HiC-Pro step. The script getSize.pl is included in this folder

```bash
perl getSize.pl Lachesis_assembly.fasta > Lachesis_assembly_chrsize.txt
```

## Mapping with HiC-Pro

The merged sample **409493** was mapped using [HiC-Pro].

This step was run using the script [`hicpro_mapping.sh`](hicpro_mapping.sh).  
The HiC-Pro configuration used is provided in [`config-hicpro_scolopes1_all.txt`](config-hicpro_scolopes1_all.txt).

## Remove interchromosomal interactions from .allValidPairs file (from the HiC-Pro output)

```bash
awk '$2 == $5' 409493.allValidPairs > 409493_intrachrom.allValidPairs
```

## Remove unplaced scaffolds from .allValidPairs file

```bash
awk '$2 == $5 && $2 ~ /Lachesis/' 409493_intrachrom.allValidPairs > 409493_intrachrom_rm_scaffolds.allValidPairs
```

## Generate the .hic file and apply KR normalisation

The filtered data was formatted into a Juicebox .hic file and KR normalised using the script [`run_hicpro2juicebox.sh`](run_hicpro2juicebox.sh).
The # hicpro2juicebox.sh script was obtained from:
# https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh





