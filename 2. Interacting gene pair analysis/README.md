# Interacting gene pair analysis

## Contents

- [Prepare files of gene pairs, orthologs, and extract their chromosomal distribution in *P. maximus*](#prepare-files-of-gene-pairs-orthologs-and-extract-their-chromosomal-distribution-in-p-maximus)
	- [Make orthology annotation file for *S. officinalis*](#make-orthology-annotation-file-for-s-officinalis)
	- [Match gene pairs to bins](#match-gene-pairs-to-bins)
	- [Check for orthologous genes](#check-for-orthologous-genes)
	- [Format interaction matrix files](#format-interaction-matrix-files)
	- [Merge interaction frequency values across species](#merge-interaction-frequency-values-across-species)
	- [Classify interactions by *P. maximus* chromosome status](#classify-interactions-by-pecten-maximus-chromosome-status)
	- [Add interaction frequency from *S. officinalis*](#add-interaction-frequency-from-s-officinalis)

This folder documents the interacting gene pair analyses. Includes gene-bin mapping, and ortholog checks between species. All steps are demonstrated using the *E. scolopes* sample 403493 at 100 kb resolution. Additional commands are provided for the species *S. officinalis* to classify orthologous genes and for downstream analyses, as no gene annotation was available for the *S. officinalis* reference genome at the time of writing this paper.

## Prepare files of gene pairs, orthologs, and extract their chromosomal distribution in **P. maximus**

**Notes on input files of orthologous genes:**

BLASTP (v2.16.0) was run with the parameters: -evalue 1E-2 -max_target_seqs 1 -outfmt 6 to classify reciprocal best hit orthologs between *E. scolopes*, O. bimaculoides* and *P. maximus*. *E. scolopes* was used as the reference (database) in all cases, except for identifying *P. maximus* orthologs of *O. bimaculoides* genes, where *O. bimaculoides* was used as the database. Text files of 1:1 orthologs in both directions were generated: **EUPgeneOBI.txt**, **EUPgenePEC.txt**, **OBIgeneEUP.txt**, and **OBIgenePEC.txt**, where: EUP = *E. scolopes*, OBI = *O. bimaculoides*, PEC = *P. maximus*. If the species abbreviation appears first in the filename, this indicates that the species' genes are in the first column of the file. For example, in EUPgeneOBI.txt, *E. scolopes* genes are in the first column, and their best reciprocal orthologs in *O. bimaculoides* are in the second. This directionality matters in the downstream analyses.

**Note on species bed files:**

Species BED files should contain four columns: chromosome, gene start, gene end, and gene name (in that order), with no additional features. 
Here, **eupsc.bed** is used as an example, corresponding to *E. scolopes*. 

### Make orthology annotation file for *S. officinalis*

*E. scolopes* proteins were mapped to the *S. officinalis* genome using miniprot. The command used creates a miniprot index at the same time. 
This step is documented in the [`index_map_miniprot.sh`](index_map_miniprot.sh) script.

The result of the miniprot is a mixture of gff and paf entries, so we filtered to keep just gff entries:

```bash
grep -v "^##PAF" sanger_sepof_eup_prot.gff > sanger_sepof_eup_prot_clean.gff
```
Extract only highest ranks hits:

```bash
grep -E 'mRNA' sanger_sepof_eup_prot_clean.gff | grep "Rank=1;" > sanger_sepof_eup_prot_best_hits.gff

```
Remove scaffolds:

```bash
grep OZ sanger_sepof_eup_prot_best_hits.gff > sanger_sepof_eup_prot_best_hits_rm_scaff.gff
```

Count number of genes left in file:

```bash
wc -l  sanger_sepof_eup_prot_best_hits_rm_scaff.gff
24431  sanger_sepof_eup_prot_best_hits_rm_scaff.gff # There are 24431 E. scolopes proteins that mapped to the S. officinalis genome that we considered orthologs.
```
The resulting **sanger_sepof_eup_prot_best_hits_rm_scaff.gff** file was then converted into bed format, with chromosome names being from the *S. officinalis* genome file, and gene names being identical to the *E. scolopes* gene names that were the best hits mapped to the *S. officinalis* genome. 

Orthology files were also made for downstream analyses: **EUPgeneSOF.txt**, **SOFgeneEUP.txt**, although gene names were identical in both columns, the files were formatted this way to be compatible with downstream scripts.


### Match gene pairs to bins

We used a custom Python script  [`check_gene_in_bin_dump_all.py`](check_gene_in_bin_dump_all.py) to identify gene pairs in bins pairs in dumped matrices. The input file is **409493_intrachrom_allchrs_KR_100000.dumped.hic.txt** is from the end of step 1 which creates a dumped matrix of all intrachromosomal interactions. This script includes any genes overlapping with the bins in the output files:

```bash
python3 check_gene_in_bin_dump_all.py eupsc.bed 409493_intrachrom_allchrs_KR_100000.dumped.hic.txt 100000
```
Where the species BED file matches the species whose interaction matrix it is. This produced an output like:

```bash
Number of bins with genes in = 26197
Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
```

This command was repeated for *S. officinalis* at 100kb resolution and *O. bimaculoides* at 50 kb resolution.

### Check for orthologous genes

The script [`check_orthos_dump.py`](check_orthos_dump.py) was then used to assess whether gene pairs in bin interactions have orthologs in other species. The species whose orthologs appear in the first column of the ortholog file (and whose abbreviation appears first in the filename) must correspond to the species from which the interaction matrix was generated.

```bash
python3 check_orthos_dump.py EUPgeneOBI.txt 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
# Example output:
# Number of reciprocal best hit orthologs for EUP and OBI = 12002 # Where EUP is E. scolopes and OBI is O. bimaculoides
# Number of EUP interactions with at least one ortholog in OBI = 1569457 # Note this file likely contains duplicate entries that have not yet been removed. As such, the reported number is not biologically meaningful at this stage
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos.txt
```

Commands were repeated for the following species comparisons:

*E. scolopes* ↔ *O. bimaculoides*

*E. scolopes* ↔ *P. maximus*

*O. bimaculoides* ↔ *E. scolopes*

*O. bimaculoides* ↔ *P. maximus*

Note: Since *S. officinalis*–*E. scolopes* orthologs were named identically to *E. scolopes* genes, it was not necessary to extract ortholog names separately in this case. *S. officinalis* interactions were added downstream.


### Format interaction matrix files

Before comparing interaction frequencies across species, run the [`format_interactions.py`](format_interactions.py) script to:
- Remove NA values
- Standardise pairwise gene interactions so that each row contains exactly one pair of interacting genes (i.e. some bins contain many genes, so all pairwise combinations are extracted)
- Remove swapped-bin duplicates (i.e. if bin1-bin2 and bin2-bin1 are both present, keep only one)

```bash
# Format E. scolopes (100 kb) with O. bimaculoides orthologs
python3 format_interactions.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos.txt
# Example output: 
# 2,565,954 interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos_form.txt

# Format E. scolopes (100 kb) with P. maximus orthologs
python3 format_interactions.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_PECorthos.txt
# Example output: 
# 2,426,852 interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_PECorthos_form.txt
```

### Merge interaction frequency values across species

Use [`eup_vs_obi_int_freq_form.py`](eup_vs_obi_int_freq_form.py) to merge gene-pair interaction frequencies between species, based on *E. scolopes* bin pairs at 100 kb resolution, handling cases where gene order differs between files.

```bash
python3 eup_vs_obi_int_freq_form.py 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos_form.txt \
                                    212489_intrachrom_allchrs_KR_50000.dumped.hic_all_genes_int_freq_EUPorthos_form.txt \
                                    EUPgeneOBI.txt
# Example output: 
# 3,185,945 matched gene-pair interactions
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq.txt
```
**Note:** the resulting file still contains duplicate gene pairs with differing interaction frequencies, along with all possible gene pair combinations across the two species. This contributes to the large number of interactions and file size. These redundancies will be resolved in later steps.


### Classify interactions by ***Pecten maximus*** chromosome status

This step uses [`syteny_by_topology_interactions.py`](syteny_by_topology_interactions.py) to classify gene-pair interactions as:
- Gene pair on the same *P. maximus* chromosome
- Gene pair on different *P. maximus* chromosomes

```bash
python3 synteny_by_topology_interactions.py EUPgenePEC.txt pmax2.bed \
       409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq.txt
# Example output: 
# Interactions on same P. maximus chromosomes: 917,352
# Interactions on different P. maximus chromosomes: 1,493,463
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_MACIs_by_topology.txt
```

### Add interaction frequency from *S. officinalis*

This step can be done without an ortholog file, since *S. officinalis* orthologs were named the same as *E. scolopes* genes. The script [`add_sof_interactions.py`](add_sof_interactions.py) was used.

```bash
python3 add_sof_interactions.py 992270_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt \
      409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_synteny_by_topology.txt
# Output: 
# S. officinalis interactions: 4,355,779
# Gene pairs with S. officinalis orthologs: 1,056,029
# Appended interactions: 4,772,919
# Output written to: 409493_intrachrom_allchrs_KR_100000_EUPvs50000_OBI_int_freq_PEC_chr_status_with_sof.txt
```
**Note**:  Swapped-bin duplicates were already removed from the second input file and so will not appear in the output of this command. However, the resulting file still contains duplicate gene pairs with differing interaction frequencies, along with all possible gene pair combinations across the three species. This contributes to the large file size. These redundancies will be resolved in later steps.



