# Interacting gene pair analysis

This folder documents the interacting gene pair analyses. Includes gene-bin mapping, and ortholog checks between species. All steps are demonstrated using the *E. scolopes* sample 403493 at 100 kb resolution. Additional commands are provided at the start for *S. officinalis* to complete orthology comparisons, since at the time of writing this paper, no gene annotation was available for the *S. officinalis* reference genome.


**Note on input files of orthologous genes:**
BLASTP (v2.16.0) was run with the parameters: -evalue 1E-2 -max_target_seqs 1 -outfmt 6 to classify reciprocal best hit orthologs between *E. scolopes*, O. bimaculoides* and *P. maximus*. *E. scolopes* was used as the reference (database) in all cases, except for identifying *P. maximus* orthologs of *O. bimaculoides* genes, where *O. bimaculoides* was used as the database. Text files of 1:1 orthologs in both directions were generated: **EUPgeneOBI.txt**, **EUPgenePEC.txt**, **OBIgeneEUP.txt**, and **OBIgenePEC.txt**, where: EUP = *E. scolopes*, OBI = *O. bimaculoides*, PEC = *P. maximus*. If the species abbreviation appears first in the filename, this indicates that the species' genes are in the first column of the file. For example, in EUPgeneOBI.txt, *E. scolopes* genes are in the first column, and their best reciprocal orthologs in *O. bimaculoides* are in the second. This directionality matters in the downstream analyses.

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

Count number of genes left in file

```bash
wc -l  sanger_sepof_eup_prot_best_hits_rm_scaff.gff
24431  sanger_sepof_eup_prot_best_hits_rm_scaff.gff #There are 24431 *. scolopes* proteins that mapped to the *S. officinalis* genome orthologs that we considered orthologs.
```
The resulting **sanger_sepof_eup_prot_best_hits_rm_scaff.gff** file was then converted into bed format, with chromosome names being from the *S. officinalis* genome file, and gene names being identical to the *E. scolopes* gene names that were the best hits mapped to the *S. officinalis* genome. 

Orthology files were also made for downstream analyses: **EUPgeneSOF.txt**, **SOFgeneEUP.txt**, although gene names were identical in both columns, the files were formatted this way to be compatible with downstream scripts.


### Match gene pairs to bins

We used a custom Python script  [`check_gene_in_bin_dump_all.py`](check_gene_in_bin_dump_all.py) to identify gene pairs in bins pairs in dumped matrices. The input file is **409493_intrachrom_allchrs_KR_100000.dumped.hic.txt** is from the end of step 1 which creates a dumped matrix of all intrachromosomal interactions. This script includes any genes overlapping with the bins in the output files:

```bash
python3 check_gene_in_bin_dump_all.py eupsc.bed 409493_intrachrom_allchrs_KR_100000.dumped.hic.txt 100000
```
Where **eupsc.bed** is the species BED file containing four columns: chromosome, gene start, gene end, and gene name (in that order), with no additional features. This produced an output like:

```bash
Number of bins with genes in = 26197
Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
```

This command was repeated for *S. officinalis* at 100kb resolution and *O. bimaculoides* at 50 kb resolution.

### Check for orthologous genes

The script [`check_orthos_dump.py`](check_orthos_dump.py) was used to assess whether gene pairs in bin interactions have orthologs in other species:

```bash
python3 check_orthos_dump.py EUPgeneOBI.txt 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
```
Where the file EUPgeneOBI.txt is a file with *E. scolopes* genes in the first column, and their *O. bimaculoides* orthologs in the second column.  The species whose orthologs appear in the first column of the ortholog file must correspond to the species from which the interaction matrix was generated.

Example output:

```bash
Number of reciprocal best hit orthologs for EUP and OBI = 12002
Number of EUP interactions with at least one ortholog in OBI = 1569457 #Note this file likely contains duplicate entries that have not yet been removed. As such, the reported number is not biologically meaningful at this stage
Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq_OBIorthos.txt
```

Commands were repeated for the following species comparisons:

*E. scolopes* ↔ *O. bimaculoides*

*E. scolopes* ↔ *P. maximus*

*O. bimaculoides* ↔ *E. scolopes*

*O. bimaculoides* ↔ *P. maximus*

Note: Since *S. officinalis*–*E. scolopes* orthologs were named identically to *E. scolopes* genes, it was not necessary to extract ortholog names separately in this case. *S. officinalis* interactions were added downstream.

