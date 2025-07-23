# Interacting gene pair analysis

This folder documents the interacting gene pair analyses. Includes  gene-bin mapping, and ortholog checks between species. All steps are demonstrated using the *E. scolopes* sample 403493 at 100 kb resolution. Additional commands are provided for *S. officinalis* to complete orthology comparisons as at the time of writing this paper, no gene annotation was available for *S. officinalis*.

Note: orthology was assessed using reciprocal best BLAST hits for *E. scolopes*, *O. bimaculoides* and *P. maximus* using BLASTP 2.16.0 with the parameters -evalue 1E-2 -max_target_seqs 1 -outfmt 6.

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
24431  sanger_sepof_eup_prot_best_hits_rm_scaff.gff #There are 24431 orthologs
```
The resulting **sanger_sepof_eup_prot_best_hits_rm_scaff.gff** file was then converted into bed format, with chromosome names being from the *S. officinalis* genome file, and gene names being identical to the *E. scolopes* gene names that were the best hits mapped to the *S. officinalis* genome.


### Match gene pairs to bins

We used a custom Python script  [`check_gene_in_bin_dump_all.py`](check_gene_in_bin_dump_all.py) to identify gene pairs in bins pairs in dumped matrices. This script includes any genes overlapping with the bins in the output files:

```bash
python3 check_gene_in_bin_dump_all.py eupsc.bed 409493_intrachrom_allchrs_KR_100000.dumped.hic.txt 100000
```

This produced an output like:

```bash
Number of bins with genes in = 26197
Output written to: 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt
```


3. Check for orthologous genes

A second script (check_orthos_dump.py) was used to assess whether gene pairs in bin interactions have orthologs in other species:

Example command:

python3 check_orthos_dump.py EUPgeneOBI.txt 409493_intrachrom_allchrs_KR_100000.dumped.hic_all_genes_int_freq.txt

Output:

Number of reciprocal best hit orthologs for EUP and OBI = 12002
Number of EUP interactions with at least one ortholog in OBI = 1569457

Commands were repeated for all relevant species comparisons:

E. scolopes ↔ O. bimaculoides

E. scolopes ↔ Pecten maximus

O. bimaculoides ↔ E. scolopes

O. bimaculoides ↔ P. maximus

Notes

The number of genes in bin interactions scales with genome size and resolution.

E. scolopes genome is roughly half the size of Octopus, justifying the use of different resolutions (e.g., 100kb in EUP ≈ 50kb in OBI).

The full resolution comparison is retained for completeness, but 50kb–25kb comparisons may be the most reliable.

Scripts Referenced

check_gene_in_bin_dump_all.py

check_orthos_dump.py

mahMbh.pl

These are included in this folder or referenced in the scripts/ subdirectory.

