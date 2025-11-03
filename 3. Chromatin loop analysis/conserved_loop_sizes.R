# Get loop sizes for conserved loops only (loop bin2 end - loop bin1 start)

# Libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Parse "<chrom>:<start>-<end>" into chrom + start
parse_loop_id <- function(x) {
  tibble(raw = x) %>%
    mutate(
      chrom = str_extract(raw, "^[^:]+"),
      start = str_extract(raw, "(?<=:)\\d+"),
      start = suppressWarnings(as.numeric(start))
    ) %>%
    select(chrom, start)
}

# 1) Read loop-size tables
eupsc <- read_tsv("eupsc_29cat_50k+100k.loopsize.tsv.genes_rm_dups",
                  col_types = cols())
sepof <- read_tsv("sepof_wt_50k+100k.loopsize.tsv.genes_rm_dups",
                  col_types = cols())
octbi <- read_tsv("octbi_wt_50k+100k.loopsize.tsv.genes_rm_dups",
                  col_types = cols())

# Standardize the columns we need
eupsc_clean <- eupsc %>%
  transmute(chromosome,
            bin1_start = as.numeric(`bin1_start (bp)`),
            loop_size_eupsc = loop_size)

sepof_clean <- sepof %>%
  transmute(chromosome,
            bin1_start = as.numeric(`bin1_start (bp)`),
            loop_size_sepof = loop_size)

octbi_clean <- octbi %>%
  transmute(chromosome,
            bin1_start = as.numeric(`bin1_start (bp)`),
            loop_size_octbi = loop_size)

# 2) eupsc ↔ sepof conserved loops
cons_es_sp <- read_tsv("eupsc_sepof_consv_loops_50k+100k.txt", col_types = cols())

cons_es_sp_parsed <- cons_es_sp %>%
  mutate(
    # parse IDs "chrom:start-end"
    eupsc_chr = parse_loop_id(eupsc_loop)$chrom,
    eupsc_bin1 = parse_loop_id(eupsc_loop)$start,
    sepof_chr = parse_loop_id(sepof_loop)$chrom,
    sepof_bin1 = parse_loop_id(sepof_loop)$start
  ) %>%
  left_join(eupsc_clean, by = c("eupsc_chr" = "chromosome", "eupsc_bin1" = "bin1_start")) %>%
  left_join(sepof_clean, by = c("sepof_chr" = "chromosome", "sepof_bin1" = "bin1_start"))

write_tsv(cons_es_sp_parsed,
          "/Users/users/Desktop/Micro-C/diff_loop_analysis/conserved_loops/eupsc_sepof_consv_loops_with_sizes.tsv")

# 3) eupsc ↔ octbi conserved loops
cons_es_ob <- read_tsv("eupsc_octbi_consv_loops_50k+100k.txt", col_types = cols())

cons_es_ob_parsed <- cons_es_ob %>%
  mutate(
    eupsc_chr = parse_loop_id(eupsc_loop)$chrom,
    eupsc_bin1 = parse_loop_id(eupsc_loop)$start,
    octbi_chr = parse_loop_id(octbi_loop)$chrom,
    octbi_bin1 = parse_loop_id(octbi_loop)$start
  ) %>%
  left_join(eupsc_clean, by = c("eupsc_chr" = "chromosome", "eupsc_bin1" = "bin1_start")) %>%
  left_join(octbi_clean, by = c("octbi_chr" = "chromosome", "octbi_bin1" = "bin1_start"))

write_tsv(cons_es_ob_parsed,
          "eupsc_octbi_consv_loops_with_sizes.tsv")

# 4) sepof ↔ octbi conserved loops
# (this file already has separate start/end columns, but we still join by chrom + start)
cons_sp_ob <- read_tsv("sepof_octbi_consv_loops_50k+100k.txt", col_types = cols())

# ensure numeric (in case they came as char)
cons_sp_ob2 <- cons_sp_ob %>%
  mutate(
    sepof_start = as.numeric(sepof_start),
    octbi_start = as.numeric(octbi_start)
  ) %>%
  left_join(sepof_clean, by = c("sepof_chrom" = "chromosome", "sepof_start" = "bin1_start")) %>%
  left_join(octbi_clean, by = c("octbi_chrom" = "chromosome", "octbi_start" = "bin1_start"))

write_tsv(cons_sp_ob2,
          "sepof_octbi_consv_loops_with_sizes.tsv")

# 5) Three-species file (eupsc + octbi + sepof in one)
# Assumes columns named exactly like your example: eupsc_loop, octbi_loop, sepof_loop
cons_all3 <- read_tsv("eupsc_octbi_consv_loops_50k+100k_with_sepof.txt", col_types = cols())

cons_all3_parsed <- cons_all3 %>%
  mutate(
    eupsc_chr = parse_loop_id(eupsc_loop)$chrom,
    eupsc_bin1 = parse_loop_id(eupsc_loop)$start,
    octbi_chr = parse_loop_id(octbi_loop)$chrom,
    octbi_bin1 = parse_loop_id(octbi_loop)$start,
    sepof_chr = parse_loop_id(sepof_loop)$chrom,
    sepof_bin1 = parse_loop_id(sepof_loop)$start
  ) %>%
  left_join(eupsc_clean, by = c("eupsc_chr" = "chromosome", "eupsc_bin1" = "bin1_start")) %>%
  left_join(octbi_clean, by = c("octbi_chr" = "chromosome", "octbi_bin1" = "bin1_start")) %>%
  left_join(sepof_clean, by = c("sepof_chr" = "chromosome", "sepof_bin1" = "bin1_start"))

write_tsv(cons_all3_parsed,
          "eupsc_octbi_sepof_consv_loops_with_sizes.tsv")

# 6) Quick diagnostics (optional)
diag <- tibble(
  dataset = c("eupsc-sepof", "eupsc-octbi", "sepof-octbi", "all-three"),
  n_rows = c(nrow(cons_es_sp_parsed), nrow(cons_es_ob_parsed), nrow(cons_sp_ob2), nrow(cons_all3_parsed)),
  n_missing_eupsc = c(sum(is.na(cons_es_sp_parsed$loop_size_eupsc)),
                      sum(is.na(cons_es_ob_parsed$loop_size_eupsc)),
                      NA_integer_,
                      sum(is.na(cons_all3_parsed$loop_size_eupsc))),
  n_missing_sepof = c(sum(is.na(cons_es_sp_parsed$loop_size_sepof)),
                      NA_integer_,
                      sum(is.na(cons_sp_ob2$loop_size_sepof)),
                      sum(is.na(cons_all3_parsed$loop_size_sepof))),
  n_missing_octbi = c(NA_integer_,
                      sum(is.na(cons_es_ob_parsed$loop_size_octbi)),
                      sum(is.na(cons_sp_ob2$loop_size_octbi)),
                      sum(is.na(cons_all3_parsed$loop_size_octbi)))
)
print(diag)

