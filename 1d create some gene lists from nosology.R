# now that the nosology has a cleaned up HGNC symbol column
# we can run some stats on the data
# should probably turn this into a Quarto document


library(tidyverse)
library(here)

# read nosology with cleaned up HGNC symbols column
nosology <- read_rds(here("data/nosology_with_HGNC.rds"))

# find entries without a gene in HGNC_symbol column
# should be the sama as nosology_no_gene above
new_nosology_no_gene <- nosology %>%
  filter(is.na(HGNC_symbol))
# it is: 34 entries with no gene
new_nosology_with_gene <- nosology %>%
  filter(!is.na(HGNC_symbol))

# find entries with a single HGNC_symbol
single_gene_symbol_pattern <- "^[A-Z0-9-]+$"
new_nosology_single_gene <- new_nosology_with_gene %>%
  filter(grepl(single_gene_symbol_pattern, HGNC_symbol))
# this is a good list for PanelAPP etc as it can be easily parsed

# find entries that done match to a single HGNC_symbol
# those should be the del/dup/reg variants
new_nosology_complex_gene <- new_nosology_with_gene %>%
  filter(!grepl(single_gene_symbol_pattern, HGNC_symbol))

# find unique values in HGNC_symbol column using dyplr
new_nosology_unique_single_genes <- new_nosology_single_gene %>%
  distinct(HGNC_symbol) |>
  # sort alphabetically
  arrange(HGNC_symbol)

# extract gene symbols from HGNC_symbol column from new_nosology_complex_gene
# requires a more complex regex pattern
complex_gene_symbol_pattern <- "(?=[A-Z])[A-Z0-9-]*" # at least one uppercase letter
new_nosology_unique_complex_genes <- new_nosology_complex_gene|>
  select (HGNC_symbol) |>
  str_extract_all(complex_gene_symbol_pattern) |>
  unlist() |>
  as_tibble() |>
  rename(HGNC_symbol = value) |>
  distinct(HGNC_symbol) |>
  arrange(HGNC_symbol)


# combine new_nosology_unique_genes and new_nosology_nonstandard_gene_symbols
# to get a list of all gene symbols in the nosology
new_nosology_all_genes <- new_nosology_unique_single_genes |>
  rows_append(new_nosology_unique_complex_genes) |>
  distinct(HGNC_symbol) |>
  arrange(HGNC_symbol)