# Compare the Nosology against the genes_to_disease.txt file provided by the HPO Team
# Entries without a match need further investigation

library(tidyverse)
library(here)


# read the nosology dataframe with HGNC symbols created in 1c
nosology <- read_rds(here("data/Nosology_2023_HGNC.rds")) |>
  add_column(AZ_Comments = "", .after = "NOS_OMIM") |>
  add_column(PMID = "", .after = "NOS_OMIM") # for key PubMed articles, especially if not in OMIM

# create backup before making changes
write_rds(nosology, here("data/nosology_old.rds"))

# read genes_to_disease.txt file, downloaded and renamed from here:
# https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-12-12/genes_to_disease.txt

genes_to_disease <- read_tsv(here("raw_data/genes_to_disease_2024-12-12.tsv")) |>
  select(c("gene_symbol", disease_id)) |>
  # remove 'OMIM:' prefix
  mutate(disease_id = str_replace(disease_id, "OMIM:", ""))


# create an antijoin to see which genes in the nosology do not have a match in mim2gene
nosology_antijoin_by_gene <- nosology |>
  filter(!is.na(HGNC_symbol)) |>
  anti_join(genes_to_disease, by = c("HGNC_symbol" = "gene_symbol"))


# fix discrepancies found in nosology_antijoin_by_gene
# we can do this with rows_update as shown below
# the advantage is that it documents clearly which changes have been made
# because of this, I can make the changes directly in the nosology object
# the downside is that it is a lot of typing, we have to specify the column to update each time,
# and I also have to switch back and forth between my view of nosology_antijoin_by_gene and this script
# to write down these changes
# see further down for how to do this with DataEditR instead
nosology <- nosology %>%
  rows_update(tibble(NOS_ID = "NOS 06-0130", PMID = "30796325", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 10-0310", NOS_OMIM = "", PMID = "33200460; 38585547; 38647386", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 13-0060", NOS_OMIM = "", PMID = "35266227", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 13-0070", NOS_OMIM = "", PMID = "35266227", AZ_Comments = "few cases, not in OMIM or MONDO"))
# add more rows here as needed

# rerun antijoin_byGene to show that these changes have been incorporated
nosology_antijoin_by_gene <- nosology |>
  filter(!is.na(HGNC_symbol)) |>
  anti_join(genes_to_disease, by = c("HGNC_symbol" = "gene_symbol"))

# lets do some more updates with DataEditR
# save the updated nosology object to a different dataframe so we can compare it to the old one
library(DataEditR)
nosology_antijoin_by_gene_updated <- data_edit(nosology_antijoin_by_gene)

# show differences between nosology and nosology_updated
library(compareDF)
library(htmlTable)

nos_dif <- compare_df(nosology_antijoin_by_gene_updated, nosology_antijoin_by_gene)
create_output_table(nos_dif)

# or even better with waldo (package by Hadley Wickham)
library(waldo)
compare(nosology_antijoin_by_gene_updated, nosology_antijoin_by_gene)

# save the edits in case we need them
waldo_nos_dif <- compare(nosology_antijoin_by_gene_updated, nosology_antijoin_by_gene)
# update nosology_updated_antijoin_byGene with the edits from nosology_updated_antijoin_byGene2
nosology_updated <- rows_update(nosology, nosology_antijoin_by_gene_updated)

# save the updated nosology object
write_rds(nosology_updated, here("data/nosology.rds"))

# things to add:
# some cases of kyphomelic dysplasia and mesomelic dysplasia Kozlowski-Reardon
# are caused by biallelic PLOD2 mutations, see PMID: 29178448

#------------------------ work on this later ------------------------

# split rows with multiple OMIM IDs or multiple MOIs
nosology_split <- nosology |>
  separate_rows(NOS_OMIM, sep = ",\\s*") |> # split at comma plus optional white space
  separate_rows(NOS_MOI, sep = ",\\s*") # split at comma plus optional white space

# find rows in nosology that match on OMIM ID to genes_to_disease
# but do not match on HGNC symbol
nosology_antijoin_byOMIM <- nosology_split |>
  # filter(!is.na(HGNC_symbol)) |>
  anti_join(genes_to_disease, by = c("NOS_OMIM" = "disease_id"))


# save the updated nosology object
