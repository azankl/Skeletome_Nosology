# Compare the Nosology against the genes_to_disease.txt file provided by the HPO Team
# Entries without a match need further investigation

library(tidyverse)
library(here)


# read the nosology dataframe with HGNC symbols created in 1c
nosology <- read_rds(here("data/nosology_with_HGNC.rds")) |>
  # separete rows with multiple OMIM IDs
  separate_rows(NOS_OMIM, sep=",\\s+") |> # split at comma plus optional white space\
  
  # add a new SkelDys-ID column as the NOS numbers are a bad idea
  # (they are linked to the current grouping, which might change in the future)
  # the SkelDys-ID is unique and random, like an OMIM_ID
  # cant just use OMIM_ID as some diseases have no OMIM_ID
  # main use of SkelDys_ID is to have an unique identifier as disease names and even gene names change
  # NOTE: decided that AD/AR versions of the same disease should have the same SkelDys_ID
  # this is how its done in OMIM and MONDO
  # but they should be separate rows in the nosology, so we can record the allelic requirements
  # functional effects, PMIDs etc separately for AD and AR
  # therefore:assign SkelDys_ID column and values first,
  mutate (id = row_number(), .before = NOS_ID) |>
  
  # then split rows with multiple MOI (they will get the same SkelDys_ID)
  separate_rows(NOS_MOI, sep=",\\s+") |> # split at comma plus optional white space
  
  # add column for G2P Mutation Consequence: https://www.ebi.ac.uk/gene2phenotype/updates_to_our_terms?ref=blog.opentargets.org
  # prefer the old G2P labels to the new GenCC ones
  # can create the new labels from the old ones if needed
  mutate (mutation_consequence = "", .after = PMID)


# read genes_to_disease.txt file from github
HPO <- read_tsv("https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt") |>
  select(c("gene_symbol", disease_id)) 

# keep OMIM only (full list includes ORPHA as well)
HPO_OMIM <- HPO |>
  filter(str_detect(disease_id, "OMIM")) |>
  mutate(disease_id = str_replace(disease_id, "OMIM:", "")) # remove 'OMIM:' prefix

# find where nosology and HPO match by OMIM ID and gene symbol
# (since nosology sometimes assigns OMIM IDs to a different disease than OMIM)
NOS_HPO_semijoin <- nosology |>
  semi_join(HPO_OMIM, by = join_by(NOS_OMIM == disease_id, HGNC_symbol == gene_symbol))

# create an antijoin to see which genes in the nosology do not have a match in HPO
NOS_HPO_antijoin <- nosology |>
 # filter(!is.na(HGNC_symbol)) |>
  anti_join(HPO_OMIM, by = join_by(NOS_OMIM == disease_id, HGNC_symbol == gene_symbol))

# describe findings here
nosology <- nosology |>
  rows_update(tibble(id = 12, AZ_Comments = "missing in HPO, opened HPO Github issue")) |>
  rows_update(tibble(id = 35, AZ_Comments = "fixed wrong OMIM ID", NOS_OMIM = "620269")) |>
  rows_update(tibble(id = 52, AZ_Comments = "not in OMIM, opened HPO Github issue")) |>
  rows_update(tibble(id = 63, AZ_Comments = "not in OMIM", PMID = "30796325")) |>
  rows_update(tibble(id = 68, AZ_Comments = "Rolland-Desbuqoius, gene missing in OMIM", PMID = "38424183", mutation_consequence = "all_missense")) |>
  rows_update(tibble(id = 67, AZ_Comments = "Silverman-Handmaker, added mutation consequence", PMID = "11279527", mutation_consequence = "absent_gene_product")) |>
  rows_update(tibble(id = 69, AZ_Comments = "added mutation consequence", PMID = "16927315", mutation_consequence = "hypomorphic")) |>
  rows_update(tibble(id = 88, AZ_Comments = "fixed wrong OMIM ID", NOS_OMIM = "620022")) |>
  #rows_delete(tibble(id = 91)) |> # remove duplicate entry for Short rib–polydactyly syndrome (SRPS), DYNC2H1-related with wrong OMIM ID 263520
  rows_update(tibble(id = 98, AZ_Comments = "gene missing in OMIM", PMID = "28370949")) |>
  rows_update(tibble(id = 99, AZ_Comments = "fixed wrong OMIM ID", NOS_OMIM = "614376")) |>
  rows_update(tibble(id = c(101,122), AZ_Comments = "removed wrong OMIM ID, referred to gene, no phenotype ID in OMIM, see PMID", NOS_OMIM = "", PMID = "29068549")) |>
  rows_update(tibble(id = c(104,129), AZ_Comments = "OMIM ID for similar phenotypes", NOS_OMIM = "617088")) |>
  rows_update(tibble(id = 105, AZ_Comments = "OMIM ID for similar phenotypes, note new gene name", NOS_OMIM = "615633" )) |>
  rows_update(tibble(id = 107, AZ_Comments = "OMIM ID for similar phenotypes, note new gene name", NOS_OMIM = "615503" )) |>
  rows_update(tibble(id = 111, AZ_Comments = "OMIM ID for similar phenotype, gene missing in OMIM", NOS_OMIM = "269860", PMID = "28370949")) |>
  rows_update(tibble(id = 114, AZ_Comments = "fixed typo in OMIM ID", NOS_OMIM = "611263")) |>
  rows_update(tibble(id = 121, AZ_Comments = "removed wrong OMIM ID, referred to gene, no phenotype ID in OMIM, see PMID", NOS_OMIM = "", PMID = "33200460" )) |>
  rows_update(tibble(id = 125, AZ_Comments = "removed wrong OMIM ID, referred to Mohr syndrome, no phenotype ID in OMIM, see PMID", NOS_OMIM ="",  PMID = "28123176")) |>

  
  
  # rows_update(tibble(id = , AZ_Comments = "", NOS_OMIM = "", PMID ="" )) |>
  
                    
# add SRPS type II caused by biallelic FUZ mutations, one case in PMID 29068549 (Cohn et al), not in OMIM, not in Nosology
  
  

# ------ stuff below is old, check before using ------

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
