#trying to make an interactive gt table
#(consider using reactable and reactablefmtr instead, has a few more options, but less nice syntax)
library(tidyverse)
library(gt)
library(here)

Nosology_2023 <- read_rds(here("data/Nosology_2023.rds")) |>
  select(NOS_ID, NOS_Name, NOS_MOI, NOS_Gene, NOS_OMIM, NOS_Notes, NOS_Group_Name)

nosology_tbl <- gt(Nosology_2023, groupname_col = "NOS_Group_Name") |>
  opt_interactive(
    active = TRUE,
    use_search = TRUE,
    use_filters = TRUE,
    use_compact_mode = TRUE,
    use_resizers = TRUE,
    use_highlight = TRUE,
    use_pagination = TRUE
  ) |>
  cols_label(
    NOS_ID = "Nosology ID",
    NOS_Name = "Disorder Name",
    NOS_MOI = "Mode of Inheritance",
    NOS_Gene = "Gene",
    NOS_OMIM = "OMIM",
    NOS_Notes = "Notes"
  )
#show the table
nosology_tbl


#highlight changes made to the table
library(waldo)
nosology <- Nosology_2023
nosology_updated <- nosology %>%
  rows_update(tibble(NOS_ID = "NOS 06-0130", PMID = "30796325", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 10-0310", NOS_OMIM = "", PMID = "33200460; 38585547; 38647386", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 13-0060", NOS_OMIM = "", PMID = "35266227", AZ_Comments = "few cases, not in OMIM or MONDO")) |>
  rows_update(tibble(NOS_ID = "NOS 13-0070", NOS_OMIM = "", PMID = "35266227", AZ_Comments = "few cases, not in OMIM or MONDO"))
nos_dif <- compare_df(nosology_updated, nosology)

nos_dif_table <-gt(nos_dif) |>
  tab_header(
    title = "Differences between Nosology and Nosology_updated"
  )
