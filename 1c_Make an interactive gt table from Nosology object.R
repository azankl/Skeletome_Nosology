#trying to make an interactive gt table
library(tidyverse)
library(gt)
library(here)


Nosology_2023 <- read_rds(here("data/Nosology_2023.rds")) |>
  select(NOS_ID, NOS_Name, NOS_MOI, NOS_Gene, NOS_OMIM, NOS_Notes, NOS_Group_Name)

nosology_tbl <- gt (Nosology_2023, groupname_col = "NOS_Group_Name") |>
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
