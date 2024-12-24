# This script performs some clean up on the Nosology_2023 dataframe,
# which was created in '1a_Convert AJMG table to Nosology object.R'

nosology <- read_rds(here("data/Nosology_2023.rds"))

# The ISDS Nosology sometimes has multiple DMIM values in the DMIM field.
# split such entries into individual rows
# the regex splits at comma plus optional white space
nosology <- nosology %>%
  separate_rows(NOS_OMIM, sep=",\\s+")

# then fix errors:
# need to decide when these errors should be fixed, before or after
# comparing and joining with other data like G2P and PanelApp

# there are multiple ways to fix errors, I will show a few examples

# Fix errors with rows_upate
# change NOS_Gene IFT40 to IFT140 for NOS 10-0470 (typo in the original table)
nosology <- nosology %>%
  rows_update(tibble (NOS_ID = "NOS 10-0470", NOS_Gene = "IFT140"))

#could also use mutate with replace to correct errors in the table
#as described here: https://stackoverflow.com/questions/36924911/how-to-assign-a-value-to-a-data-frame-filtered-by-dplyr

#Fix errors:
#NOS 40-0170 is LADD syndrome 3 (FGF10 associated), but the nosology gives the DMIM number for LADD syndrome 1 (FGFR2 associated)
#NOS 40-0160 is LADD syndrome 2 (FGFR3 associated), but the nosology gives the DMIM number for LADD syndrome 1 (FGFR2 associated)

nosology <- nosology %>%
  mutate(NOS_OMIM = replace(NOS_OMIM, NOS_ID=="NOS 40-0170", "620193")) %>%
  mutate(NOS_OMIM = replace(NOS_OMIM, NOS_ID=="NOS 40-0160", "620192"))

#alternatively, replace directly using base R syntax as described here: https://sparkbyexamples.com/r-programming/replace-values-in-r/
#nosology$NOS_OMIM[nosology$NOS_ID=="NOS 40-0170"] <- "620193"

#More errors I found in Nosology:
#NOS 10-0480 Cranioectodermal dysplasia Levin-Sensenbrenner is listed as DMIM 614009, should be DMIM 614099.
nosology <- nosology %>%
  mutate(NOS_OMIM = replace(NOS_OMIM, NOS_ID=="NOS 10-0480", "614099"))

#NOS 34-0190 (Non-syndromic midline craniosynostosis, RUNX2-related) is more a risk factor than a Mendelian disorder, see PMID: 32360898, should be removed
nosology <- nosology %>%
  filter(NOS_ID != "NOS 34-0190")


#Export cleaned nosology
write_rds (nosology, here("data/nosology_cleaned.rds"))
#write_csv(nosology, here("data/nosology_cleaned.csv"))
#write_tsv(nosology, here("data/nosology_cleaned.tsv"))