# Add HGNC gene symbols to the nosology
# in the process, we also explore some issues with the NOS_Gene column

library(tidyverse)
library(tidyxl)
library(here)
library(glue)
library(DataEditR)
library(forgts)
library(waldo)
library(compareDF)
library(diffobj)
library(limma) # needed for alias2SymbolTable function
# org.Hs.eg.db needs to be installed for limma, make sure its up-to-date


# read the original nosology dataframe (uncleaned) and
# add HGNC gene symbols via alias2SymbolTable function from limma package
# NOTE: use alias2SymbolTable instead of alias2Symbol to avoid errors
# alias2SymbolTable returns NA for gene symbols that cannot be matched
# NOTE: alias2Table issues a warning for gene symbols that cannot be matched
nosology <- read_rds(here("data/Nosology_2023.rds")) |>
  select(
    NOS_ID,
    NOS_Name,
    NOS_MOI,
    NOS_Gene,
    NOS_OMIM,
    NOS_Notes,
    NOS_Group_Name
  ) |>
  mutate(
    HGNC_symbol = alias2SymbolTable(NOS_Gene, species = "Hs"),
    .after = "NOS_Gene"
  ) |>
  add_column(
    AZ_Comments = "", # for comments from AZ
    PMID = "", # for PMID of key articles
    .after = "NOS_OMIM"
  ) 

# find rows in nosology where NOS_Gene is empty
nosology_no_gene <- nosology %>%
  filter(NOS_Gene == "")
nosology_with_gene <- nosology %>%
  filter(NOS_Gene != "")
# there are 34 entries in the nosology with no gene
# reviewed nosology_no_gene and found that they really have no gene
# so no need to fix these


# find rows in nosology where NOS_Gene does not look like a gene symbol
# HGNC gene symbols only consist of uppercase letters, numbers and hyphens
gene_symbol_pattern <- "^[A-Z0-9-]+$"
nosology_complex_gene <- nosology_with_gene %>%
  filter(!grepl(gene_symbol_pattern, NOS_Gene))

# fix entries in nosology_complex_gene using DataEditR
# then comment out not to overwrite the fixed data!
# nosology_complex_gene_fixed <- data_edit(nosology_complex_gene)

# DataEditR tries to guess column_type, which can create problems down the track
# here, DataEditR guessed that PMID was integer, not character
# so converting PMID column from integer to character
# nosology_complex_gene_fixed <- nosology_complex_gene_fixed |>
#  mutate( PMID = as.character(PMID))
#  saveRDS(nosology_complex_gene_fixed, file = here("data/nosology_complex_gene_fixed.rds"))

# document the changes I made to nosology_complex_gene_fixed
compareDF_diff_complex <- compare_df(nosology_complex_gene_fixed, nosology_complex_gene)
create_output_table(
  compareDF_diff_complex,
  output = "xlsx", 
  file_name = here("data/compareDF_diff_complex.xlsx"))
forgts(here("data/compareDF_diff_complex.xlsx"))

# and write the changes into the main nosology object
nosology <- nosology |>
  rows_update(nosology_complex_gene_fixed)

# Nosology does not list all Fanconi genes
# fixing that below
# Fanconi data from OMIM https://omim.org/phenotypicSeries/PS227650
# copy-pasted with with datapasta package
# manually labelled columns of interest (NOS_Name, NOS_MOI, NOS_OMIM, NOS_Gene)
Fanconi <- tribble(
         ~V1,                                        ~NOS_Name,   ~NOS_MOI, ~V4,     ~NOS_OMIM,      ~NOS_Gene,     ~V7,
   "1p36.22", "?Fanconi anemia, complementation group V",  "AR",  3L, 617243L, "MAD2L2", 604094L,
    "1q32.1",  "Fanconi anemia, complementation group T",  "AR",  3L, 616435L,  "UBE2T", 610538L,
    "2p16.1",  "Fanconi anemia, complementation group L",  "AR",  3L, 614083L,   "PHF9", 608111L,
    "3p25.3", "Fanconi anemia, complementation group D2",  "AR",  3L, 227646L, "FANCD2", 613984L,
   "6p21.31",  "Fanconi anemia, complementation group E",  "AR",  3L, 600901L,  "FANCE", 613976L,
    "7q36.1", "?Fanconi anemia, complementation group U",  "AR",  3L, 617247L,  "XRCC2", 600375L,
    "9p13.3",  "Fanconi anemia, complementation group G",  "AR",  3L, 614082L,  "XRCC9", 602956L,
   "9q22.32",  "Fanconi anemia, complementation group C",  "AR",  3L, 227645L,  "FANCC", 613899L,
   "11p14.3",  "Fanconi anemia, complementation group F",  "AR",  3L, 603467L,  "FANCF", 613897L,
   "13q13.1", "Fanconi anemia, complementation group D1",  "AR",  3L, 605724L,  "BRCA2", 600185L,
   "15q15.1",  "Fanconi anemia, complementation group R",  "AD",  3L, 617244L,  "RAD51", 179617L,
   "15q26.1",  "Fanconi anemia, complementation group I",  "AR",  3L, 609053L,  "FANCI", 611360L,
   "16p13.3",  "Fanconi anemia, complementation group P",  "AR",  3L, 613951L,   "SLX4", 613278L,
  "16p13.12",  "Fanconi anemia, complementation group Q",  "AR",  3L, 615272L,  "ERCC4", 133520L,
   "16p12.2",  "Fanconi anemia, complementation group N",    NA,  3L, 610832L,  "PALB2", 610355L,
   "16q23.1", "?Fanconi anemia, complementation group W",  "AR",  3L, 617784L,  "RFWD3", 614151L,
   "16q24.3",  "Fanconi anemia, complementation group A",  "AR",  3L, 227650L,  "FANCA", 607139L,
  "17q21.31",  "Fanconi anemia, complementation group S",  "AR",  3L, 617883L,  "BRCA1", 113705L,
     "17q22",  "Fanconi anemia, complementation group O",  "AR",  3L, 613390L, "RAD51C", 602774L,
   "17q23.2",  "Fanconi anemia, complementation group J",    NA,  3L, 609054L,  "BRIP1", 605882L,
    "Xp22.2",  "Fanconi anemia, complementation group B", "XLR",  3L, 300514L,  "FANCB", 300515L
  ) |>
  select(starts_with("NOS")) |>
  mutate(NOS_ID = "NOS 38-0420", .before = NOS_Name) |>
  mutate(NOS_Group_Name = "Limb hypoplasia—reduction defects group") |>
  mutate(NOS_OMIM = as.character(NOS_OMIM)) |>
  mutate(
    HGNC_symbol = alias2SymbolTable(NOS_Gene, species = "Hs"),
    .after = "NOS_Gene") |>
  add_column(
    AZ_Comments = "", # for comments from AZ
    PMID = "", # for PMID of key articles
    .after = "NOS_OMIM"
  ) |>
  # manually update some data missing from OMIM
  rows_update(tibble (NOS_Gene = c("PALB2","BRIP1"),
                      NOS_MOI = c("AR","AR"),
                      AZ_Comments = c("added AR", "added AR"),
                      PMID = c("17200671","16116423")
                      )
              )

# no need to document the changes as the code above is very clear
# just saving the Fanconi object, just in case
# saveRDS(Fanconi, here("data/Fanconi_fixed.rds"))

# add Fanconi data to nosology
nosology <- nosology |>
  rows_delete(Fanconi) |> # remove old Fanconi data
  rows_insert(Fanconi)


# find rows in nosology where NOS_Gene looks like a gene symbol
nosology_clean_gene <- nosology %>%
  filter(grepl(gene_symbol_pattern, NOS_Gene))

# check that length of nosology matches the sum of the above subsets
# glue("Length of nosology: {nrow(nosology)}")
# glue("Sum of subsets: {nrow(nosology_no_gene) + nrow(nosology_complex_gene) + nrow(nosology_clean_gene)}")

# find rows in nosology where NOS_GENE and HGNC_symbol do not match
# i.e. the NOS_Gene was likely not an HGNC symbol
nosology_HGNC_mismatch <- nosology_clean_gene %>%
  filter(NOS_Gene != HGNC_symbol)
# there were 21 genes where NOS_Gene used non-HGNC symbols

# find rows in nosology where HGNCSymbol is NA
# this includes typos in NOS_Gene, missing gene in NOS_Gene and
# complex entries with multiple genes in NOS_Gene
# nosology_NA_HGNC <- nosology %>%
#  filter(is.na(HGNC_symbol))

# find rows with missing HGNC symbols where NOS_Gene looks like a gene symbol
# and thus Alias2SymbolTable should have been able to find a match
nosology_NA_HGNC_diff <- nosology_clean_gene |>
  filter(is.na(HGNC_symbol))

# fix these errors using DataEditR
# now commented out to avoid overriding the fixed data
# nosology_NA_HGNC_diff_fixed <- data_edit(nosology_NA_HGNC_diff)
# saveRDS(nosology_NA_HGNC_diff_fixed, file = here("data/nosology_NA_HGNC_diff_fixed.rds"))
# there were 7 changes: 5 typos, the HOXD cluster and MAENLI lncRNA

# show differences between original and fixed dataframe
# using waldo
nosology_NA_HGNC_diff <- as.data.frame(nosology_NA_HGNC_diff)
waldo_NA_HGNC_diff <- compare(
  nosology_NA_HGNC_diff,
  nosology_NA_HGNC_diff_fixed
  )
saveRDS(waldo_NA_HGNC_diff, file = here("data/waldo_NA_HGNC_diff.rds"))

# using compareDF
compareDF_NA_HGNC_diff <- compare_df(
  nosology_NA_HGNC_diff_fixed,
  nosology_NA_HGNC_diff,
  keep_unchanged_cols = FALSE
)
compareDF_NA_HGNC_diff_html <- create_output_table(compareDF_NA_HGNC_diff)
saveRDS(compareDF_NA_HGNC_diff, file = here("data/compareDF_NA_HGNC_diff.rds"))
# try exporting to xls and then reimporting with forgts
# to get a gt table of the differences
library(forgts)
create_output_table(
  compareDF_NA_HGNC_diff,
  output = "xlsx", 
  file_name = here("data/compareDF_NA_HGNC_diff.xlsx"))
forgts(here("data/compareDF_NA_HGNC_diff.xlsx"))
# this looks great!


# using diffobj
# not working, not sure why
# but other diff methods work
diffobj_NA_HGNC_diff <- diffPrint(
  nosology_NA_HGNC_diff,
  nosology_NA_HGNC_diff_fixed
)
saveRDS(diffobj_NA_HGNC_diff, file = here("data/diffobj_NA_HGNC_diff.rds"))


# update main nosology with these changes
nosology <- nosology %>%
  rows_update(nosology_NA_HGNC_diff_fixed)

# save the updated nosology object
write_rds(nosology, here("data/nosology.rds"))

# now that gene column is completely updated
# we can run some stats on the data
# ...
