#Add HGNC gene symbols to the nosology 
#in the process, we also explore some issues with the NOS_Gene column

library(tidyverse)
library(here)
library(glue)

# needed for alias2SymbolTable function, 
# also needs up-to-date org.Hs.eg.db package installed
library(limma) 


#read the original nosology dataframe (uncleaned)
nosology <- read_rds(here("data/Nosology_2023.rds")) |>
  select(NOS_ID, NOS_Name, NOS_MOI, NOS_Gene, NOS_OMIM, NOS_Notes, NOS_Group_Name) |>
  add_column(HGNC_symbol = "", .after = "NOS_Gene") #add a HGNC column

#find rows in nosology where NOS_Gene is empty
nosology_NoGene <- nosology %>%
  filter(NOS_Gene == "")
nosology_with_gene <- nosology %>%
  filter(NOS_Gene != "")

#find rows in nosology where NOS_Gene does not look like a gene symbol
#(ie. does not consist of only uppercase letters and numbers)
gene_symbol_pattern <- "^[A-Z0-9]+$"
nosology_complex_Gene <- nosology_with_gene %>%
  filter(!grepl(gene_symbol_pattern, NOS_Gene))

#find rows in nosology where NOS_Gene looks like a gene symbol
nosology_cleanGene <- nosology %>%
  filter(grepl(gene_symbol_pattern, NOS_Gene))
  
#check that length of nosology matches the sum of the above subsets
#print out the calculation and the result
glue("Length of nosology: {nrow(nosology)}")
glue("Sum of subsets: {nrow(nosology_NoGene) + nrow(nosology_complex_Gene) + nrow(nosology_cleanGene)}")


#add HGNC gene symbols via alias2SymbolTable function from limma package
#NOTE: use alias2SymbolTable instead of alias2Symbol to avoid errors
#alias2SymbolTable returns NA for gene symbols that cannot be matched
#NOTE: alias2Table issues a warning for gene symbols that cannot be matched
nosology <- nosology |>
  mutate(HGNC_symbol = alias2SymbolTable(NOS_Gene, species = "Hs"))

#find rows in nosology where HGNCSymbol is NA
nosology_NA_HGNC <- nosology %>%
  filter(is.na(HGNC_symbol))

# find rows with missing HGNC symbols where NOS_Gene looks like a gene symbol
# (and thus Alias2SymbolTable should have been able to find a match)
nosology_NA_HGNC_diff <- nosology_NA_HGNC |>
  filter(grepl(gene_symbol_pattern, NOS_Gene))


# find rows in nosology where NOS_GENE and HGNC_symbol do not match
# i.e. the NOS_Gene was likely not an HGNC symbol
nosology_HGNC_mismatch <- nosology %>%
  filter(NOS_Gene != HGNC_symbol)

#save the updated nosology dataframe
saveRDS(nosology, here("data/Nosology_2023_HGNC.rds"))
