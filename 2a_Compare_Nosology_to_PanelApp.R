#Compare Nosology to PanelApp data


library(tidyverse)
library(here)

#read the cleaned-up nosology, created in '1 Convert AJMG Table to Nosology object'
nosology <- read_rds(here("data/nosology_cleaned.rds"))

#read PanelAppAU Skeletal Dysplasia data
PanAppAU <- read_tsv(here("raw_data/Skeletal dysplasia_PanelApp_AU_Version_0.298.tsv"))

#select gene column from nosology and PanelApp
nosology_genes <- nosology %>% select(NOS_Gene)
PanAppAU_genes <- PanAppAU %>% select(`Gene Symbol`)

#find the genes that are in both datasets
common_genes <- intersect(nosology_genes$NOS_Gene, PanAppAU_genes$`Gene Symbol`)

#find the genes that are in PanelApp but not in nosology
PanAppAU_only_genes <- setdiff(PanAppAU_genes$`Gene Symbol`, nosology_genes$NOS_Gene)

#find the genes that are in nosology but not in PanelApp
nosology_only_genes <- setdiff(nosology_genes$NOS_Gene, PanAppAU_genes$`Gene Symbol`)

#the above only gives lists of genes, which is not very informative
#need disease name, OMIM number etc., so will use anitjoin below

#antijoin nosology with PanelApp on gene symbol
nosology_antijoin_byGene <- anti_join(nosology, PanAppAU, by = join_by ('NOS_Gene' == 'Gene Symbol'))
PanAppAU_antijoin_byGene <- anti_join(PanAppAU, nosology, by = join_by ('Gene Symbol'== 'NOS_Gene'))

#rearrange columns in PanAppAU_antijoin_byGene
PanAppAU_antijoin_byGene <- PanAppAU_antijoin_byGene %>% 
  select(`Gene Symbol`, `Phenotypes`, `Omim` ) |>
  mutate(OMIM = str_extract_all(Phenotypes, "\\d{6}")) 
  #unnest(c('Gene Symbol', 'Phenotypes', 'OMIM'))

  
  

# #read PanelAppUK Skeletal Dysplasia data
# PanAppUK <- read_tsv(here("raw_data/Skeletal dysplasia_PanelApp_UK_Version_7.5.tsv"))
# 
# #antijoin PanelAppAU with PanelAppUK on gene symbol
# PanAppAUvsUK_antijoin_byGene <- anti_join(PanAppAU, PanAppUK, by = join_by ('Gene Symbol' == 'Gene Symbol')) |>
#   select(`Gene Symbol`, `Phenotypes`, `Omim` )


