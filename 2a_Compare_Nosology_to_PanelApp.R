#Compare Nosology to PanelApp data


library(tidyverse)
library(here)

#read the cleaned-up nosology, created in '1 Convert AJMG Table to Nosology object'
nosology <- read_rds(here("data/nosology_cleaned.rds"))

#read PanelApp Skeletal Dysplasia data
PanApp <- read_tsv(here("raw_data/Skeletal dysplasia_PanelApp_AU_Version_0.298.tsv"))

#select gene column from nosology and PanelApp
nosology_genes <- nosology %>% select(NOS_Gene)
PanApp_genes <- PanApp %>% select(`Gene Symbol`)

#find the genes that are in both datasets
common_genes <- intersect(nosology_genes$NOS_Gene, PanApp_genes$`Gene Symbol`)

#find the genes that are in PanelApp but not in nosology
PanApp_only_genes <- setdiff(PanApp_genes$`Gene Symbol`, nosology_genes$NOS_Gene)

#find the genes that are in nosology but not in PanelApp
nosology_only_genes <- setdiff(nosology_genes$NOS_Gene, PanApp_genes$`Gene Symbol`)

nosology_only_genes
