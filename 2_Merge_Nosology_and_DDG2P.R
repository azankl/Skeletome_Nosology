#Merge Nosology with DDG2P data
#needs updating!

#uses the Nosology object created in 'Convert AJMG Table to Nosology object'
#IMPORTANT: DDG2P_1_12_2021.csv is the last version of DDG2P that used G2Ps older nomenclature that includes functional consequence terms like 'activating'
#So plan 2 step approach:
#First merge Nosology to lastest DDG2P, to match the most diseases and also to have a starting point for the Skeletal G2P project
#Then merge DDG2P_1_12_2021 to this combined table, to get the older nomenclature where it exists


library(tidyverse)
library(here)

#read the cleaned-up nosology, created in '1 Convert AJMG Table to Nosology object'
nosology <- read_rds(here("data/nosology_cleaned.rds"))

#download latest version of DDG2P or use a stored snapshot
#use snapshots for development as new data might be brake things
#but try with latest version from time to time
#ddg2p <- read_csv("https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz")
ddg2p <- read_csv(here("raw_data/DDG2P_12_5_2023.csv"))

#take note of some importing problems, these are actually errors in ddg2p:
problems(ddg2p)

#filter ddg2p for entries that have a matching DMIM number in Nosology
skg2p_semiJoin_byDMIM <- semi_join(ddg2p, nosology, by = join_by ('disease mim' == 'NOS_OMIM'))
#save as the first version of skg2p to publish online:
write_csv(skg2p_semiJoin_byDMIM,here("data/skg2p_v1.csv"), col_names = TRUE)

#creating an anti_join to see which entries in nosology are NOT in DDG2P
nosology_antijoin_byDMIM <- anti_join(nosology, ddg2p, by = join_by ('NOS_OMIM' == 'disease mim'))

#some entries in NOS_OMIM are empty (no OMIM) or cross-references to other OMIM numbers ("see xxxxxx")
#what to do with these?


