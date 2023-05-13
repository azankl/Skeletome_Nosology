#Merge Nosology with DDG2P data
#uses the Nosology object created in 'Convert AJMG Table to Nosology object'
#IMPORTANT: DDG2P_1_12_2021.csv is the last version of DDG2P that used G2Ps older nomenclature that includes functional consequence terms like 'activating'
#So plan 2 step approach:
#First merge Nosology to lastest DDG2P, to match the most diseases and also to have a starting point for the Skeletal G2P project
#Then merge DDG2P_1_12_2021 to this combined table, to get the older nomenclature where it exists


here::i_am("Merge Nosology to DDG2P.R")

library(tidyverse)
library(here)

nosology <- read_rds("Nosology_2023.rds")
#download latest version of DDG2P or use a stored snapshot
#use snapshots for development as new data might be brake things
#but try with latest version from time to time
#ddg2p <- read_csv("https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz")
ddg2p <- read_csv(here("DDG2P_12_5_2023.csv"))

#take note of some importing problems:
problems(ddg2p)

skg2p_innerJoin <- inner_join(ddg2p, nosology, by = join_by ('disease mim' == 'NOS_OMIM'))
skg2p_semiJoin <- semi_join(ddg2p, nosology, by = join_by ('disease mim' == 'NOS_OMIM'))

#use mutate with replace to correct errors in the table
#as described here: https://stackoverflow.com/questions/36924911/how-to-assign-a-value-to-a-data-frame-filtered-by-dplyr
#need to decide if better to correct errors before of after joining (I think better before) 

#also do join by DMIM and check for rows where genes dont match, that must be an error
