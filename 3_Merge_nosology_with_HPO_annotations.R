#Trying to merge the Nosology with the HPO annotations file
#to explore the state of HPO annotations for skeletal dysplasias
#eg. see how many skeletal dysplasias have no (or very few) HPO annotations
#needs more work?

library(tidyverse)
library(here)

#read the cleaned-up nosology, created in '1 Convert AJMG Table to Nosology object'
nosology <- read_rds(here("data/nosology_cleaned.rds"))

#read HPO annotation file (originally downloaded from https://hpo.jax.org/app/data/annotations)
hpoa <- read_tsv(here("raw_data/HPOA_2023-10-09.tsv"), skip = 4, col_names = TRUE)

#remove ORPHA and DECIPHER entries
hpoa_OMIM <- filter(hpoa, str_detect(database_id, "OMIM.*"))
  
#remove 'OMIM' prefix
hpoa_OMIM <- hpoa_OMIM %>%
  mutate(database_id = str_replace(database_id,"OMIM:",""))

#keep only HPO annotations for diseases that have an OMIM match in the nosology
skeldys_hpoa <- filter(hpoa_OMIM, database_id %in% nosology$NOS_OMIM, )

hpo_count <- skeldys_hpoa %>%
  add_count(database_id, sort = TRUE)


 
