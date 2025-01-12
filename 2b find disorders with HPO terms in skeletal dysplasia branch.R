#Find all disorders with HPO terms in the Abnormality of the skeletal system (HP:0000924) branch
#to see if these are all skeletal dysplasias and should be added to the nosology

library(tidyverse)
library(here)
library(ontologyIndex)
data(hpo)


#read the nosology, created in '1a Convert AJMG Table to Nosology object'
nosology <- read_rds(here("data/nosology.rds"))

#read HPO annotation file (originally downloaded from https://hpo.jax.org/app/data/annotations)
hpoa <- read_tsv(here("raw_data/HPOA_2024-12-12.tsv"), skip = 4, col_names = TRUE)

#remove ORPHA and DECIPHER entries
hpoa_OMIM <- filter(hpoa, str_detect(database_id, "OMIM.*"))

#get all descendants of Abnormality of the skeletal system (HP:0000924) 
Abnormality_of_the_skeletal_system <- get_descendants(hpo, "HP:0000924")

#filter hpoa for all rows where the HPO term is in the Abnormality_of_the_skeletal_system
#then group by database_id and count the number of rows
hpoa_skeletal <- hpoa_OMIM %>%
  filter(hpo_id %in% Abnormality_of_the_skeletal_system) %>%
  select(database_id, disease_name) %>%
  add_count(database_id) %>%
  distinct() %>%
  arrange(desc(n))
#this shows that for example Schaaf-Yang syndrome has 14 skeletal HPO terms,
#but this is not really a skeletal dysplasia, not listed in Nosology or PanelApp
#so this is not a good way of defining skeletal dysplasias
#also, hpoa_skeletal contains almost 4000 disorders, which is clearly too many to be useful












