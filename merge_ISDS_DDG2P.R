#Trying to combine the ISDSnosology.csv with the DDG2P.csv
#using 'join'

library(readr)
library(here)
library (dplyr)
library (tidyr)
library (magrittr)
library (purrr)
library(readxl)
library(stringr)

#Import the ISDS 2019 Table
ISDS_2019_Table <- read_excel("ISDS 2019 Table.xlsx", col_names = c("Name", "MOI","Gene", "OMIM", "ORPHA","Notes"), skip=1)

#extract the Group name and Group Number from the table:
#find rows that start with a number and dot in the Name column, indicating a group heading.
#Then extract group name and number from this row 
pattern <- "^([1-9]+)\\.\\s(.*$)"
grouphits<-unlist(ISDS_2019_Table[,"Name"]) %>%
  str_match(pattern)
colnames(grouphits) <- c("Match","Number","Name")
#rows that did not match the group pattern above (ie. are disease rows) will return NA in grouphits
#fill those rows with the group name and number extracted previously
for (i in seq_along(grouphits[,'Match'])) {
  groupName = ifelse(!is.na(grouphits[i,'Match']), grouphits[i,'Name'], groupName)
  groupNumber = ifelse(!is.na(grouphits[i,'Match']), grouphits[i,'Number'], groupNumber)
  grouphits[i,2:3] = c(groupNumber, groupName)
}
#add new column to ISDS table with group number and group name                      
ISDS <- cbind(ISDS_2019_Table, grouphits[,c('Number','Name')])

#still need to clean up the table a bit, eg remove rows with Group names in it
#this is unfinished


#the parts below seem to come from mere_ISDS_DDG2P_old and make use of "ISDS2019nosology.csv"
#the below wont work
nosology <- nosology_raw %>%
  select(-c(1,2)) %>% #remove unnecessary columns
  drop_na(GroupN) %>% #removes headers in table that contain no useful data
  rename(DMIM = Omim, GENE = Gene) #unify names for join



#read the DDG2P list:
ddg2p<-read_csv(here("DDG2P_21_7_2022.csv"), col_names = TRUE) %>%
  rename(DMIM = "disease mim", GENE = "gene symbol")

#creating an inner_join shows where ISDS and DDG2P match on DMIM
#but also reveals inconsistencies where other columns dont add up
#for example DMIM 184255 is associated with Gene FN1 in DDG2P, but not in ISDS
#can be a useful view to identify these discrepancies
innerjoin_byDMIM <-inner_join(nosology,ddg2p, by = "DMIM")
write_csv(innerjoin_byDMIM,here("innerjoinByDMIM.csv"), col_names = TRUE)

#tried some other joins to explore how ISDS and DDG2P differ:

#creating a semi_join to see how many entries in ISDS have a match in DDG2P (based on DMIM)
#more on joins here: https://r4ds.had.co.nz/relational-data.html
semijoin_byDMIM <-semi_join(nosology,ddg2p, by = "DMIM")
write_csv(semijoin_byDMIM,here("semijoinByDMIM.csv"), col_names = TRUE)

#creating an anti_join to see how many entries in ISDS are NOT in DDG2P 
#more on joins here: https://r4ds.had.co.nz/relational-data.html
antijoin_byDMIM <-anti_join(nosology,ddg2p, by = "DMIM")
write_csv(antijoin_byDMIM,here("antijoinByDMIM.csv"), col_names = TRUE)

#see which ISDS entries are not in DDG2P, but their gene is in DDG2P (suggests mapping errors based on DMIM)
noDMIM_butGene<- inner_join(antijoin_byDMIM, ddg2p, by = "GENE")

#find out duplicates created by  innerjoin (the reason why innerjoin has more entries than semijoin)
dups_innerjoin <- innerjoin_byDMIM %>%
  filter(duplicated(DMIM))

#find number of named disorders in ISDS (ie dont count subentries without a name)
nosology_named <- nosology %>%
  filter(!is.na(Group))
#interestingly, the result is 463, while the paper says 461

#find number of unique genes in ISDS
n_distinct(nosology$GENE)
#the result is 438 while the paper says 437

#find number of Orphanet only entities in ISDS
nosology_orphaonly <- nosology %>%
  filter ((is.na(DMIM) & !is.na(Orph)))

#find disorders without a gene in ISDS (has a name but no gene)
nosology_nogene <- nosology %>%
  filter((is.na(GENE) & !is.na(Group)))

#creating a left join, which keeps all ISDS entries
leftjoinedByDMIM <- left_join(nosology, ddg2p, by = "DMIM")
#differences found
#48 FMD: wrong DMIM in DDG2P
#60 Larsen:
write_csv(leftjoinedByDMIM,here("leftjoinedByDMIM.csv"), col_names = TRUE)



#testing on how to update the table:
#
#this does not work because DMIM values are not unique (some DMIM values occur more than once)
#leftjoinedByDMIM_updated <- leftjoinedByDMIM %>%
#  rows_update(tibble(DMIM = 616583, Group = "Stanescu"))
#
#this works, but need to replace the NA values first (because we access rows by matching DMIM, and we dont know if NA is a match, so its throws an error)
 leftjoinedByDMIM_updated <- leftjoinedByDMIM %>%
   mutate(across(everything(),as.character)) %>%
    mutate(across(everything(), replace_na, "0"))
 leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="616583",c("Group","GENE.x", "inheri") ] <- tibble("Stanescu","COL2A1","AD") 
#nice feature in RStudio: type the column names in the example above without quotation marks, RStudio will offer auto-completion of column names and add the required quotation marks!

#replacing values with mutate and replace should also work, but is more verbose
#needs mutate/replace statements for each variable to be updated, cannot group variables with c() and tibble () as above
#more on mutate/replace here: https://community.rstudio.com/t/mutate-and-replace-question/55235/4
 
 
 
#adding rows is easier than updating existing rows, see add.row() or tribble()
#or create the new rows in a spreadsheet and then import them and merge them to the main table in R
#this way it remains traceable what was added

 
#find entries that have the same DMIM but not the same gene
joined_diffGene <- leftjoinedByDMIM %>%
  filter(GENE.x != GENE.y) %>%
  select(-c(5,6)) #drop Orpha and GroupN

#fixing errors seen in joined_diffGene Table

#error in DDG2P data, emailed to David Fitzpatrick:
#The disease MIM for MAP3K7-associated Frontometaphyseal Dysplasia should be 617137, not 305620 (the latter is for FLNA-associated FMD).

#wrong gene name in ISDS
leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="614078","GENE.x"] <- leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="614078","GENE.y"]
#wrong DMIM and name in G2P (this fixes it in the joined table, but better to correct in DDG2P and then reimport)
leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="617662",c(7:22)] <- leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$GENE.y=="GZF1",c(7:22)]
leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="617662", "disease name"] <- leftjoinedByDMIM_updated[leftjoinedByDMIM_updated$DMIM=="617662","Group"]
leftjoinedByDMIM_updated <- leftjoinedByDMIM_updated %>%
  filter(!(DMIM=="150250" & GENE.y=="GZF1"))



#trying to create a clean ISDS list from original ISDS table
ISDS2019nosology <- read_csv("ISDS2019nosology.csv")


#Use RNotebook to document the merging process, will help with future publication

#look at:
#Nosology entries with no gene
#Nosology entries with no DMIM
#try inner join on DMIM to get an idea how much already is covered in G2P

#Goal
#1: get a list of all skeletal dysplasias already in DDG2P
#to create the first Skeletal G2P list.
#
#2: get a list of all skeletal dysplasias that should be included in Notion
#with DDG2P data filled in where available
#this list is bigger as it includes:
#disorders with a gene but not yet included in DDG2P
#disorders without a gene

#Thoughts on how to solve the many inconsistencies:
#Start with disorders in ISDS that have a gene
#use MONDO as the identifier (not DMIM)\
#create new MONDO-IDs for ISDS disorders with a gene, if does not already exist
#...more to come here