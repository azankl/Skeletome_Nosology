# ## Convert AJMG Table to Nosology Object
# 
# This document deals with converting the Nosology table published in AJMG
# into an R tibble (dataframe).
# 
# ASF provided the "Nosology 2023 table Revised Jan 31.xls' document
# stored in this project, but as he was not sure if this was the final version,
# I decided to scrape the table straight from the published version using rvest.
# 
# The 2023 Nosology is published in the AJMG here:
# https://onlinelibrary.wiley.com/doi/10.1002/ajmg.a.63132
# 
# I could not read from this URL directly using `read_html`
# probably because of some paywall issue.
# I therefore opened the URL in Chrome, and saved the page as "AJMG.html"
# using Chrome \> File \> Save Page As... with Format "Web Page, HTML Only".
# 
# I then opened `AJMG.html` with Chrome and used the
# rvest SelectorGadget (https://rvest.tidyverse.org/articles/selectorgadget.html)
# to identify the selector for the table
# (which turned out to be called `".article-section__table"`).
# 
# I then extracted the table from AJMG.html with the code below.
# 
# Note that some quotation marks in the Notes column got garbled. Probably not
# worth fixing as they are just Notes, but need to check nothing else got garbled.
# 
# Selecting the xpath as described [here]
# (https://www.r-bloggers.com/2015/01/using-rvest-to-scrape-an-html-table/)
# also works, but some of the code in the blog post does not work,
# must be an older version of rvest.

library(tidyverse)
library(here)
library(rvest)

######################################################
# Step 1: read the AJMG table from the saved html file

# read the AJMG table from the saved html file.
# see comments above for finding the correct html_element
# and comment below for finding the correct encoding
AJMG_table<- read_html (here("raw_data/AJMG.html"), encoding = "UTF-8") %>%
  html_element(".article-section__table") %>%
  html_table()

# Below helped to find the correct encoding.
# Using the default read_html (no encoding specified) created some weird characters.
html_encoding_guess(AJMG_table)

#######################################################
# Step 2: add Nosology Group-Number and Group-Name as extra columns to each disease
# and then remove the Group headings, NOS-numbers and comments

#this regex pattern matches rows that start with "Group", followed by whitespace and a number, and captures this number
pattern <- "Group\\s([1-9]+)"
#extract rows where Group Number/Disorder Number column matches this pattern
grouphits<- str_match(unlist(AJMG_table[,1]), pattern)
colnames(grouphits)<-c("Match","GroupNumber")
#add new column to AJMG_table with group number                      
AJMG_table <- add_column(AJMG_table,NOS_Group_No = grouphits[,'GroupNumber'])
#add new column for NOS_Group_Name (empty for now)
AJMG_table <- add_column(AJMG_table, NOS_Group_Name = NA)
#rows that did not match the group pattern above (ie. are disease rows) will have NA in NOS_Group_No
#fill those rows with the group number extracted previously and
#fill NOS_Group_Name with Names from 'Name of Group/Disorder'
for (i in seq_along(unlist(AJMG_table[,"NOS_Group_No"]))) {
  groupNumber = ifelse(!is.na(AJMG_table[i,"NOS_Group_No"]), AJMG_table[i,"NOS_Group_No"], groupNumber)
  groupName = ifelse(!is.na(AJMG_table[i,"NOS_Group_No"]), AJMG_table[i,2], groupName)
  AJMG_table[i,c("NOS_Group_No","NOS_Group_Name")] = c(groupNumber,groupName)
}
#now some clean up:
#remove rows without NOS-number (ie Group Headings and Group Comments)
AJMG_table <- AJMG_table %>%
  filter (str_detect(`Group number/number of disorder`,"NOS"))
#decided to keep the NOS-numbers column for now, as it allows easy sorting according to Nosology Groups
#rearrange an rename columns
Nosology_2023 <- AJMG_table %>% 
  rename (
    "NOS_ID" = `Group number/number of disorder`,
    "NOS_Name" = `Name of group/name of disorder`,
    "NOS_MOI" = `Inheritance`,
    "NOS_Gene" = `Gene or locus`,
    "NOS_OMIM" = `MIM No.`,
    "NOS_Notes" = `Notes`
  )

#######################################################
# Step 3: save the Nosology_2023 object to file

# Option 1: as a serialised RDS object. To read back in, use read_rds()
write_rds (Nosology_2023, here("data/Nosology_2023.rds"))

# Option 2: as a csv or tsv file.
## csv is easier to view in RStudio, as the text includes commas, write_csv puts quotation marks around those fields, which might add complexity in downstream applications.
# tsv does not need to work around commas in the text, but cannot be imported directly into RStudio Viewer.
# both csv and tsv versions display correctly in Easy CSV Editor app.

#write_tsv(Nosology_2023, here("data/Nosology_2023.tsv"))
#write_csv(Nosology_2023, here("data/Nosology_2023.csv"))

######################################################
# Step 4: clean up the Nosology_2023 object

# The ISDS Nosology sometimes has multiple DMIM values in the DMIM field.
# split such entries into individual rows
# the regex splits at comma plus optional white space
nosology <- Nosology_2023 %>%
  separate_rows(NOS_OMIM, sep=",\\s+") %>%
  add_column(Comment = "NA") #add new Comment column to nosology dataframe

write_rds (nosology, here("data/nosology.rds"))

