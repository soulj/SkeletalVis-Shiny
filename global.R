library(feather)
library(readr)

foldChangeTable <- read_feather("data/foldChangeTable.RDS")

#load the accessions
accessions <- read.delim("data/accessions.txt",stringsAsFactors = F)
accessions$combined <- paste0(accessions$accession,"_",accessions$comparison)

#Load the experiment table
expTable<-read_csv(file = "data/expTable.csv")

#load the signature lists for the response comparisons
upSigs <- readRDS(file="data/foldChangeListUp.pval.RDS")
downSigs <- readRDS(file="data/foldChangeListDown.pval.RDS")

#load the chrDir lists
chrDirsList <- readRDS("data/chrDirs.RDS")

#load the human2otherspecies homology table
human2otherspecies <- as.data.frame(read_feather("data/human2otherspecies.feather"))
