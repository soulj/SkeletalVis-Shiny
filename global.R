library(feather)
library(readr)
library(enrichR)

foldChangeTable <- read_feather("foldChangeTable.feather")
pvalTable <- read_feather("pvalTable.feather")

#load the accessions
accessions <- read.delim("accessions.txt",stringsAsFactors = F)
accessions$combined <- paste0(accessions$accession,"_",accessions$comparison)

#Load the experiment table
expTable<-read_csv(file = "data/expTable.csv")

#load the signature lists for the response comparisons
upSigs <- readRDS(file="foldChangeListUp.pval.RDS")
downSigs <- readRDS(file="foldChangeListDown.pval.RDS")

#load the chrDir lists
chrDirsList <- readRDS("chrDirs.RDS")

#load the human2otherspecies homology table
human2otherspecies <- as.data.frame(read_feather("data/human2otherspecies.feather"))

#get the enrichR databases
databases <- sort(listEnrichrDbs()[,1])
