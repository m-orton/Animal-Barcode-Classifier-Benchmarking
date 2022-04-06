# install.packages("foreach")
library(foreach)
# For sequence alignments we need the biostrings (DNAStringSet function) 
# and muscle libraries, as follows:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings) 
# BiocManager::install("DECIPHER")
library(DECIPHER)
# install.packages("dplyr")
library(dplyr)
# install.packages("ggplot2")
library(ggplot2)
# For rBLAST
# install.packages("remotes")
# library(remotes)
# remotes::install_github("mhahsler/rBLAST")
library(rBLAST)
# For plotting
# install.packages("viridis")
library(viridis)
# BiocManager::install("rRDP")
library(rRDP)
# install.packages("data.table")
library(data.table)
# install.packages("ROCR")
library(ROCR)
# install.packages("ape")
library(ape)
# install.packages("caret")
library(caret)
# install.packages("ggnewscale")
library(ggnewscale)
# install.packages("RColorBrewer")
library(RColorBrewer)
# install.packages("cowplot")
library(cowplot)
# install.packages("ggpubr")
library(ggpubr)
# install.packages(c("taxize", "usethis"))
library(taxize)
library(usethis)
# install.packages("rfishbase")
library(rfishbase)
# For getting fasta files I need in FASTA
# writeXStringSet(dna, ".fas", format="FASTA")

########
# Groups this pipeline is testing currently:

# Amphibibia - Family, genus and species level benchmarking 

# Mammalia - Family, genus and species level benchmarking 

# Aves - Family, genus and species level benchmarking 

# Gastropoda - Family, genus and species level benchmarking 

# Araneae - Family, genus and species level benchmarking 

# Hymenoptera - Family, genus and species level benchmarking 

#########
# Dataset download

# Diptera
dnaNames <- strsplit(names(dna), "_")
genus <- sapply(dnaNames, `[`, 1)
species <- sapply(dnaNames, `[`, 2)
species <- paste(genus, "_", species, sep="")

# diptera <- data.frame(species=dfInitial$species_name, genus=dfInitial$genus_name, family=dfInitial$family_name)
# write.csv(diptera, "dipteraTaxonomy.csv")

dipTax$species <-  gsub(" ", "_", dipTax$species)
dfLep$species <- dfLep$sp
dfLep2 <- merge(dfLep, dipTax, by=c("species"), all.x = TRUE, all.y = FALSE)
naFam <- which(is.na(dfLep2$family))
dfLep2 <- dfLep2[-naFam,]
dfLep2 <- dfLep2[!duplicated(dfLep2[c("id")]),]
dfLep2$id <- 1:nrow(dfLep2)
dfLep <- dfLep2
dfLep$names <- dfLep$family

# Fish
dna <- readDNAStringSet('fishDNAUnfiltered.fas', format="FASTA")

dnaNames <- strsplit(names(dna), "_")
genus <- sapply(dnaNames, `[`, 3)
species <- sapply(dnaNames, `[`, 4)
species <- unique(species)
fish <- species(species)
family <- data.frame(species=fish$Species, genus=fish$Genus, family=fish$FamCode)
family$species <-  gsub(" ", "_", family$species)
species <- sapply(dnaNames, `[`, 4)
species <- gsub(" ", "_", species)

dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$genus <- genus
dfLep$species <- species
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

dfLep2 <- merge(dfLep, family, by=c("species"), all.x = TRUE, all.y = FALSE)
naFam <- which(is.na(dfLep2$family))
dfLep2 <- dfLep2[-naFam,]
dfLep2$id <- 1:nrow(dfLep2)
dfLep2$names <- dfLep2$family
dfLep2$genus <- dfLep2$genus.x
dfLep <- dfLep2

# Amphibia
dna <- readDNAStringSet("amphibiaDNAUnfiltered.fas", format="FASTA")

# Mammalia
dna <- readDNAStringSet("mammaliaDNAUnfiltered.fas", format="FASTA")

# Gastropoda
dna <- readDNAStringSet("gastropodaDNAUnfiltered.fas", format="FASTA")

# Aves
eBird <- read_csv("https://api.ebird.org/v2/ref/taxonomy/ebird?fmt=csv")
eBird$SCIENTIFIC_NAME <- gsub(" ", "_", eBird$SCIENTIFIC_NAME)
family_genus_ebird <- data.frame(family=eBird$FAMILY_SCI_NAME)
ebirdgenus <- strsplit(family_species$species, "_")
family_genus_ebird$genus <- sapply(ebirdgenus, `[`, 1)
avesTaxonomy <- read_tsv("aves_data.tsv")
family_genus <- data.frame(family=avesTaxonomy$family_name, genus=avesTaxonomy$genus_name)

dna <- readDNAStringSet("Birds_COI.fas", format="FASTA")

dnaNames <- strsplit(names(dna), "_")
genus <- sapply(dnaNames, `[`, 3)
species <- sapply(dnaNames, `[`, 4)
species <- gsub(" ", "_", species)

dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$genus <- genus
dfLep$species <- species

# Remove duplicate sequences
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

family <- which(dfLep$genus %in% family_genus$genus)
family <- family_genus[family,]
family <- unique(family)

dfLep2 <- merge(dfLep, family, by=c("genus"), all.x = TRUE, all.y = FALSE)
naFam <- which(is.na(dfLep2$family))
dfnaFam <- dfLep2[naFam,]
dfnaFam <- dfnaFam[,1:4]
dfLep2 <- dfLep2[-naFam,]

dfnaFam <- merge(dfnaFam, family_genus_ebird, by=c("genus"), all.x = TRUE, all.y = FALSE)
dfnaFam <- unique(dfnaFam)
naFam <- which(is.na(dfnaFam$family))
dfnaFam <- dfnaFam[-naFam,]

dfLep2 <- rbind(dfLep2, dfnaFam)
dfLep <- dfLep2
dfLep$id <- 1:nrow(dfLep)
dfLep$names <- dfLep$family

# Araneae
# Renaming dataset with full taxonomy
names(dna) <- paste(dfInitial$family_name, dfInitial$genus_name, dfInitial$species_name, sep=";")

dnaNames <- strsplit(names(dna), ";")

# For family level
family <- sapply(dnaNames, `[`, 1)
genus <- sapply(dnaNames, `[`, 2)
species <- sapply(dnaNames, `[`, 3)

# Start of stratified sampling
dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$names <- family
dfLep$genus <- genus
dfLep$species <- species

# Remove duplicate sequences
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

# Gastropoda
dna <- readDNAStringSet("gastropodaDNAUnfiltered.fas", format="FASTA")

# Hymenoptera
writeXStringSet(dna, "hymenopteraDNAUnfiltered.fas", format="FASTA")
dna <- readDNAStringSet("hymenopteraDNAUnfiltered.fas", format="FASTA")

# Lepidoptera
writeXStringSet(dna, "lepDNAUnfiltered.fas", format="FASTA")
dna <- readDNAStringSet("lepDNAUnfiltered.fas", format="FASTA")

# Birds
writeXStringSet(dna, "lepDNAUnfiltered.fas", format="FASTA")
dna <- readDNAStringSet("lepDNAUnfiltered.fas", format="FASTA")

#########
# Start of pipeline - Same process for all groups (repeated for each group until ROC and MCC plots)

#########
# For porter datasets
dnaNames <- names(dna)

dnaNames <- strsplit(dnaNames, " ")

# For family level
familyID <- sapply(dnaNames, `[`, 2)
familyID <- strsplit(familyID, ";")
familyID <- sapply(familyID, `[`, 7)

# Start of stratified sampling
dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$names <- familyID

# Remove duplicate sequences
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

# Picked 5 largest families
# largestFam <- tail(sort(table(dfLep$names)), 5)
# dfLepSub <- foreach(i=1:nrow(largestFam)) %do% which(dfLep$names == names(largestFam)[i])
# dfLepSub <- do.call(c, dfLepSub)
# dfLepSub <- dfLep[dfLepSub,]

#########

# Total num of sequences used
nrow(dfLep)

# Number of families
length(table(unique(dfLep$names)))

# Break down by family identifier (or other taxonomic rank for full sequence datasets)
nameList <- lapply(unique(dfLep$names), function(x) 
  dfLep[dfLep$names == x,])

# Conversion to DNAStringSet and naming with family and sequence id
nameList2 <- sapply( nameList , function(x) DNAStringSet( x$dna ) )

for (i in seq(from = 1, to = length(nameList), by = 1)) {
  names(nameList2[[i]]) <- paste(nameList[[i]]$names, nameList[[i]]$id, sep=";")
}

# Removal of gaps in sequences before testing
nameList2 <- foreach(i=1:length(nameList2)) %do% RemoveGaps(nameList2[[i]],
                                                            removeGaps = "all",
                                                            processors = 1)

# Create 3 folds for cross-validation testing with the package caret
allFolds <- foreach(i=1:length(nameList2)) %do% createFolds(nameList2[[i]], k=3, list = TRUE)

# Fold 1
fold1 <- foreach(i=1:length(nameList2)) %do% nameList2[[i]][allFolds[[i]]$Fold1]
fold1 <- do.call(c, fold1)
# Fold 2
fold2 <- foreach(i=1:length(nameList2)) %do% nameList2[[i]][allFolds[[i]]$Fold2]
fold2 <- do.call(c, fold2)
# Fold 3
fold3 <- foreach(i=1:length(nameList2)) %do% nameList2[[i]][allFolds[[i]]$Fold3]
fold3 <- do.call(c, fold3)

#######
# IDTAXA Testing (Murali at al., 2019)

# Save true taxon and ids for testing sequences
fam_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

fam_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

fam_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

###########
# For lep to cut down on larger families/genus/species
df_f1 <- data.frame(nucleotides=fold1, name=fam_f1, id=famids_f1)
df_f1 <- data.table(df_f1, key="name")
df_f1 <- df_f1[, head(.SD, 350), by=name]

df_f2 <- data.frame(nucleotides=fold2, name=fam_f2, id=famids_f2)
df_f2 <- data.table(df_f2, key="name")
df_f2 <- df_f2[, head(.SD, 350), by=name]

df_f3 <- data.frame(nucleotides=fold3, name=fam_f3, id=famids_f3)
df_f3 <- data.table(df_f3, key="name")
df_f3 <- df_f3[, head(.SD, 350), by=name]

fold1 <- DNAStringSet(df_f1$nucleotides)
fold2 <- DNAStringSet(df_f2$nucleotides)
fold3 <- DNAStringSet(df_f3$nucleotides)

names(fold1) <- paste(df_f1$name, df_f1$id, sep=";")
names(fold2) <- paste(df_f2$name, df_f2$id, sep=";")
names(fold3) <- paste(df_f3$name, df_f3$id, sep=";")

# Save true taxon and ids for testing sequences
fam_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

fam_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

fam_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

############
# Folds 1+2 for testing
fold1_2 <- append(fold1, fold2)

# Folds 2+3
fold2_3 <- append(fold2, fold3)

# Folds 1+3
fold1_3 <- append(fold1, fold3)

# Combined for training sequences
fam_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 1)
famids_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 2)

fam_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 1)
famids_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 2)

fam_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 1)
famids_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 2)

# For mammal dataset only
# names(fold1) <- paste(fam_f1, famids_f1, sep=";")
# names(fold2) <- paste(fam_f2, famids_f2, sep=";")
# names(fold3) <- paste(fam_f3, famids_f3, sep=";")
# names(fold1_2) <- paste(fam_f1_2, famids_f1_2, sep=";")
# names(fold2_3) <- paste(fam_f2_3, famids_f2_3, sep=";")
# names(fold1_3) <- paste(fam_f1_3, famids_f1_3, sep=";")

#########
# Family Testing - 1st run - F1_F2/F3
taxonomy <- paste("Root",fam_f1_2, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_2, taxonomy)

#############
# classify the test sequences
ids <- IdTaxa(fold3, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- fam_f3
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confF1_2 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confF1_2$outcome <- confF1_2$trueTaxa == confF1_2$assignedTaxa
mccID1_2 <- prediction(confF1_2$familyConf, confF1_2$outcome)
mccID1_2 <- performance(mccID1_2,"mat")
mccID1_2 <- data.frame(mccID1_2@y.values, mccID1_2@x.values)
colnames(mccID1_2) <- c("MCC Score", "Confidence Cutoff")
mccID1_2 <- mccID1_2[2:(nrow(mccID1_2)-1),]
optMCC <- which(mccID1_2$`MCC Score` == max(mccID1_2$`MCC Score`))
mccID1_2 <- mccID1_2$`Confidence Cutoff`[optMCC]
confF1_2$TP <- confF1_2$trueTaxa == confF1_2$assignedTaxa & confF1_2$familyConf >= mccID1_2
confF1_2$TN <- confF1_2$trueTaxa != confF1_2$assignedTaxa & confF1_2$familyConf < mccID1_2
confF1_2$FP <- confF1_2$trueTaxa != confF1_2$assignedTaxa & confF1_2$familyConf >= mccID1_2
confF1_2$FN <- confF1_2$trueTaxa == confF1_2$assignedTaxa & confF1_2$familyConf < mccID1_2

###########
# Family Testing - 2nd run - F2_F3/F1
taxonomy <- paste("Root", fam_f2_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold2_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold1, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- fam_f1
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confF2_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confF2_3$outcome <- confF2_3$trueTaxa == confF2_3$assignedTaxa
mccID2_3 <- prediction(confF2_3$familyConf, confF2_3$outcome)
mccID2_3 <- performance(mccID2_3,"mat")
mccID2_3 <- data.frame(mccID2_3@y.values, mccID2_3@x.values)
colnames(mccID2_3) <- c("MCC Score", "Confidence Cutoff")
mccID2_3 <- mccID2_3[2:(nrow(mccID2_3)-1),]
optMCC <- which(mccID2_3$`MCC Score` == max(mccID2_3$`MCC Score`))
mccID2_3 <- mccID2_3$`Confidence Cutoff`[optMCC]
confF2_3$TP <- confF2_3$trueTaxa == confF2_3$assignedTaxa & confF2_3$familyConf >= mccID2_3
confF2_3$TN <- confF2_3$trueTaxa != confF2_3$assignedTaxa & confF2_3$familyConf < mccID2_3
confF2_3$FP <- confF2_3$trueTaxa != confF2_3$assignedTaxa & confF2_3$familyConf >= mccID2_3
confF2_3$FN <- confF2_3$trueTaxa == confF2_3$assignedTaxa & confF2_3$familyConf < mccID2_3

###########
# Family Testing - 3rd run - F1_3/F2
taxonomy <- paste("Root", fam_f1_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold2, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- fam_f2
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confF1_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confF1_3$outcome <- confF1_3$trueTaxa == confF1_3$assignedTaxa
mccID1_3 <- prediction(confF1_3$familyConf, confF1_3$outcome)
mccID1_3 <- performance(mccID1_3,"mat")
mccID1_3 <- data.frame(mccID1_3@y.values, mccID1_3@x.values)
colnames(mccID1_3) <- c("MCC Score", "Confidence Cutoff")
mccID1_3 <- mccID1_3[2:(nrow(mccID1_3)-1),]
optMCC <- which(mccID1_3$`MCC Score` == max(mccID1_3$`MCC Score`))
mccID1_3 <- mccID1_3$`Confidence Cutoff`[optMCC]
confF1_3$TP <- confF1_3$trueTaxa == confF1_3$assignedTaxa & confF1_3$familyConf >= mccID1_3
confF1_3$TN <- confF1_3$trueTaxa != confF1_3$assignedTaxa & confF1_3$familyConf < mccID1_3
confF1_3$FP <- confF1_3$trueTaxa != confF1_3$assignedTaxa & confF1_3$familyConf >= mccID1_3
confF1_3$FN <- confF1_3$trueTaxa == confF1_3$assignedTaxa & confF1_3$familyConf < mccID1_3

#############
# write.csv(confF1_2, "aranFamConfF1_2.csv")
# write.csv(confF1_3, "aranFamConfF1_3.csv")
# write.csv(confF2_3, "aranFamConfF2_3.csv")

# Accuracy metrics IDTAXA
avgIncNumFamAmph1_2 <- length(which(confF1_2$FP == TRUE)) + length(which(confF1_2$TN == TRUE))
incAssignFamAmph1_2 <- avgIncNumFamAmph1_2 / nrow(confF1_2)
incAssignFamAmph1_2 * 100

avgIncNumFamAmph2_3 <- length(which(confF2_3$FP == TRUE)) + length(which(confF2_3$TN == TRUE))
incAssignFamAmph2_3 <- avgIncNumFamAmph2_3 / nrow(confF2_3)
incAssignFamAmph2_3 * 100

avgIncNumFamAmph1_3 <- length(which(confF1_3$FP == TRUE)) + length(which(confF1_3$TN == TRUE))
incAssignFamAmph1_3 <- avgIncNumFamAmph1_3 / nrow(confF1_3)
incAssignFamAmph1_3 * 100

amphFam <- rbind(confF1_2, confF2_3, confF1_3)

avgIncNumFamAmph <- length(which(amphFam$FP == TRUE)) + length(which(amphFam$TN == TRUE))
incAssignFamAvgAmph <- avgIncNumFamAmph / nrow(amphFam)
incAssignFamAvgAmph * 100

#########
# ROC Curve Setup

# Family - Amph
amphFam$outcome <- as.integer(as.logical(amphFam$trueTaxa == amphFam$assignedTaxa))

#######
# BLAST Tests - Data Import and database creation
# Amphibia BLAST
writeXStringSet(fold1_2, "fold1_2ampF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3ampF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3ampF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampF.fas"))

#######
# Mammalia BLAST

writeXStringSet(fold1_2, "fold1_2mamF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3mamF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3mamF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2mamF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1mamF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamF.fas"))

########
# Gastropda BLAST
writeXStringSet(fold1_2, "fold1_2gasF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3gasF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3gasF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasF.fas"))

########
# Araneae BLAST
writeXStringSet(fold1_2, "fold1_2aranF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3aranF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3aranF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranF.fas"))

########
# Hymenoptera BLAST
writeXStringSet(fold1_2, "fold1_2hymF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3hymF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3hymF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymF.fas"))

########
# Lep BLAST
writeXStringSet(fold1_2, "fold1_2lepF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3lepF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3lepF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepF.fas"))

########
# Aves BLAST
writeXStringSet(fold1_2, "fold1_2avesF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3avesF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3avesF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesF.fas"))

########
# Fish BLAST
writeXStringSet(fold1_2, "fold1_2fishF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3fishF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3fishF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishF.fas"))

########
# Dip BLAST
writeXStringSet(fold1_2, "fold1_2dipF.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3dipF.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3dipF.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipF.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipF.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipF.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipF.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipF.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipF.fas"))

##########
# BLAST predictions
# Use the predict function to query against custom BLAST database
predictF2F3 <- predict(db2_3, fold1)

predictF1F3 <- predict(db1_3, fold2)

predictF1F2 <- predict(db1_2, fold3)

##########
# Family predictions - all runs
# F1F2
familyPredictF1F2 <- strsplit(as.character(predictF1F2$SubjectID), ";")
predictF1F2$assignedTaxa <- sapply(familyPredictF1F2, function(x) x[1])

# F1F3
familyPredictF1F3 <- strsplit(as.character(predictF1F3$SubjectID), ";")
predictF1F3$assignedTaxa <- sapply(familyPredictF1F3, function(x) x[1])

# F2F3
familyPredictF2F3 <- strsplit(as.character(predictF2F3$SubjectID), ";")
predictF2F3$assignedTaxa <- sapply(familyPredictF2F3, function(x) x[1])

# Separate query id from name
familyTrue1 <- strsplit(as.character(predictF1F2$QueryID), ";")
predictF1F2$trueTaxa <- sapply(familyTrue1, function(x) x[1])
predictF1F2$QueryID <- sapply(familyTrue1, function(x) x[2])

familyTrue2 <- strsplit(as.character(predictF1F3$QueryID), ";")
predictF1F3$trueTaxa <- sapply(familyTrue2, function(x) x[1])
predictF1F3$QueryID <- sapply(familyTrue2, function(x) x[2])

familyTrue3 <- strsplit(as.character(predictF2F3$QueryID), ";")
predictF2F3$trueTaxa <- sapply(familyTrue3, function(x) x[1])
predictF2F3$QueryID <- sapply(familyTrue3, function(x) x[2])

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST1 <- data.table(predictF1F2, key="QueryID")
dfBLAST1 <- dfBLAST1[, head(.SD, 1), by=QueryID]
blastF1_2 <- dfBLAST1

# MCC determination for optimal values and TP/TN/FP/FN
blastF1_2$outcome <- blastF1_2$trueTaxa == blastF1_2$assignedTaxa
mccb1_2 <- prediction(blastF1_2$Perc.Ident, blastF1_2$outcome)
mccb1_2 <- performance(mccb1_2,"mat")
mccb1_2 <- data.frame(mccb1_2@y.values, mccb1_2@x.values)
colnames(mccb1_2) <- c("MCC Score", "Confidence Cutoff")
mccb1_2 <- mccb1_2[2:(nrow(mccb1_2)-1),]
optMCCb <- which(mccb1_2$`MCC Score` == max(mccb1_2$`MCC Score`))
mccb1_2 <- mccb1_2$`Confidence Cutoff`[optMCCb]
blastF1_2$TP <- blastF1_2$trueTaxa == blastF1_2$assignedTaxa & blastF1_2$Perc.Ident >= mccb1_2
blastF1_2$TN <- blastF1_2$trueTaxa != blastF1_2$assignedTaxa & blastF1_2$Perc.Ident < mccb1_2
blastF1_2$FP <- blastF1_2$trueTaxa != blastF1_2$assignedTaxa & blastF1_2$Perc.Ident >= mccb1_2
blastF1_2$FN <- blastF1_2$trueTaxa == blastF1_2$assignedTaxa & blastF1_2$Perc.Ident < mccb1_2

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST2 <- data.table(predictF1F3, key="QueryID")
dfBLAST2 <- dfBLAST2[, head(.SD, 1), by=QueryID]
blastF1_3 <- dfBLAST2

# MCC determination for optimal values and TP/TN/FP/FN
blastF1_3$outcome <- blastF1_3$trueTaxa == blastF1_3$assignedTaxa
mccb1_3 <- prediction(blastF1_3$Perc.Ident, blastF1_3$outcome)
mccb1_3 <- performance(mccb1_3,"mat")
mccb1_3 <- data.frame(mccb1_3@y.values, mccb1_3@x.values)
colnames(mccb1_3) <- c("MCC Score", "Confidence Cutoff")
mccb1_3 <- mccb1_3[2:(nrow(mccb1_3)-1),]
optMCCb <- which(mccb1_3$`MCC Score` == max(mccb1_3$`MCC Score`))
mccb1_3 <- mccb1_3$`Confidence Cutoff`[optMCCb]
blastF1_3$TP <- blastF1_3$trueTaxa == blastF1_3$assignedTaxa & blastF1_3$Perc.Ident >= mccb1_3
blastF1_3$TN <- blastF1_3$trueTaxa != blastF1_3$assignedTaxa & blastF1_3$Perc.Ident < mccb1_3
blastF1_3$FP <- blastF1_3$trueTaxa != blastF1_3$assignedTaxa & blastF1_3$Perc.Ident >= mccb1_3
blastF1_3$FN <- blastF1_3$trueTaxa == blastF1_3$assignedTaxa & blastF1_3$Perc.Ident < mccb1_3

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST3 <- data.table(predictF2F3, key="QueryID")
dfBLAST3 <- dfBLAST3[, head(.SD, 1), by=QueryID]
blastF2_3 <- dfBLAST3

# MCC determination for optimal values and TP/TN/FP/FN
blastF2_3$outcome <- blastF2_3$trueTaxa == blastF2_3$assignedTaxa
mccb2_3 <- prediction(blastF2_3$Perc.Ident, blastF2_3$outcome)
mccb2_3 <- performance(mccb2_3,"mat")
mccb2_3 <- data.frame(mccb2_3@y.values, mccb2_3@x.values)
colnames(mccb2_3) <- c("MCC Score", "Confidence Cutoff")
mccb2_3 <- mccb2_3[2:(nrow(mccb2_3)-1),]
optMCCb <- which(mccb2_3$`MCC Score` == max(mccb2_3$`MCC Score`))
mccb2_3 <- mccb2_3$`Confidence Cutoff`[optMCCb]
blastF2_3$TP <- blastF2_3$trueTaxa == blastF2_3$assignedTaxa & blastF2_3$Perc.Ident >= mccb2_3
blastF2_3$TN <- blastF2_3$trueTaxa != blastF2_3$assignedTaxa & blastF2_3$Perc.Ident < mccb2_3
blastF2_3$FP <- blastF2_3$trueTaxa != blastF2_3$assignedTaxa & blastF2_3$Perc.Ident >= mccb2_3
blastF2_3$FN <- blastF2_3$trueTaxa == blastF2_3$assignedTaxa & blastF2_3$Perc.Ident < mccb2_3

rm(predictF1F2)
rm(predictF1F3)
rm(predictF2F3)
rm(familyPredictF1F2)
rm(familyPredictF1F3)
rm(familyPredictF2F3)
rm(familyTrue1)
rm(familyTrue2)
rm(familyTrue3)

#############
# write.csv(blastF1_2, "aranFamblastF1_2.csv")
# write.csv(blastF1_3, "aranFamblastF1_3.csv")
# write.csv(blastF2_3, "aranFamblastF2_3.csv")

# Accuracy metrics BLAST
avgIncNumFamAmph1_2b <- length(which(blastF1_2$FP == TRUE)) + length(which(blastF1_2$TN == TRUE))
incAssignFamAmph1_2b <- avgIncNumFamAmph1_2b / nrow(blastF1_2)
incAssignFamAmph1_2b * 100

avgIncNumFamAmph2_3b <- length(which(blastF2_3$FP == TRUE)) + length(which(blastF2_3$TN == TRUE))
incAssignFamAmph2_3b <- avgIncNumFamAmph2_3b / nrow(blastF2_3)
incAssignFamAmph2_3b * 100

avgIncNumFamAmph1_3b <- length(which(blastF1_3$FP == TRUE)) + length(which(blastF1_3$TN == TRUE))
incAssignFamAmph1_3b <- avgIncNumFamAmph1_3b / nrow(blastF1_3)
incAssignFamAmph1_3b * 100

amphFamB <- rbind(blastF1_2, blastF1_3, blastF2_3)

avgIncNumFamAmphB <- length(which(amphFamB$FP == TRUE)) + length(which(amphFamB$TN == TRUE))
incAssignFamAvgAmphB <- avgIncNumFamAmphB / nrow(amphFamB)
incAssignFamAvgAmphB * 100

#########
# ROC Curve Setup

# Family - Amph
amphFamB$outcome <- as.integer(as.logical(amphFamB$trueTaxa == amphFamB$assignedTaxa))

####################
# RDP Classifier 

# First write out all files to fasta
writeXStringSet(fold1, "fold1amp_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2amp_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3amp_dec8F.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1mam_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2mam_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3mam_dec8F.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1gasF.fas", format = "FASTA", compress = FALSE)
writeXStringSet(fold2, "fold2gasF.fas", format = "FASTA", compress = FALSE)
writeXStringSet(fold3, "fold3gasF.fas", format = "FASTA", compress = FALSE)

writeXStringSet(fold1, "fold1aran_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2aran_dec8F.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3aran_dec8F.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1hymF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2hymF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3hymF.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1lepF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2lepF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3lepF.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1avesF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2avesF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3avesF.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1fishF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2fishF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3fishF.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1dipF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2dipF.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3dipF.fas", format = "FASTA", compress = TRUE)

# Reformat prediction data
# Family

# Converting to old variable name for ease of pipeline running
predictRDPf1_2ampF <- predictRDPf1_2dipF
predictRDPf1_3ampF <- predictRDPf1_3dipF
predictRDPf2_3ampF <- predictRDPf2_3dipF
customRDPf1_2ampF <- customRDPf1_2dipF
customRDPf1_3ampF <- customRDPf1_3dipF
customRDPf2_3ampF <- customRDPf2_3dipF

predictRDPf1_2ampF$assignedTaxa <- predictRDPf1_2ampF$Genus # RDP shifted the taxononmy 1 over so genus is now family according to the labeling
predictRDPf1_3ampF$assignedTaxa <- predictRDPf1_3ampF$Genus
predictRDPf2_3ampF$assignedTaxa <- predictRDPf2_3ampF$Genus
predictRDPf1_2ampF$assignedTaxa[is.na(predictRDPf1_2ampF$assignedTaxa)] <- "Unknown"
predictRDPf1_3ampF$assignedTaxa[is.na(predictRDPf1_3ampF$assignedTaxa)] <- "Unknown"
predictRDPf2_3ampF$assignedTaxa[is.na(predictRDPf2_3ampF$assignedTaxa)] <- "Unknown"

predictRDPf1_2ampF$trueTaxa <- fam_f3
predictRDPf1_3ampF$trueTaxa <- fam_f2
predictRDPf2_3ampF$trueTaxa <- fam_f1

# Convert to percentage
predictRDPf1_2ampF$confF <- customRDPf1_2ampF$Genus * 100  
predictRDPf1_3ampF$confF <- customRDPf1_3ampF$Genus * 100
predictRDPf2_3ampF$confF <- customRDPf2_3ampF$Genus * 100

RDPF1_2F <- predictRDPf1_2ampF

# MCC determination for optimal values and TP/TN/FP/FN
RDPF1_2F$outcome <- RDPF1_2F$trueTaxa == RDPF1_2F$assignedTaxa
mccr1_2 <- prediction(RDPF1_2F$confF, RDPF1_2F$outcome)
mccr1_2 <- performance(mccr1_2,"mat")
mccr1_2 <- data.frame(mccr1_2@y.values, mccr1_2@x.values)
colnames(mccr1_2) <- c("MCC Score", "Confidence Cutoff")
mccr1_2 <- mccr1_2[2:(nrow(mccr1_2)-1),]
optMCCr <- which(mccr1_2$`MCC Score` == max(mccr1_2$`MCC Score`))
mccr1_2 <- mccr1_2$`Confidence Cutoff`[optMCCr]
RDPF1_2F$TP <- RDPF1_2F$trueTaxa == RDPF1_2F$assignedTaxa & RDPF1_2F$confF >= mccr1_2
RDPF1_2F$TN <- RDPF1_2F$trueTaxa != RDPF1_2F$assignedTaxa & RDPF1_2F$confF < mccr1_2
RDPF1_2F$FP <- RDPF1_2F$trueTaxa != RDPF1_2F$assignedTaxa & RDPF1_2F$confF >= mccr1_2
RDPF1_2F$FN <- RDPF1_2F$trueTaxa == RDPF1_2F$assignedTaxa & RDPF1_2F$confF < mccr1_2

# 2nd run
RDPF1_3F <- predictRDPf1_3ampF

RDPF1_3F$outcome <- RDPF1_3F$trueTaxa == RDPF1_3F$assignedTaxa
mccr1_3 <- prediction(RDPF1_3F$confF, RDPF1_3F$outcome)
mccr1_3 <- performance(mccr1_3,"mat")
mccr1_3 <- data.frame(mccr1_3@y.values, mccr1_3@x.values)
colnames(mccr1_3) <- c("MCC Score", "Confidence Cutoff")
mccr1_3 <- mccr1_3[2:(nrow(mccr1_3)-1),]
optMCCr <- which(mccr1_3$`MCC Score` == max(mccr1_3$`MCC Score`))
mccr1_3 <- mccr1_3$`Confidence Cutoff`[optMCCr]
RDPF1_3F$TP <- RDPF1_3F$trueTaxa == RDPF1_3F$assignedTaxa & RDPF1_3F$confF >= mccr1_3
RDPF1_3F$TN <- RDPF1_3F$trueTaxa != RDPF1_3F$assignedTaxa & RDPF1_3F$confF < mccr1_3
RDPF1_3F$FP <- RDPF1_3F$trueTaxa != RDPF1_3F$assignedTaxa & RDPF1_3F$confF >= mccr1_3
RDPF1_3F$FN <- RDPF1_3F$trueTaxa == RDPF1_3F$assignedTaxa & RDPF1_3F$confF < mccr1_3

# 3rd run
RDPF2_3F <- predictRDPf2_3ampF

# MCC determination for optimal values and TP/TN/FP/FN
RDPF2_3F$outcome <- RDPF2_3F$trueTaxa == RDPF2_3F$assignedTaxa
mccr2_3 <- prediction(RDPF2_3F$confF, RDPF2_3F$outcome)
mccr2_3 <- performance(mccr2_3,"mat")
mccr2_3 <- data.frame(mccr2_3@y.values, mccr2_3@x.values)
colnames(mccr2_3) <- c("MCC Score", "Confidence Cutoff")
mccr2_3 <- mccr2_3[2:(nrow(mccr2_3)-1),]
optMCCr <- which(mccr2_3$`MCC Score` == max(mccr2_3$`MCC Score`))
mccr2_3 <- mccr2_3$`Confidence Cutoff`[optMCCr]
RDPF2_3F$TP <- RDPF2_3F$trueTaxa == RDPF2_3F$assignedTaxa & RDPF2_3F$confF >= mccr2_3
RDPF2_3F$TN <- RDPF2_3F$trueTaxa != RDPF2_3F$assignedTaxa & RDPF2_3F$confF < mccr2_3
RDPF2_3F$FP <- RDPF2_3F$trueTaxa != RDPF2_3F$assignedTaxa & RDPF2_3F$confF >= mccr2_3
RDPF2_3F$FN <- RDPF2_3F$trueTaxa == RDPF2_3F$assignedTaxa & RDPF2_3F$confF < mccr2_3

#############
# write.csv(RDPF1_2F, "aranFamRDPF1_2F.csv")
# write.csv(RDPF1_3F, "aranFamRDPF1_3F.csv")
# write.csv(RDPF2_3F, "aranFamRDPF2_3F.csv")

# Accuracy metrics RDP
avgIncNumFamAmph1_2r <- length(which(RDPF1_2F$FP == TRUE)) + length(which(RDPF1_2F$TN == TRUE))
incAssignFamAmph1_2r <- avgIncNumFamAmph1_2r / nrow(RDPF1_2F)
incAssignFamAmph1_2r * 100

avgIncNumFamAmph2_3r <- length(which(RDPF2_3F$FP == TRUE)) + length(which(RDPF2_3F$TN == TRUE))
incAssignFamAmph2_3r <- avgIncNumFamAmph2_3r / nrow(RDPF2_3F)
incAssignFamAmph2_3r * 100

avgIncNumFamAmph1_3r <- length(which(RDPF1_3F$FP == TRUE)) + length(which(RDPF1_3F$TN == TRUE))
incAssignFamAmph1_3r <- avgIncNumFamAmph1_3r / nrow(RDPF1_3F)
incAssignFamAmph1_3r * 100

amphFamRDP <- rbind(RDPF1_2F, RDPF1_3F, RDPF2_3F)

avgIncNumFamAmphR <- length(which(amphFamRDP$FP == TRUE)) + length(which(amphFamRDP$TN == TRUE))
incAssignFamAvgAmphR <- avgIncNumFamAmphR / nrow(amphFamRDP)
incAssignFamAvgAmphR * 100

#########
# ROC Curve Setup

# Family 
amphFamRDP$outcome <- as.integer(as.logical(amphFamRDP$trueTaxa == amphFamRDP$assignedTaxa))
amphFam$outcome <- as.integer(as.logical(amphFam$trueTaxa == amphFam$assignedTaxa))
