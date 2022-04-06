# install.packages("foreach")
library(foreach)
# For sequence alignments we need the biostrings (DNAStringSet function) 
# and muscle libraries, as follows:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings) 
# BiocManager::install("DECIPHER")
library(DECIPHER)
library(dplyr)
library(ggplot2)
# For rBLAST
# install.packages("remotes")
# library(remotes)
# remotes::install_github("mhahsler/rBLAST")
library(rBLAST)
# For plotting
# install.packages("viridis")
library(viridis)
library(ggplot2)
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
library(RColorBrewer)
library(cowplot)
# install.packages("ggpubr")
library(ggpubr)

# For getting fasta files I need in FASTA
# writeXStringSet(dna, ".fas", format="FASTA")

########
# Groups this pipeline is testing currently:

# Amphibibia - Family, genus and species level benchmarking 

# Mammalia - Family, genus and species level benchmarking 

# Gastropoda - Family, genus and species level benchmarking 

# Araneae - Family, genus and species level benchmarking 

# Hymenoptera - Family, genus and species level benchmarking 

#########
# Dataset download

# Amphibia
dna <- readDNAStringSet("amphibiaDNAUnfiltered.fas", format="FASTA")

# Mammalia
dna <- readDNAStringSet("mammaliaDNAUnfiltered.fas", format="FASTA")

# Aves
genus <- sapply(dnaNames, `[`, 3)
species <- sapply(dnaNames, `[`, 4)
species <- gsub(" ", "_", species)

names(dna) <- paste(genus, species, sep=";")

dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$gen <- genus
dfLep$sp <- species
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

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
dfLep$gen <- genus
dfLep$sp <- species

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

#########
# Start of pipeline - Same process for all groups (repeated for each group until ROC and MCC plots)

#########
# For porter datasets
dnaNames <- names(dna)

dnaNames <- strsplit(dnaNames, " ")

# For species level
spID <- sapply(dnaNames, `[`, 2)
spID <- strsplit(spID, ";")
spID <- sapply(spID, `[`, 9)

# Start of stratified sampling
dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$names <- familyID

# Remove duplicate sequences
dfLep$sp <- spID
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

# Picked 5 largest families
# largestFam <- tail(sort(table(dfLep$names)), 5)
# dfLepSub <- foreach(i=1:nrow(largestFam)) %do% which(dfLep$names == names(largestFam)[i])
# dfLepSub <- do.call(c, dfLepSub)
# dfLepSub <- dfLep[dfLepSub,]

########

# Number of species
length(table(unique(dfLep$sp)))

# Break down by family identifier (or other taxonomic rank for full sequence datasets)
nameList <- lapply(unique(dfLep$sp), function(x) 
  dfLep[dfLep$sp == x,])

# Conversion to DNAStringSet and naming with family and sequence id
nameList2 <- sapply( nameList , function(x) DNAStringSet( x$dna ) )

for (i in seq(from = 1, to = length(nameList), by = 1)) {
  names(nameList2[[i]]) <- paste(nameList[[i]]$sp, nameList[[i]]$id, sep=";")
}

# Removal of gaps in sequences before testing
nameList2 <- foreach(i=1:length(nameList2)) %do% RemoveGaps(nameList2[[i]],
                                                            removeGaps = "all",
                                                            processors = 1)

#############
# For lep to reduce number of species
speciesNum <- sapply( nameList2 , function(x) length( x@ranges@NAMES ) )
speciesLarge <- which(speciesNum >= 4)
nameList3 <- nameList2[c(speciesLarge)]
allFolds <- foreach(i=1:length(nameList3)) %do% createFolds(nameList3[[i]], k=3, list = TRUE)

# Fold 1
fold1 <- foreach(i=1:length(nameList3)) %do% nameList3[[i]][allFolds[[i]]$Fold1]
fold1 <- do.call(c, fold1)
# Fold 2
fold2 <- foreach(i=1:length(nameList3)) %do% nameList3[[i]][allFolds[[i]]$Fold2]
fold2 <- do.call(c, fold2)
# Fold 3
fold3 <- foreach(i=1:length(nameList3)) %do% nameList3[[i]][allFolds[[i]]$Fold3]
fold3 <- do.call(c, fold3)

# Save true taxon and ids for testing sequences
fam_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

fam_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

fam_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

###########
# For lep to cut down on larger genus
df_f1 <- data.frame(nucleotides=fold1, name=fam_f1, id=famids_f1)
df_f1 <- data.table(df_f1, key="name")
df_f1 <- df_f1[, head(.SD, 2), by=name]

df_f2 <- data.frame(nucleotides=fold2, name=fam_f2, id=famids_f2)
df_f2 <- data.table(df_f2, key="name")
df_f2 <- df_f2[, head(.SD, 2), by=name]

df_f3 <- data.frame(nucleotides=fold3, name=fam_f3, id=famids_f3)
df_f3 <- data.table(df_f3, key="name")
df_f3 <- df_f3[, head(.SD, 2), by=name]

fold1 <- DNAStringSet(df_f1$nucleotides)
fold2 <- DNAStringSet(df_f2$nucleotides)
fold3 <- DNAStringSet(df_f3$nucleotides)

names(fold1) <- paste(df_f1$name, df_f1$id, sep=";")
names(fold2) <- paste(df_f2$name, df_f2$id, sep=";")
names(fold3) <- paste(df_f3$name, df_f3$id, sep=";")

#############

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

# Folds 1+2 for testing
fold1_2 <- append(fold1, fold2)

# Folds 2+3
fold2_3 <- append(fold2, fold3)

# Folds 1+3
fold1_3 <- append(fold1, fold3)

#######
# IDTAXA Testing (Murali at al., 2019)
# 1st run - F1/F2+F3, 2nd run - F2/F1+F3, 3rd run - F3/F1+F2

# Save true taxon and ids for testing sequences
sp_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

sp_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

sp_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

# Combined for training sequences
sp_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 1)
famids_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 2)

sp_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 1)
famids_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 2)

sp_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 1)
famids_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 2)

#########
# Species Testing - 1st run - F1_F2/F3
taxonomy <- paste("Root",sp_f1_2, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_2, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold3, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- sp_f3
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confS1_2 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confS1_2$outcome <- confS1_2$trueTaxa == confS1_2$assignedTaxa
mccID1_2 <- prediction(confS1_2$familyConf, confS1_2$outcome)
mccID1_2 <- performance(mccID1_2,"mat")
mccID1_2 <- data.frame(mccID1_2@y.values, mccID1_2@x.values)
colnames(mccID1_2) <- c("MCC Score", "Confidence Cutoff")
mccID1_2 <- mccID1_2[2:(nrow(mccID1_2)-1),]
optMCC <- which(mccID1_2$`MCC Score` == max(mccID1_2$`MCC Score`))
mccID1_2 <- mccID1_2$`Confidence Cutoff`[optMCC]
confS1_2$TP <- confS1_2$trueTaxa == confS1_2$assignedTaxa & confS1_2$familyConf >= mccID1_2
confS1_2$TN <- confS1_2$trueTaxa != confS1_2$assignedTaxa & confS1_2$familyConf < mccID1_2
confS1_2$FP <- confS1_2$trueTaxa != confS1_2$assignedTaxa & confS1_2$familyConf >= mccID1_2
confS1_2$FN <- confS1_2$trueTaxa == confS1_2$assignedTaxa & confS1_2$familyConf < mccID1_2

###########
# Species Testing - 2nd run - F2_F3/F1
taxonomy <- paste("Root", sp_f2_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold2_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold1, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- sp_f1
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confS2_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confS2_3$outcome <- confS2_3$trueTaxa == confS2_3$assignedTaxa
mccID2_3 <- prediction(confS2_3$familyConf, confS2_3$outcome)
mccID2_3 <- performance(mccID2_3,"mat")
mccID2_3 <- data.frame(mccID2_3@y.values, mccID2_3@x.values)
colnames(mccID2_3) <- c("MCC Score", "Confidence Cutoff")
mccID2_3 <- mccID2_3[2:(nrow(mccID2_3)-1),]
optMCC <- which(mccID2_3$`MCC Score` == max(mccID2_3$`MCC Score`))
mccID2_3 <- mccID2_3$`Confidence Cutoff`[optMCC]
confS2_3$TP <- confS2_3$trueTaxa == confS2_3$assignedTaxa & confS2_3$familyConf >= mccID2_3
confS2_3$TN <- confS2_3$trueTaxa != confS2_3$assignedTaxa & confS2_3$familyConf < mccID2_3
confS2_3$FP <- confS2_3$trueTaxa != confS2_3$assignedTaxa & confS2_3$familyConf >= mccID2_3
confS2_3$FN <- confS2_3$trueTaxa == confS2_3$assignedTaxa & confS2_3$familyConf < mccID2_3

###########
# Species Testing - 3rd run - F1_3/F2
taxonomy <- paste("Root", sp_f1_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold2, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- sp_f2
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confS1_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confS1_3$outcome <- confS1_3$trueTaxa == confS1_3$assignedTaxa
mccID1_3 <- prediction(confS1_3$familyConf, confS1_3$outcome)
mccID1_3 <- performance(mccID1_3,"mat")
mccID1_3 <- data.frame(mccID1_3@y.values, mccID1_3@x.values)
colnames(mccID1_3) <- c("MCC Score", "Confidence Cutoff")
mccID1_3 <- mccID1_3[2:(nrow(mccID1_3)-1),]
optMCC <- which(mccID1_3$`MCC Score` == max(mccID1_3$`MCC Score`))
mccID1_3 <- mccID1_3$`Confidence Cutoff`[optMCC]
confS1_3$TP <- confS1_3$trueTaxa == confS1_3$assignedTaxa & confS1_3$familyConf >= mccID1_3
confS1_3$TN <- confS1_3$trueTaxa != confS1_3$assignedTaxa & confS1_3$familyConf < mccID1_3
confS1_3$FP <- confS1_3$trueTaxa != confS1_3$assignedTaxa & confS1_3$familyConf >= mccID1_3
confS1_3$FN <- confS1_3$trueTaxa == confS1_3$assignedTaxa & confS1_3$familyConf < mccID1_3

#############
# Accuracy metrics IDTAXA
avgIncNumFamAmph1_2 <- length(which(confS1_2$FP == TRUE)) + length(which(confS1_2$TN == TRUE))
incAssignFamAmph1_2 <- avgIncNumFamAmph1_2 / nrow(confS1_2)
incAssignFamAmph1_2 * 100

avgIncNumFamAmph2_3 <- length(which(confS2_3$FP == TRUE)) + length(which(confS2_3$TN == TRUE))
incAssignFamAmph2_3 <- avgIncNumFamAmph2_3 / nrow(confS2_3)
incAssignFamAmph2_3 * 100

avgIncNumFamAmph1_3 <- length(which(confS1_3$FP == TRUE)) + length(which(confS1_3$TN == TRUE))
incAssignFamAmph1_3 <- avgIncNumFamAmph1_3 / nrow(confS1_3)
incAssignFamAmph1_3 * 100

amphSp <- rbind(confS1_2, confS2_3, confS1_3)

avgIncNumSpAmph <- length(which(amphSp$FP == TRUE)) + length(which(amphSp$TN == TRUE))
incAssignSpAvgAmph <- avgIncNumSpAmph / nrow(amphSp)
incAssignSpAvgAmph * 100

#########
# ROC Curve Setup

# Species - Amph
amphSp$outcome <- as.integer(as.logical(amphSp$trueTaxa == amphSp$assignedTaxa))

#######
# BLAST Tests - Data Import and database creation
# Amphibia BLAST
writeXStringSet(fold1_2, "fold1_2ampS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3ampS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3ampS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampS.fas"))

#######
# Mammalia BLAST

writeXStringSet(fold1_2, "fold1_2mamS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3mamS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3mamS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2mamS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2mamS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamS.fas"))

########
# Gastropda BLAST
writeXStringSet(fold1_2, "fold1_2gasS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3gasS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3gasS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasS.fas"))

########
# Araneae BLAST
names(fold1) <- gsub("\\s", "_", names(fold1))
names(fold2) <- gsub("\\s", "_", names(fold2))
names(fold3) <- gsub("\\s", "_", names(fold3))
names(fold1_2) <- gsub("\\s", "_", names(fold1_2))
names(fold1_3) <- gsub("\\s", "_", names(fold1_3))
names(fold2_3) <- gsub("\\s", "_", names(fold2_3))

writeXStringSet(fold1_2, "fold1_2aranS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3aranS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3aranS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranS.fas"))

########
# Diptera BLAST
writeXStringSet(fold1_2, "fold1_2dipS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3dipS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3dipS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipS.fas"))

#########
# Anthophila BLAST
writeXStringSet(fold1_2, "fold1_2anthoS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3anthoS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3anthoS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2anthoS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2anthoS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3anthoS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3anthoS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3anthoS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3anthoS.fas"))

########
# Fish BLAST
writeXStringSet(fold1_2, "fold1_2fishS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3fishS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3fishS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishS.fas"))

###########
# Hymenoptera BLAST
writeXStringSet(fold1_2, "fold1_2hymS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3hymS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3hymS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymS.fas"))

###########
# Lep BLAST
writeXStringSet(fold1_2, "fold1_2lepS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3lepS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3lepS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepS.fas"))

###########
# Aves BLAST
writeXStringSet(fold1_2, "fold1_2avesS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3avesS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3avesS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesS.fas"))

###########
names(fold1_3) <- gsub(" ", "_", names(fold1_3))
names(fold1_2) <- gsub(" ", "_", names(fold1_2))
names(fold2_3) <- gsub(" ", "_", names(fold2_3))

# Bees BLAST 
writeXStringSet(fold1_2, "fold1_2beeS.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3beeS.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3beeS.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2beeS.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2beeS.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3beeS.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3beeS.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3beeS.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3beeS.fas"))

##########
# BLAST predictions
# Use the predict function to query against custom BLAST database
predictF2F3 <- predict(db2_3, fold1)

predictF1F3 <- predict(db1_3, fold2)

predictF1F2 <- predict(db1_2, fold3)

###########
# Species predictions - all runs

names(fold3) <- gsub(" ", "_", names(fold3))
names(fold2) <- gsub(" ", "_", names(fold2))
names(fold1) <- gsub(" ", "_", names(fold1))

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
blastS1_2 <- dfBLAST1

blastS1_2$outcome <- blastS1_2$trueTaxa == blastS1_2$assignedTaxa
mccb1_2 <- prediction(blastS1_2$Perc.Ident, blastS1_2$outcome)
mccb1_2 <- performance(mccb1_2,"mat")
mccb1_2 <- data.frame(mccb1_2@y.values, mccb1_2@x.values)
colnames(mccb1_2) <- c("MCC Score", "Confidence Cutoff")
mccb1_2 <- mccb1_2[2:(nrow(mccb1_2)-1),]
optMCCb <- which(mccb1_2$`MCC Score` == max(mccb1_2$`MCC Score`))
mccb1_2 <- mccb1_2$`Confidence Cutoff`[optMCCb]
blastS1_2$TP <- blastS1_2$trueTaxa == blastS1_2$assignedTaxa & blastS1_2$Perc.Ident >= mccb1_2
blastS1_2$TN <- blastS1_2$trueTaxa != blastS1_2$assignedTaxa & blastS1_2$Perc.Ident < mccb1_2
blastS1_2$FP <- blastS1_2$trueTaxa != blastS1_2$assignedTaxa & blastS1_2$Perc.Ident >= mccb1_2
blastS1_2$FN <- blastS1_2$trueTaxa == blastS1_2$assignedTaxa & blastS1_2$Perc.Ident < mccb1_2

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST2 <- data.table(predictF1F3, key="QueryID")
dfBLAST2 <- dfBLAST2[, head(.SD, 1), by=QueryID]
blastS1_3 <- dfBLAST2

# MCC determination for optimal values and TP/TN/FP/FN
blastS1_3$outcome <- blastS1_3$trueTaxa == blastS1_3$assignedTaxa
mccb1_3 <- prediction(blastS1_3$Perc.Ident, blastS1_3$outcome)
mccb1_3 <- performance(mccb1_3,"mat")
mccb1_3 <- data.frame(mccb1_3@y.values, mccb1_3@x.values)
colnames(mccb1_3) <- c("MCC Score", "Confidence Cutoff")
mccb1_3 <- mccb1_3[2:(nrow(mccb1_3)-1),]
optMCCb <- which(mccb1_3$`MCC Score` == max(mccb1_3$`MCC Score`))
mccb1_3 <- mccb1_3$`Confidence Cutoff`[optMCCb]
blastS1_3$TP <- blastS1_3$trueTaxa == blastS1_3$assignedTaxa & blastS1_3$Perc.Ident >= mccb1_3
blastS1_3$TN <- blastS1_3$trueTaxa != blastS1_3$assignedTaxa & blastS1_3$Perc.Ident < mccb1_3
blastS1_3$FP <- blastS1_3$trueTaxa != blastS1_3$assignedTaxa & blastS1_3$Perc.Ident >= mccb1_3
blastS1_3$FN <- blastS1_3$trueTaxa == blastS1_3$assignedTaxa & blastS1_3$Perc.Ident < mccb1_3

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST3 <- data.table(predictF2F3, key="QueryID")
dfBLAST3 <- dfBLAST3[, head(.SD, 1), by=QueryID]
blastS2_3 <- dfBLAST3

# MCC determination for optimal values and TP/TN/FP/FN
blastS2_3$outcome <- blastS2_3$trueTaxa == blastS2_3$assignedTaxa
mccb2_3 <- prediction(blastS2_3$Perc.Ident, blastS2_3$outcome)
mccb2_3 <- performance(mccb2_3,"mat")
mccb2_3 <- data.frame(mccb2_3@y.values, mccb2_3@x.values)
colnames(mccb2_3) <- c("MCC Score", "Confidence Cutoff")
mccb2_3 <- mccb2_3[2:(nrow(mccb2_3)-1),]
optMCCb <- which(mccb2_3$`MCC Score` == max(mccb2_3$`MCC Score`))
mccb2_3 <- mccb2_3$`Confidence Cutoff`[optMCCb]
blastS2_3$TP <- blastS2_3$trueTaxa == blastS2_3$assignedTaxa & blastS2_3$Perc.Ident >= mccb2_3
blastS2_3$TN <- blastS2_3$trueTaxa != blastS2_3$assignedTaxa & blastS2_3$Perc.Ident < mccb2_3
blastS2_3$FP <- blastS2_3$trueTaxa != blastS2_3$assignedTaxa & blastS2_3$Perc.Ident >= mccb2_3
blastS2_3$FN <- blastS2_3$trueTaxa == blastS2_3$assignedTaxa & blastS2_3$Perc.Ident < mccb2_3

rm(familyPredictF1F2)
rm(familyPredictF1F3)
rm(familyPredictF2F3)
rm(familyTrue1)
rm(familyTrue2)
rm(familyTrue3)
#############
# Accuracy metrics BLAST
avgIncNumFamAmph1_2b <- length(which(blastS1_2$FP == TRUE)) + length(which(blastS1_2$TN == TRUE))
incAssignFamAmph1_2b <- avgIncNumFamAmph1_2b / nrow(blastS1_2)
incAssignFamAmph1_2b * 100

avgIncNumFamAmph2_3b <- length(which(blastS2_3$FP == TRUE)) + length(which(blastS2_3$TN == TRUE))
incAssignFamAmph2_3b <- avgIncNumFamAmph2_3b / nrow(blastS2_3)
incAssignFamAmph2_3b * 100

avgIncNumFamAmph1_3b <- length(which(blastS1_3$FP == TRUE)) + length(which(blastS1_3$TN == TRUE))
incAssignFamAmph1_3b <- avgIncNumFamAmph1_3b / nrow(blastS1_3)
incAssignFamAmph1_3b * 100

amphSpB <- rbind(blastS1_2, blastS1_3, blastS2_3)

avgIncNumSpAmphB <- length(which(amphSpB$FP == TRUE)) + length(which(amphSpB$TN == TRUE))
incAssignSpAvgAmphB <- avgIncNumSpAmphB / nrow(amphSpB)
incAssignSpAvgAmphB * 100

#########
# ROC Curve Setup

# Species - Amph
amphSpB$outcome <- as.integer(as.logical(amphSpB$trueTaxa == amphSpB$assignedTaxa))

####################
# RDP Classifier 

# First write out all files to fasta
writeXStringSet(fold1, "fold1ampS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2ampS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3ampS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1mamS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2mamS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3mamS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1gasS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2gasS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3gasS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1aranS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2aranS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3aranS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1dipS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2dipS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3dipS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1anthoS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2anthoS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3anthoS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1fishS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2fishS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3fishS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1hymS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2hymS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3hymS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1lepS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2lepS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3lepS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1avesS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2avesS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3avesS.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1beeS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2beeS.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3beeS.fas", format = "FASTA", compress = TRUE)

# Reformat prediction data
predictRDPf1_2ampS <- predictRDPf1_2beeS
predictRDPf1_3ampS <- predictRDPf1_3beeS
predictRDPf2_3ampS <- predictRDPf2_3beeS
customRDPf1_2ampS <- customRDPf1_2beeS
customRDPf1_3ampS <- customRDPf1_3beeS
customRDPf2_3ampS <- customRDPf2_3beeS

# Species
predictRDPf1_2ampS$assignedTaxaS <-  gsub(" ", "_", predictRDPf1_2ampS$Genus) # RDP shifted the taxonomy 1 over so genus is now species according to the labeling
predictRDPf1_3ampS$assignedTaxaS <-  gsub(" ", "_", predictRDPf1_3ampS$Genus)
predictRDPf2_3ampS$assignedTaxaS <-  gsub(" ", "_", predictRDPf2_3ampS$Genus)
predictRDPf1_2ampS$assignedTaxaS[is.na(predictRDPf1_2ampS$assignedTaxaS)] <- "Unknown"
predictRDPf1_3ampS$assignedTaxaS[is.na(predictRDPf1_3ampS$assignedTaxaS)] <- "Unknown"
predictRDPf2_3ampS$assignedTaxaS[is.na(predictRDPf2_3ampS$assignedTaxaS)] <- "Unknown"

# If species is not underscored
sp_f3 <- gsub(" ", "_", sp_f3)
sp_f2 <- gsub(" ", "_", sp_f2)
sp_f1 <- gsub(" ", "_", sp_f1)

predictRDPf1_2ampS$trueTaxaS <- sp_f3
predictRDPf1_3ampS$trueTaxaS <- sp_f2
predictRDPf2_3ampS$trueTaxaS <- sp_f1

predictRDPf1_2ampS$confS <- customRDPf1_2ampS$Genus * 100  
predictRDPf1_3ampS$confS <- customRDPf1_3ampS$Genus * 100
predictRDPf2_3ampS$confS <- customRDPf2_3ampS$Genus * 100
RDPF1_2S <- predictRDPf1_2ampS

# MCC determination for optimal values and TP/TN/FP/FN
RDPF1_2S$outcome <- RDPF1_2S$trueTaxaS == RDPF1_2S$assignedTaxaS
mccr1_2 <- prediction(RDPF1_2S$confS, RDPF1_2S$outcome)
mccr1_2 <- performance(mccr1_2,"mat")
mccr1_2 <- data.frame(mccr1_2@y.values, mccr1_2@x.values)
colnames(mccr1_2) <- c("MCC Score", "Confidence Cutoff")
mccr1_2 <- mccr1_2[2:(nrow(mccr1_2)-1),]
optMCCr <- which(mccr1_2$`MCC Score` == max(mccr1_2$`MCC Score`))
mccr1_2 <- mccr1_2$`Confidence Cutoff`[optMCCr]
RDPF1_2S$TP <- RDPF1_2S$trueTaxaS == RDPF1_2S$assignedTaxaS & RDPF1_2S$confS >= mccr1_2
RDPF1_2S$TN <- RDPF1_2S$trueTaxaS != RDPF1_2S$assignedTaxaS & RDPF1_2S$confS < mccr1_2
RDPF1_2S$FP <- RDPF1_2S$trueTaxaS != RDPF1_2S$assignedTaxaS & RDPF1_2S$confS >= mccr1_2
RDPF1_2S$FN <- RDPF1_2S$trueTaxaS == RDPF1_2S$assignedTaxaS & RDPF1_2S$confS < mccr1_2

# 2nd run
RDPF1_3S <- predictRDPf1_3ampS

# MCC determination for optimal values and TP/TN/FP/FN
RDPF1_3S$outcome <- RDPF1_3S$trueTaxaS == RDPF1_3S$assignedTaxaS
mccr1_3 <- prediction(RDPF1_3S$confS, RDPF1_3S$outcome)
mccr1_3 <- performance(mccr1_3,"mat")
mccr1_3 <- data.frame(mccr1_3@y.values, mccr1_3@x.values)
colnames(mccr1_3) <- c("MCC Score", "Confidence Cutoff")
mccr1_3 <- mccr1_3[2:(nrow(mccr1_3)-1),]
optMCCr <- which(mccr1_3$`MCC Score` == max(mccr1_3$`MCC Score`))
mccr1_3 <- mccr1_3$`Confidence Cutoff`[optMCCr]
RDPF1_3S$TP <- RDPF1_3S$trueTaxaS == RDPF1_3S$assignedTaxaS & RDPF1_3S$confS >= mccr1_3
RDPF1_3S$TN <- RDPF1_3S$trueTaxaS != RDPF1_3S$assignedTaxaS & RDPF1_3S$confS < mccr1_3
RDPF1_3S$FP <- RDPF1_3S$trueTaxaS != RDPF1_3S$assignedTaxaS & RDPF1_3S$confS >= mccr1_3
RDPF1_3S$FN <- RDPF1_3S$trueTaxaS == RDPF1_3S$assignedTaxaS & RDPF1_3S$confS < mccr1_3

# 3rd run
# TP Values 
RDPF2_3S <- predictRDPf2_3ampS

RDPF2_3S$outcome <- RDPF2_3S$trueTaxaS == RDPF2_3S$assignedTaxaS
mccr2_3 <- prediction(RDPF2_3S$confS, RDPF2_3S$outcome)
mccr2_3 <- performance(mccr2_3,"mat")
mccr2_3 <- data.frame(mccr2_3@y.values, mccr2_3@x.values)
colnames(mccr2_3) <- c("MCC Score", "Confidence Cutoff")
mccr2_3 <- mccr2_3[2:(nrow(mccr2_3)-1),]
optMCCr <- which(mccr2_3$`MCC Score` == max(mccr2_3$`MCC Score`))
mccr2_3 <- mccr2_3$`Confidence Cutoff`[optMCCr]
RDPF2_3S$TP <- RDPF2_3S$trueTaxaS == RDPF2_3S$assignedTaxaS & RDPF2_3S$confS >= mccr2_3
RDPF2_3S$TN <- RDPF2_3S$trueTaxaS != RDPF2_3S$assignedTaxaS & RDPF2_3S$confS < mccr2_3
RDPF2_3S$FP <- RDPF2_3S$trueTaxaS != RDPF2_3S$assignedTaxaS & RDPF2_3S$confS >= mccr2_3
RDPF2_3S$FN <- RDPF2_3S$trueTaxaS == RDPF2_3S$assignedTaxaS & RDPF2_3S$confS < mccr2_3

#############
# Accuracy metrics RDP

avgIncNumFamAmph1_2r <- length(which(RDPF1_2S$FP == TRUE)) + length(which(RDPF1_2S$TN == TRUE))
incAssignFamAmph1_2r <- avgIncNumFamAmph1_2r / nrow(RDPF1_2S)
incAssignFamAmph1_2r * 100

avgIncNumFamAmph2_3r <- length(which(RDPF2_3S$FP == TRUE)) + length(which(RDPF2_3S$TN == TRUE))
incAssignFamAmph2_3r <- avgIncNumFamAmph2_3r / nrow(RDPF2_3S)
incAssignFamAmph2_3r * 100

avgIncNumFamAmph1_3r <- length(which(RDPF1_3S$FP == TRUE)) + length(which(RDPF1_3S$TN == TRUE))
incAssignFamAmph1_3r <- avgIncNumFamAmph1_3r / nrow(RDPF1_3S)
incAssignFamAmph1_3r * 100

amphSpRDP <- rbind(RDPF1_2S, RDPF1_3S, RDPF2_3S)

avgIncNumSpAmphR <- length(which(amphSpRDP$FP == TRUE)) + length(which(amphSpRDP$TN == TRUE))
incAssignSpAvgAmphR <- avgIncNumSpAmphR / nrow(amphSpRDP)
incAssignSpAvgAmphR * 100

#########
# ROC Curve Setup

# Species - Amph
amphSpRDP$outcomeS <- as.integer(as.logical(amphSpRDP$trueTaxaS == amphSpRDP$assignedTaxaS))
