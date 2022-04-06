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

# Anthophila
dnaNames <- strsplit(names(dna), "_")
genus <- sapply(dnaNames, `[`, 3)
species <- sapply(dnaNames, `[`, 4)

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

# Diptera/Anthophila/Fish
dnaNames <- strsplit(names(dna), "_")

genus <- sapply(dnaNames, `[`, 3)
species <- sapply(dnaNames, `[`, 4)

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

# For genus level
genID <- sapply(dnaNames, `[`, 2)
genID <- strsplit(genID, ";")
genID <- sapply(genID, `[`, 8)

# Start of stratified sampling
dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"

# Remove duplicate sequences
dfLep$gen <- genID
dfLep <- dfLep[!duplicated(dfLep[c("dna")]),]
dfLep$id <- 1:nrow(dfLep)

# Picked 5 largest families
# largestFam <- tail(sort(table(dfLep$names)), 5)
# dfLepSub <- foreach(i=1:nrow(largestFam)) %do% which(dfLep$names == names(largestFam)[i])
# dfLepSub <- do.call(c, dfLepSub)
# dfLepSub <- dfLep[dfLepSub,]

#########

# Fish
dfLep$gen <- dfLep$genus
dfLep$sp <- dfLep$species

# Total num of sequences used
nrow(dfLep)

# Number of genus
length(table(unique(dfLep$gen)))

# Break down by family identifier (or other taxonomic rank for full sequence datasets)
nameList <- lapply(unique(dfLep$gen), function(x) 
  dfLep[dfLep$gen == x,])

# Conversion to DNAStringSet and naming with family and sequence id
nameList2 <- sapply( nameList , function(x) DNAStringSet( x$dna ) )

for (i in seq(from = 1, to = length(nameList), by = 1)) {
  names(nameList2[[i]]) <- paste(nameList[[i]]$gen, nameList[[i]]$id, sep=";")
}

# Removal of gaps in sequences before testing
nameList2 <- foreach(i=1:length(nameList2)) %do% RemoveGaps(nameList2[[i]],
                                                            removeGaps = "all",
                                                            processors = 1)
#############
# For lep to reduce number of genuses
genusNum <- sapply( nameList2 , function(x) length( x@ranges@NAMES ) )
genusLarge <- which(genusNum > 6)
nameList3 <- nameList2[c(genusLarge)]
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
df_f1 <- df_f1[, head(.SD, 6), by=name]

df_f2 <- data.frame(nucleotides=fold2, name=fam_f2, id=famids_f2)
df_f2 <- data.table(df_f2, key="name")
df_f2 <- df_f2[, head(.SD, 6), by=name]

df_f3 <- data.frame(nucleotides=fold3, name=fam_f3, id=famids_f3)
df_f3 <- data.table(df_f3, key="name")
df_f3 <- df_f3[, head(.SD, 6), by=name]

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
gen_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

gen_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

gen_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

# Combined for training sequences
gen_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 1)
famids_f1_2 <- sapply(strsplit(names(fold1_2), ";"), `[`, 2)

gen_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 1)
famids_f2_3 <- sapply(strsplit(names(fold2_3), ";"), `[`, 2)

gen_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 1)
famids_f1_3 <- sapply(strsplit(names(fold1_3), ";"), `[`, 2)

#########
# Genus Testing - 1st run - F1_F2/F3
taxonomy <- paste("Root",gen_f1_2, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_2, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold3, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- gen_f3
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confG1_2 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confG1_2$outcome <- confG1_2$trueTaxa == confG1_2$assignedTaxa
mccID1_2 <- prediction(confG1_2$familyConf, confG1_2$outcome)
mccID1_2 <- performance(mccID1_2,"mat")
mccID1_2 <- data.frame(mccID1_2@y.values, mccID1_2@x.values)
colnames(mccID1_2) <- c("MCC Score", "Confidence Cutoff")
mccID1_2 <- mccID1_2[2:(nrow(mccID1_2)-1),]
optMCC <- which(mccID1_2$`MCC Score` == max(mccID1_2$`MCC Score`))
mccID1_2 <- mccID1_2$`Confidence Cutoff`[optMCC]
confG1_2$TP <- confG1_2$trueTaxa == confG1_2$assignedTaxa & confG1_2$familyConf >= mccID1_2
confG1_2$TN <- confG1_2$trueTaxa != confG1_2$assignedTaxa & confG1_2$familyConf < mccID1_2
confG1_2$FP <- confG1_2$trueTaxa != confG1_2$assignedTaxa & confG1_2$familyConf >= mccID1_2
confG1_2$FN <- confG1_2$trueTaxa == confG1_2$assignedTaxa & confG1_2$familyConf < mccID1_2

###########
# Genus Testing - 2nd run - F2_F3/F1
taxonomy <- paste("Root", gen_f2_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold2_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold1, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- gen_f1
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confG2_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confG2_3$outcome <- confG2_3$trueTaxa == confG2_3$assignedTaxa
mccID2_3 <- prediction(confG2_3$familyConf, confG2_3$outcome)
mccID2_3 <- performance(mccID2_3,"mat")
mccID2_3 <- data.frame(mccID2_3@y.values, mccID2_3@x.values)
colnames(mccID2_3) <- c("MCC Score", "Confidence Cutoff")
mccID2_3 <- mccID2_3[2:(nrow(mccID2_3)-1),]
optMCC <- which(mccID2_3$`MCC Score` == max(mccID2_3$`MCC Score`))
mccID2_3 <- mccID2_3$`Confidence Cutoff`[optMCC]
confG2_3$TP <- confG2_3$trueTaxa == confG2_3$assignedTaxa & confG2_3$familyConf >= mccID2_3
confG2_3$TN <- confG2_3$trueTaxa != confG2_3$assignedTaxa & confG2_3$familyConf < mccID2_3
confG2_3$FP <- confG2_3$trueTaxa != confG2_3$assignedTaxa & confG2_3$familyConf >= mccID2_3
confG2_3$FN <- confG2_3$trueTaxa == confG2_3$assignedTaxa & confG2_3$familyConf < mccID2_3

###########
# Genus Testing - 3rd run - F1_3/F2
taxonomy <- paste("Root", gen_f1_3, sep="; ")

# Training with train dataset
trainingSet <- LearnTaxa(fold1_3, taxonomy)

# classify the test sequences
ids <- IdTaxa(fold2, trainingSet, strand="top", type = "extended")

confidenceValues <- sapply(ids, `[`, 2)
confDataframe = do.call("rbind", lapply(confidenceValues, "[", 2))

confDataframe <- as.data.frame(confDataframe)
names(confDataframe)[1] <- paste("familyConf")
confDataframe$trueTaxa <- gen_f2
confDataframe$assignedTaxa <- do.call("rbind", lapply(sapply(ids, `[`, 1), "[", 2))
confG1_3 <- confDataframe

# MCC determination for optimal values and TP/TN/FP/FN
confG1_3$outcome <- confG1_3$trueTaxa == confG1_3$assignedTaxa
mccID1_3 <- prediction(confG1_3$familyConf, confG1_3$outcome)
mccID1_3 <- performance(mccID1_3,"mat")
mccID1_3 <- data.frame(mccID1_3@y.values, mccID1_3@x.values)
colnames(mccID1_3) <- c("MCC Score", "Confidence Cutoff")
mccID1_3 <- mccID1_3[2:(nrow(mccID1_3)-1),]
optMCC <- which(mccID1_3$`MCC Score` == max(mccID1_3$`MCC Score`))
mccID1_3 <- mccID1_3$`Confidence Cutoff`[optMCC]
confG1_3$TP <- confG1_3$trueTaxa == confG1_3$assignedTaxa & confG1_3$familyConf >= mccID1_3
confG1_3$TN <- confG1_3$trueTaxa != confG1_3$assignedTaxa & confG1_3$familyConf < mccID1_3
confG1_3$FP <- confG1_3$trueTaxa != confG1_3$assignedTaxa & confG1_3$familyConf >= mccID1_3
confG1_3$FN <- confG1_3$trueTaxa == confG1_3$assignedTaxa & confG1_3$familyConf < mccID1_3

#############
# Accuracy metrics IDTAXA
avgIncNumFamAmph1_2 <- length(which(confG1_2$FP == TRUE)) + length(which(confG1_2$TN == TRUE))
incAssignFamAmph1_2 <- avgIncNumFamAmph1_2 / nrow(confG1_2)
incAssignFamAmph1_2 * 100

avgIncNumFamAmph2_3 <- length(which(confG2_3$FP == TRUE)) + length(which(confG2_3$TN == TRUE))
incAssignFamAmph2_3 <- avgIncNumFamAmph2_3 / nrow(confG2_3)
incAssignFamAmph2_3 * 100

avgIncNumFamAmph1_3 <- length(which(confG1_3$FP == TRUE)) + length(which(confG1_3$TN == TRUE))
incAssignFamAmph1_3 <- avgIncNumFamAmph1_3 / nrow(confG1_3)
incAssignFamAmph1_3 * 100

amphGen <- rbind(confG1_2, confG2_3, confG1_3)

avgIncNumGenAmph <- length(which(amphGen$FP == TRUE)) + length(which(amphGen$TN == TRUE))
incAssignGenAvgAmph <- avgIncNumGenAmph / nrow(amphGen)
incAssignGenAvgAmph * 100

#######
# BLAST Tests - Data Import and database creation
# Amphibia BLAST
writeXStringSet(fold1_2, "fold1_2ampG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3ampG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3ampG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2ampG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3ampG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3ampG.fas"))

#######
# Mammalia BLAST

writeXStringSet(fold1_2, "fold1_2mamG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3mamG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3mamG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2mamG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2mamG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3mamG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3mamG.fas"))

########
# Gastropda BLAST
writeXStringSet(fold1_2, "fold1_2gasG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3gasG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3gasG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2gasG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3gasG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3gasG.fas"))

########
# Araneae BLAST
writeXStringSet(fold1_2, "fold1_2aranG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3aranG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3aranG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2aranG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3aranG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3aranG.fas"))

#########
# Diptera BLAST
writeXStringSet(fold1_2, "fold1_2dipG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3dipG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3dipG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipG.fas"))

#########
# Diptera BLAST
writeXStringSet(fold1_2, "fold1_2anthoG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3anthoG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3anthoG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2anthoG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2anthoG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3anthoG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3anthoG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3anthoG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3anthoG.fas"))

########
# Fish BLAST
writeXStringSet(fold1_2, "fold1_2fishG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3fishG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3fishG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2fishG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3fishG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3fishG.fas"))

########
# Hymenoptera BLAST
writeXStringSet(fold1_2, "fold1_2hymG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3hymG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3hymG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2hymG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3hymG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3hymG.fas"))

########
# Lepidoptera BLAST
writeXStringSet(fold1_2, "fold1_2lepG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3lepG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3lepG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2lepG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3lepG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3lepG.fas"))

########
# Aves BLAST
writeXStringSet(fold1_2, "fold1_2avesG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3avesG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3avesG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2avesG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3avesG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3avesG.fas"))

########
# Dip BLAST
writeXStringSet(fold1_2, "fold1_2dipG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3dipG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3dipG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2dipG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3dipG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3dipG.fas"))

########
# Bee BLAST
writeXStringSet(fold1_2, "fold1_2beeG.fas", format = "FASTA")
writeXStringSet(fold1_3, "fold1_3beeG.fas", format = "FASTA")
writeXStringSet(fold2_3, "fold2_3beeG.fas", format = "FASTA")

# Make a custom BLAST database - 1st fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2beeG.fas"), dbtype = "nucl")
db1_2 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_2beeG.fas"))

# 2nd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3beeG.fas"), dbtype = "nucl")
db1_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold1_3beeG.fas"))

# 3rd fold
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3beeG.fas"), dbtype = "nucl")
db2_3 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "fold2_3beeG.fas"))

##########
# BLAST predictions

# Genus predictions - all runs
predictF2F3 <- predict(db2_3, fold1)

predictF1F3 <- predict(db1_3, fold2)

predictF1F2 <- predict(db1_2, fold3)

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

# For bees only
mccb1_2 <- mccID1_2 
mccb1_3 <- mccID1_3 
mccb2_3 <- mccID2_3 

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST1 <- data.table(predictF1F2, key="QueryID")
dfBLAST1 <- dfBLAST1[, head(.SD, 1), by=QueryID]
blastG1_2 <- dfBLAST1

# MCC determination for optimal values and TP/TN/FP/FN
blastG1_2$outcome <- blastG1_2$trueTaxa == blastG1_2$assignedTaxa
mccb1_2 <- prediction(blastG1_2$Perc.Ident, blastG1_2$outcome)
mccb1_2 <- performance(mccb1_2,"mat")
mccb1_2 <- data.frame(mccb1_2@y.values, mccb1_2@x.values)
colnames(mccb1_2) <- c("MCC Score", "Confidence Cutoff")
mccb1_2 <- mccb1_2[2:(nrow(mccb1_2)-1),]
optMCCb <- which(mccb1_2$`MCC Score` == max(mccb1_2$`MCC Score`))
mccb1_2 <- mccb1_2$`Confidence Cutoff`[optMCCb]
blastG1_2$TP <- blastG1_2$trueTaxa == blastG1_2$assignedTaxa & blastG1_2$Perc.Ident >= mccb1_2
blastG1_2$TN <- blastG1_2$trueTaxa != blastG1_2$assignedTaxa & blastG1_2$Perc.Ident < mccb1_2
blastG1_2$FP <- blastG1_2$trueTaxa != blastG1_2$assignedTaxa & blastG1_2$Perc.Ident >= mccb1_2
blastG1_2$FN <- blastG1_2$trueTaxa == blastG1_2$assignedTaxa & blastG1_2$Perc.Ident < mccb1_2

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST2 <- data.table(predictF1F3, key="QueryID")
dfBLAST2 <- dfBLAST2[, head(.SD, 1), by=QueryID]
blastG1_3 <- dfBLAST2

# MCC determination for optimal values and TP/TN/FP/FN
blastG1_3$outcome <- blastG1_3$trueTaxa == blastG1_3$assignedTaxa
mccb1_3 <- prediction(blastG1_3$Perc.Ident, blastG1_3$outcome)
mccb1_3 <- performance(mccb1_3,"mat")
mccb1_3 <- data.frame(mccb1_3@y.values, mccb1_3@x.values)
colnames(mccb1_3) <- c("MCC Score", "Confidence Cutoff")
mccb1_3 <- mccb1_3[2:(nrow(mccb1_3)-1),]
optMCCb <- which(mccb1_3$`MCC Score` == max(mccb1_3$`MCC Score`))
mccb1_3 <- mccb1_3$`Confidence Cutoff`[optMCCb]
blastG1_3$TP <- blastG1_3$trueTaxa == blastG1_3$assignedTaxa & blastG1_3$Perc.Ident >= mccb1_3
blastG1_3$TN <- blastG1_3$trueTaxa != blastG1_3$assignedTaxa & blastG1_3$Perc.Ident < mccb1_3
blastG1_3$FP <- blastG1_3$trueTaxa != blastG1_3$assignedTaxa & blastG1_3$Perc.Ident >= mccb1_3
blastG1_3$FN <- blastG1_3$trueTaxa == blastG1_3$assignedTaxa & blastG1_3$Perc.Ident < mccb1_3

# For each set of queries, take the top match (matches are already sorted according to percent identity, alignment length, ecore etc.)
dfBLAST3 <- data.table(predictF2F3, key="QueryID")
dfBLAST3 <- dfBLAST3[, head(.SD, 1), by=QueryID]
blastG2_3 <- dfBLAST3

# MCC determination for optimal values and TP/TN/FP/FN
blastG2_3$outcome <- blastG2_3$trueTaxa == blastG2_3$assignedTaxa
mccb2_3 <- prediction(blastG2_3$Perc.Ident, blastG2_3$outcome)
mccb2_3 <- performance(mccb2_3,"mat")
mccb2_3 <- data.frame(mccb2_3@y.values, mccb2_3@x.values)
colnames(mccb2_3) <- c("MCC Score", "Confidence Cutoff")
mccb2_3 <- mccb2_3[2:(nrow(mccb2_3)-1),]
optMCCb <- which(mccb2_3$`MCC Score` == max(mccb2_3$`MCC Score`))
mccb2_3 <- mccb2_3$`Confidence Cutoff`[optMCCb]
blastG2_3$TP <- blastG2_3$trueTaxa == blastG2_3$assignedTaxa & blastG2_3$Perc.Ident >= mccb2_3
blastG2_3$TN <- blastG2_3$trueTaxa != blastG2_3$assignedTaxa & blastG2_3$Perc.Ident < mccb2_3
blastG2_3$FP <- blastG2_3$trueTaxa != blastG2_3$assignedTaxa & blastG2_3$Perc.Ident >= mccb2_3
blastG2_3$FN <- blastG2_3$trueTaxa == blastG2_3$assignedTaxa & blastG2_3$Perc.Ident < mccb2_3

rm(predictF1F2)
rm(predictF1F3)
rm(predictF2F3)
rm(familyPredictF1F2)
rm(familyPredictF1F3)
rm(familyPredictF2F3)
rm(familyTrue1)
rm(familyTrue2)
rm(familyTrue3)

# Accuracy metrics BLAST
avgIncNumFamAmph1_2b <- length(which(blastG1_2$FP == TRUE)) + length(which(blastG1_2$TN == TRUE))
incAssignFamAmph1_2b <- avgIncNumFamAmph1_2b / nrow(blastG1_2)
incAssignFamAmph1_2b * 100

avgIncNumFamAmph2_3b <- length(which(blastG2_3$FP == TRUE)) + length(which(blastG2_3$TN == TRUE))
incAssignFamAmph2_3b <- avgIncNumFamAmph2_3b / nrow(blastG2_3)
incAssignFamAmph2_3b * 100

avgIncNumFamAmph1_3b <- length(which(blastG1_3$FP == TRUE)) + length(which(blastG1_3$TN == TRUE))
incAssignFamAmph1_3b <- avgIncNumFamAmph1_3b / nrow(blastG1_3)
incAssignFamAmph1_3b * 100

amphGenB <- rbind(blastG1_2, blastG1_3, blastG2_3)

avgIncNumGenAmphB <- length(which(amphGenB$FP == TRUE)) + length(which(amphGenB$TN == TRUE))
incAssignGenAvgAmphB <- avgIncNumGenAmphB / nrow(amphGenB)
incAssignGenAvgAmphB * 100

#########
# ROC Curve Setup

# Genus - Amph
amphGenB$outcome <- as.integer(as.logical(amphGenB$trueTaxa == amphGenB$assignedTaxa))

####################
# RDP Classifier 

# First write out all files to fasta
writeXStringSet(fold1, "fold1ampG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2ampG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3ampG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1mamG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2mamG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3mamG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1gasG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2gasG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3gasG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1aranG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2aranG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3aranG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1dipG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2dipG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3dipG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1anthoG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2anthoG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3anthoG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1fishG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2fishG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3fishG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1hymG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2hymG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3hymG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1lepG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2lepG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3lepG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1avesG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2avesG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3avesG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1dipG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2dipG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3dipG.fas", format = "FASTA", compress = TRUE)

writeXStringSet(fold1, "fold1beeG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "fold2beeG.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "fold3beeG.fas", format = "FASTA", compress = TRUE)

# Reformat prediction data
# Converting to old variable name for ease of pipeline running
predictRDPf1_2ampG <- predictRDPf1_2beeG
predictRDPf1_3ampG <- predictRDPf1_3beeG
predictRDPf2_3ampG <- predictRDPf2_3beeG
customRDPf1_2ampG <- customRDPf1_2beeG
customRDPf1_3ampG <- customRDPf1_3beeG
customRDPf2_3ampG <- customRDPf2_3beeG

# Genus
predictRDPf1_2ampG$Genus[is.na(predictRDPf1_2ampG$Genus)] <- "Unknown"
predictRDPf1_3ampG$Genus[is.na(predictRDPf1_3ampG$Genus)] <- "Unknown"
predictRDPf2_3ampG$Genus[is.na(predictRDPf2_3ampG$Genus)] <- "Unknown"

predictRDPf1_2ampG$trueTaxaG <- gen_f3
predictRDPf1_3ampG$trueTaxaG <- gen_f2
predictRDPf2_3ampG$trueTaxaG <- gen_f1

predictRDPf1_2ampG$confG <- customRDPf1_2ampG$Genus * 100  
predictRDPf1_3ampG$confG <- customRDPf1_3ampG$Genus * 100
predictRDPf2_3ampG$confG <- customRDPf2_3ampG$Genus * 100

# For bees only
mccr1_2 <- mccID1_2 
mccr1_3 <- mccID1_3 
mccr2_3 <- mccID2_3 

# 1st run
# TP Values 
RDPF1_2G <- predictRDPf1_2ampG

# MCC determination for optimal values and TP/TN/FP/FN
RDPF1_2G$outcome <- RDPF1_2G$trueTaxaG == RDPF1_2G$Genus
mccr1_2 <- prediction(RDPF1_2G$confG, RDPF1_2G$outcome)
mccr1_2 <- performance(mccr1_2,"mat")
mccr1_2 <- data.frame(mccr1_2@y.values, mccr1_2@x.values)
colnames(mccr1_2) <- c("MCC Score", "Confidence Cutoff")
mccr1_2 <- mccr1_2[2:(nrow(mccr1_2)-1),]
optMCCr <- which(mccr1_2$`MCC Score` == max(mccr1_2$`MCC Score`))
mccr1_2 <- mccr1_2$`Confidence Cutoff`[optMCCr]
RDPF1_2G$TP <- RDPF1_2G$trueTaxaG == RDPF1_2G$Genus & RDPF1_2G$confG >= mccr1_2
RDPF1_2G$TN <- RDPF1_2G$trueTaxaG != RDPF1_2G$Genus & RDPF1_2G$confG < mccr1_2
RDPF1_2G$FP <- RDPF1_2G$trueTaxaG != RDPF1_2G$Genus & RDPF1_2G$confG >= mccr1_2
RDPF1_2G$FN <- RDPF1_2G$trueTaxaG == RDPF1_2G$Genus & RDPF1_2G$confG < mccr1_2

# 2nd run
RDPF1_3G <- predictRDPf1_3ampG

RDPF1_3G$outcome <- RDPF1_3G$trueTaxaG == RDPF1_3G$Genus
mccr1_3 <- prediction(RDPF1_3G$confG, RDPF1_3G$outcome)
mccr1_3 <- performance(mccr1_3,"mat")
mccr1_3 <- data.frame(mccr1_3@y.values, mccr1_3@x.values)
colnames(mccr1_3) <- c("MCC Score", "Confidence Cutoff")
mccr1_3 <- mccr1_3[2:(nrow(mccr1_3)-1),]
optMCCr <- which(mccr1_3$`MCC Score` == max(mccr1_3$`MCC Score`))
mccr1_3 <- mccr1_3$`Confidence Cutoff`[optMCCr]
RDPF1_3G$TP <- RDPF1_3G$trueTaxaG == RDPF1_3G$Genus & RDPF1_3G$confG >= mccr1_3
RDPF1_3G$TN <- RDPF1_3G$trueTaxaG != RDPF1_3G$Genus & RDPF1_3G$confG < mccr1_3
RDPF1_3G$FP <- RDPF1_3G$trueTaxaG != RDPF1_3G$Genus & RDPF1_3G$confG >= mccr1_3
RDPF1_3G$FN <- RDPF1_3G$trueTaxaG == RDPF1_3G$Genus & RDPF1_3G$confG < mccr1_3

# 3rd run
# TP Values 
RDPF2_3G <- predictRDPf2_3ampG

# MCC determination for optimal values and TP/TN/FP/FN
RDPF2_3G$outcome <- RDPF2_3G$trueTaxaG == RDPF2_3G$Genus
mccr2_3 <- prediction(RDPF2_3G$confG, RDPF2_3G$outcome)
mccr2_3 <- performance(mccr2_3,"mat")
mccr2_3 <- data.frame(mccr2_3@y.values, mccr2_3@x.values)
colnames(mccr2_3) <- c("MCC Score", "Confidence Cutoff")
mccr2_3 <- mccr2_3[2:(nrow(mccr2_3)-1),]
optMCCr <- which(mccr2_3$`MCC Score` == max(mccr2_3$`MCC Score`))
mccr2_3 <- mccr2_3$`Confidence Cutoff`[optMCCr]
RDPF2_3G$TP <- RDPF2_3G$trueTaxaG == RDPF2_3G$Genus & RDPF2_3G$confG >= mccr2_3
RDPF2_3G$TN <- RDPF2_3G$trueTaxaG != RDPF2_3G$Genus & RDPF2_3G$confG < mccr2_3
RDPF2_3G$FP <- RDPF2_3G$trueTaxaG != RDPF2_3G$Genus & RDPF2_3G$confG >= mccr2_3
RDPF2_3G$FN <- RDPF2_3G$trueTaxaG == RDPF2_3G$Genus & RDPF2_3G$confG < mccr2_3

#############
# Accuracy metrics RDP
avgIncNumFamAmph1_2r <- length(which(RDPF1_2G$FP == TRUE)) + length(which(RDPF1_2G$TN == TRUE))
incAssignFamAmph1_2r <- avgIncNumFamAmph1_2r / nrow(RDPF1_2G)
incAssignFamAmph1_2r * 100

avgIncNumFamAmph2_3r <- length(which(RDPF2_3G$FP == TRUE)) + length(which(RDPF2_3G$TN == TRUE))
incAssignFamAmph2_3r <- avgIncNumFamAmph2_3r / nrow(RDPF2_3G)
incAssignFamAmph2_3r * 100

avgIncNumFamAmph1_3r <- length(which(RDPF1_3G$FP == TRUE)) + length(which(RDPF1_3G$TN == TRUE))
incAssignFamAmph1_3r <- avgIncNumFamAmph1_3r / nrow(RDPF1_3G)
incAssignFamAmph1_3r * 100

amphGenRDP <- rbind(RDPF1_2G, RDPF1_3G, RDPF2_3G)

avgIncNumGenAmphR <- length(which(amphGenRDP$FP == TRUE)) + length(which(amphGenRDP$TN == TRUE))
incAssignGenAvgAmphR <- avgIncNumGenAmphR / nrow(amphGenRDP)
incAssignGenAvgAmphR * 100

#########
# ROC Curve Setup

# Genus - Amph
amphGenRDP$outcomeG <- as.integer(as.logical(amphGenRDP$trueTaxaG == amphGenRDP$Genus))
amphGen$outcomeG <- as.integer(as.logical(amphGen$trueTaxa == amphGen$assignedTaxa))
