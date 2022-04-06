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
library(ape)
# install.packages("caret")
library(caret)
# install.packages("ggnewscale")
library(ggnewscale)
# install.packages("RColorBrewer")
library(RColorBrewer)
# install.packages("cowplot")
library(profvis)
# install.packages("bench")
library(bench)
# install.packages("microbenchmark")
library(microbenchmark)
# install.packages("profmem")
library(profmem)

# Using mammalia for performance benchmarks
dna <- readDNAStringSet("mammaliaDNAUnfiltered.fas", format="FASTA")

dnaNames <- names(dna)

dnaNames <- strsplit(dnaNames, " ")

# For family level
familyID <- sapply(dnaNames, `[`, 2)
familyID <- strsplit(familyID, ";")
familyID <- sapply(familyID, `[`, 7)

dfLep <- data.frame(as.character(dna))
colnames(dfLep)[1] <- "dna"
dfLep$names <- familyID
dfLep$id <- 1:nrow(dfLep)

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

# Fold names
fam_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 1)
famids_f1 <- sapply(strsplit(names(fold1), ";"), `[`, 2)

fam_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 1)
famids_f2 <- sapply(strsplit(names(fold2), ";"), `[`, 2)

fam_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 1)
famids_f3 <- sapply(strsplit(names(fold3), ";"), `[`, 2)

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

##########
# Timing

# IDTAXA
sample <- sample(fold1_2, 100)
taxonomy <- paste("Root",names(sample), sep="; ")
trainingSet <- LearnTaxa(sample, taxonomy)

# Performance of IDTAXA
mbm100 = microbenchmark(
  ids <- IdTaxa(fold3, trainingSet), times=3
)
mbm100 <- round(mean(mbm100$time) / 1000000000, 2)

sample <- sample(fold1_2, 500)
taxonomy <- paste("Root",names(sample), sep="; ")
trainingSet <- LearnTaxa(sample, taxonomy)

mbm500 = microbenchmark(
  ids <- IdTaxa(fold3, trainingSet), times=3
)
mbm500 <- round(mean(mbm500$time) / 1000000000, 2)

sample <- sample(fold1_2, 1000)
taxonomy <- paste("Root",names(sample), sep="; ")
trainingSet <- LearnTaxa(sample, taxonomy)

mbm1000 = microbenchmark(
  ids <- IdTaxa(fold3, trainingSet), times=3
)
mbm1000 <- round(mean(mbm1000$time) / 1000000000, 2)

sample <- sample(fold1_2, 2500)
taxonomy <- paste("Root",names(sample), sep="; ")
trainingSet <- LearnTaxa(sample, taxonomy)

mbm2500 = microbenchmark(
  ids <- IdTaxa(fold3, trainingSet), times=3
)
mbm2500 <- round(mean(mbm2500$time) / 1000000000, 2)

sample <- sample(fold1_2, 5000)
taxonomy <- paste("Root",names(sample), sep="; ")
trainingSet <- LearnTaxa(sample, taxonomy)

mbm5000 = microbenchmark(
  ids <- IdTaxa(fold3, trainingSet), times=3
)
mbm5000 <- round(mean(mbm5000$time) / 1000000000, 2)

##########
# Performance of BLAST
sample <- sample(fold1_2, 100)
writeXStringSet(sample, "test100.fas", format = "FASTA")
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "test100.fas"), dbtype = "nucl")
db1_2_100 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "test100.fas"))

mbm100b = microbenchmark(
  predictF1F2 <- predict(db1_2_100, fold3), times=3
)
mbm100b <- round(mean(mbm100b$time) / 1000000000, 2)

sample <- sample(fold1_2, 500)
writeXStringSet(sample, "test500.fas", format = "FASTA")
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "test500.fas"), dbtype = "nucl")
db1_2_500 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "test500.fas"))

mbm500b = microbenchmark(
  predictF1F2 <- predict(db1_2_500, fold3), times=3
)
mbm500b <- round(mean(mbm500b$time) / 1000000000, 2)

sample <- sample(fold1_2, 1000)
writeXStringSet(sample, "test1000.fas", format = "FASTA")
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "test1000.fas"), dbtype = "nucl")
db1_2_1000 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "test1000.fas"))

mbm1000b = microbenchmark(
  predictF1F2 <- predict(db1_2_1000, fold3), times=3
)
mbm1000b <- round(mean(mbm1000b$time) / 1000000000, 2)

sample <- sample(fold1_2, 2500)
writeXStringSet(sample, "test2500.fas", format = "FASTA")
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "test2500.fas"), dbtype = "nucl")
db1_2_2500 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "test2500.fas"))

mbm2500b = microbenchmark(
  predictF1F2 <- predict(db1_2_2500, fold3), times=3
)
mbm2500b <- round(mean(mbm2500b$time) / 1000000000, 2)

sample <- sample(fold1_2, 5000)
writeXStringSet(sample, "test5000.fas", format = "FASTA")
makeblastdb(file.path("F:/Dropbox/BenchmarkingStudy", "test5000.fas"), dbtype = "nucl")
db1_2_5000 <- blast(file.path("F:/Dropbox/BenchmarkingStudy", "test5000.fas"))

mbm5000b = microbenchmark(
  predictF1F2 <- predict(db1_2_5000, fold3), times=3
)
mbm5000b <- round(mean(mbm5000b$time) / 1000000000, 2)

##########
# Performance of RDP
writeXStringSet(fold1, "testf1.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold2, "testf2.fas", format = "FASTA", compress = TRUE)
writeXStringSet(fold3, "testf3.fas", format = "FASTA", compress = TRUE)

##########
# Performance of RDP
RDPTrainNamesf1_2 <- paste(famids_f1_2, " ", "Root", ";", "Animalia", ";", "Chordata", ";", 
                           "Mammalia", ";", "NA", ";", "NA", ";", gen_f1_2, sep="")
rdpf1_2 <- fold1_2
names(rdpf1_2) <- RDPTrainNamesf1_2
sample <- sample(rdpf1_2, 100)
test100 <- trainRDP(sample, dir = "test100") 

mbm100r = microbenchmark(
  predict(test100, fold3, confidence=0), times=3
)
mbm100r <- round(mean(mbm100r$time) / 1000000000, 2)

sample <- sample(rdpf1_2, 500)
test500 <- trainRDP(sample, dir = "test500") 

mbm500r = microbenchmark(
  predict(test500, fold3, confidence=0), times=3 
)
mbm500r <- round(mean(mbm500r$time) / 1000000000, 2)

sample <- sample(rdpf1_2, 1000)
test1000 <- trainRDP(sample, dir = "test1000") 

mbm1000r = microbenchmark(
  predict(test1000, fold3, confidence=0), times=3
)
mbm1000r <- round(mean(mbm1000r$time) / 1000000000, 2)

sample <- sample(rdpf1_2, 2500)
test2500 <- trainRDP(sample, dir = "test2500") 

mbm2500r = microbenchmark(
  predict(test2500, fold3, confidence=0), times=3
)
mbm2500r <- round(mean(mbm2500r$time) / 1000000000, 2)

sample <- sample(rdpf1_2, 5000)
test5000 <- trainRDP(sample, dir = "test5000") 

mbm5000r = microbenchmark(
  predict(test5000, fold3, confidence=0), times=3
)
mbm5000r <- round(mean(mbm5000r$time) / 1000000000, 2)

##########
# Graphing 
pal <- wes_palette("Zissou1", 5, type = "discrete")
pal <- c(pal[1], pal[3], pal[5])

PerformanceValuesSep29$SE <- sd(PerformanceValuesSep29$`Average Processing Time (s)`) / sqrt(PerformanceValuesSep29$`Number of Sequences`)

perfGraph <- ggplot(data = PerformanceValuesSep29, aes(`Number of Sequences`, `Average Processing Time (s)`)) +
  geom_line(aes(col=`Software Tool`, linetype = `Software Tool`), size = 1, alpha = 1) +
  geom_point(aes(col=`Software Tool`, shape = `Software Tool`), size = 3) +
  scale_colour_manual(values=pal) +
  theme(text = element_text(size=18)) +
  ylim(0,2000)
  
  
