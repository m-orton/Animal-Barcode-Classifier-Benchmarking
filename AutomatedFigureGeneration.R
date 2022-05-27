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
library(ggrepel)
# install.packages("wesanderson")
library(wesanderson)
# install.packages("ggsci")
library("ggsci")
library(pROC)
library(gridExtra)
library(grid)
# install.packages("lemon")
library(lemon)

#########
lines <- c("Family" = "solid", "Genus" = "dotted", "Species" = "dashed")
pal <- wes_palette("Zissou1", 100, type = "continuous")
# linesbees <- c("Genus" = "solid", "Species" = "dotted")

# Plotting of ROC curve IDTAXA
amphFam$outcomeG <- amphGen$outcomeG
amphFam$confG <- amphGen$familyConf
amphFam$outcomeS <- amphSp$outcome
amphFam$confS <- amphSp$familyConf
amphFam$conf <- amphFam$familyConf

plotROClistID <- list(Family=roc(amphFam$outcome, amphFam$conf, ci=TRUE, auc=TRUE), 
                     Genus=roc(amphFam$outcomeG, amphFam$confG, ci=TRUE, auc=TRUE), 
                     Species=roc(amphFam$outcomeS, amphFam$confS, ci=TRUE, auc=TRUE))
#########
# For bees
# amphFam <- amphGen
# amphFamB <- amphGenB
# amphRDP <- amphGenRDP

# amphFam$outcomeG <- amphGen$outcomeG
# amphFam$confG <- amphGen$familyConf
# amphFam$outcomeS <- amphSp$outcome
# amphFam$confS <- amphSp$familyConf
# amphFam$conf <- amphFam$familyConf

# plotROClistID <- list(Genus=roc(amphFam$outcomeG, amphFam$confG, ci=TRUE, auc=TRUE), 
#                       Species=roc(amphFam$outcomeS, amphFam$confS, ci=TRUE, auc=TRUE))

########
# Aves & Diptera

amphFam$conf <- amphFam$familyConf
plotROClistID <- list(Family=roc(amphFam$outcome, amphFam$conf, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for family
coordsF <- coords(plotROClistID[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsF$threshold[1] <- 0
coordsF$threshold[nrow(coordsF)] <- 100
coordsF = coordsF[order(coordsF[,'tpr'],-coordsF[,'fpr']),]
coordsF = coordsF[!duplicated(coordsF$tpr),]

# Grouping and values for legend
coordsF$group <- as.character(paste("Family"," (AUC=", round(as.numeric(plotROClistID$Family$auc), 2), ")", sep=""))

########
# Aves & Diptera
plotROClistID <- list(Genus=roc(amphGen$outcome, amphGen$familyConf, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for genus
coordsG <- coords(plotROClistID[[2]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsG$threshold[1] <- 0
coordsG$threshold[nrow(coordsG)] <- 100
coordsG = coordsG[order(coordsG[,'tpr'],-coordsG[,'fpr']),]
coordsG = coordsG[!duplicated(coordsG$tpr),]

coordsG$group <- as.character(paste("Genus"," (AUC=", round(as.numeric(plotROClistID$Genus$auc), 2), ")", sep=""))

########
# Aves & Diptera
# plotROClistID <- list(Species=roc(amphSp$outcome, amphSp$familyConf, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for species
coordsS <- coords(plotROClistID[[3]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsS$threshold[1] <- 0
coordsS$threshold[nrow(coordsS)] <- 100
coordsS = coordsS[order(coordsS[,'tpr'],-coordsS[,'fpr']),]
coordsS = coordsS[!duplicated(coordsS$tpr),]

coordsS$group <- as.character(paste("Species"," (AUC=", round(as.numeric(plotROClistID$Species$auc), 2), ")", sep=""))

coordsAllID <- rbind(coordsF, coordsG, coordsS) 

# Bees
# coordsAllID <- rbind(coordsG, coordsS) 

# Using plotROC
# ROC curves for IDTAXA
gpROC_ID <- ggplot() + 
  geom_line(data = coordsAllID, aes(`fpr`, `tpr`,  col=`threshold`), size = 1.3, alpha = 1) +
  scale_color_gradientn(name = paste("Confidence", "Threshold", sep="\n"), colours = pal) +
  xlab("False Positive Rate") + ylab("True Positive Rate")  +
  ggtitle("IDTAXA") +
  theme(text = element_text(size=16),
        legend.position = "none") +
  facet_grid(~ group, scales = "free_x", space = "free_x") 
gpROC_ID <- gpROC_ID + annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + theme(panel.spacing = unit(1, "lines"))

##############
# MCC Curves - IDTAXA

# Prediction dataframes
# Generate a prediction variable using percent identity as the cutoff, correct/incorrect as the outcome (1/0)
predIDTAX <- prediction(amphFam$familyConf, amphFam$outcome)
predIDTAXg <- prediction(amphGen$familyConf, amphGen$outcome)
# predIDTAXg <- prediction(amphGen$familyConf, amphGen$outcomeG)
predIDTAXs <- prediction(amphSp$familyConf, amphSp$outcome)

# Family
perfMCCIDf <- performance(predIDTAX,"mat")
perfMCCIDf <- data.frame(perfMCCIDf@y.values, perfMCCIDf@x.values)
colnames(perfMCCIDf) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDf$`Taxonomic Rank` <- "Family"

# Genus
perfMCCIDg <- performance(predIDTAXg,"mat")
perfMCCIDg <- data.frame(perfMCCIDg@y.values, perfMCCIDg@x.values)
colnames(perfMCCIDg) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDg$`Taxonomic Rank` <- "Genus"

# Species
perfMCCIDs <- performance(predIDTAXs,"mat")
perfMCCIDs <- data.frame(perfMCCIDs@y.values, perfMCCIDs@x.values)
colnames(perfMCCIDs) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDs$`Taxonomic Rank` <- "Species"

perfMCCID <- rbind(perfMCCIDf, perfMCCIDg, perfMCCIDs)

#########
# For bees
# perfMCCID <- rbind(perfMCCIDg, perfMCCIDs)

#########
# MCC curves for IDTAXA
idtaxaPlotMCC1 <- ggplot(data = perfMCCID, aes(`Confidence Cutoff`,`MCC Score`, col=`Taxonomic Rank`)) +
  geom_line(aes(linetype=`Taxonomic Rank`), size = 1.1, alpha = 1) +
  scale_linetype_manual(values = lines) +
  ggtitle("IDTAXA") +
  labs(x = "Confidence Threshold", 
       y = "Avg. MCC") +
  ylim(0, 1.01) +
  scale_color_locuszoom() +
  theme(text = element_text(size=16),
        legend.position = "none")

#########
# BLAST ROC
# Plotting of ROC curve IDTAXA
amphFamB$outcomeG <- amphGenB$outcome[1:nrow(amphFamB)]
amphFamB$confG <- amphGenB$Perc.Ident[1:nrow(amphFamB)]
amphFamB$outcomeS <- amphSpB$outcome[1:nrow(amphFamB)]
amphFamB$confS <- amphSpB$Perc.Ident[1:nrow(amphFamB)]
amphFamB$conf <- amphFamB$Perc.Ident[1:nrow(amphFamB)]

plotROClistB <- list(Family=roc(amphFamB$outcome, amphFamB$conf, ci=TRUE, auc=TRUE), 
                     Genus=roc(amphFamB$outcomeG, amphFamB$confG, ci=TRUE, auc=TRUE), 
                     Species=roc(amphFamB$outcomeS, amphFamB$confS, ci=TRUE, auc=TRUE))

# plotROClistB <- list(Family=roc(amphFamB$outcome, amphFamB$Perc.Ident, ci=TRUE, auc=TRUE))

# Bees
# plotROClistB <- list(Genus=roc(amphFamB$outcomeG, amphFamB$confG, ci=TRUE, auc=TRUE), 
#                      Species=roc(amphFamB$outcomeS, amphFamB$confS, ci=TRUE, auc=TRUE))

# Extract coordinate values for family
coordsF <- coords(plotROClistB[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsF$threshold[1] <- 0
coordsF$threshold[nrow(coordsF)] <- 100
coordsF = coordsF[order(coordsF[,'tpr'],-coordsF[,'fpr']),]
coordsF = coordsF[!duplicated(coordsF$tpr),]

# Grouping and values for legend
coordsF$group <- as.character(paste("Family"," (AUC=", round(as.numeric(plotROClistB$Family$auc), 2), ")", sep=""))

########
# Aves & Diptera
plotROClistB <- list(Genus=roc(amphGenB$outcome, amphGenB$Perc.Ident, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for genus
coordsG <- coords(plotROClistB[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsG$threshold[1] <- 0
coordsG$threshold[nrow(coordsG)] <- 100
coordsG = coordsG[order(coordsG[,'tpr'],-coordsG[,'fpr']),]
coordsG = coordsG[!duplicated(coordsG$tpr),]

coordsG$group <- as.character(paste("Genus"," (AUC=", round(as.numeric(plotROClistB$Genus$auc), 2), ")", sep=""))

########
# Aves & Diptera
plotROClistB <- list(Species=roc(amphGenB$outcome, amphGenB$Perc.Ident, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for species
coordsS <- coords(plotROClistB[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsS$threshold[1] <- 0
coordsS$threshold[nrow(coordsS)] <- 100
coordsS = coordsS[order(coordsS[,'tpr'],-coordsS[,'fpr']),]
coordsS = coordsS[!duplicated(coordsS$tpr),]

coordsS$group <- as.character(paste("Species"," (AUC=", round(as.numeric(plotROClistB$Species$auc), 2), ")", sep=""))

coordsAllB <- rbind(coordsF, coordsG, coordsS) 

# Bees
# coordsAllB <- rbind(coordsG, coordsS) 

# Plotting of ROC curve for BLAST 
gpROC_B <- ggplot() + 
  geom_line(data = coordsAllB, aes(`fpr`, `tpr`,  col=`threshold`), size = 1.3, alpha = 1) +
  scale_color_gradientn(name = paste("Confidence", "Threshold", sep="\n"), colours = pal) +
  xlab("False Positive Rate") +
  ggtitle("TOPBLASTHIT") +
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_grid(~ group, scales = "free_x", space = "free_x")

gpROC_B <- gpROC_B + annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + theme(panel.spacing = unit(1, "lines"))

##############
# MCC Curves - BLAST
# Generate a prediction variable using percent identity as the cutoff, correct/incorrect as the outcome (1/0)
predIDTAXB <- prediction(amphFamB$Perc.Ident, amphFamB$outcome)
predIDTAXgB <- prediction(amphGenB$Perc.Ident, amphGenB$outcome)
predIDTAXsB <- prediction(amphSpB$Perc.Ident, amphSpB$outcome)

# Family
perfMCCIDfB <- performance(predIDTAXB,"mat")
perfMCCIDfB <- data.frame(perfMCCIDfB@y.values, perfMCCIDfB@x.values)
colnames(perfMCCIDfB) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDfB$`Taxonomic Rank` <- "Family"

# Genus
perfMCCIDgB <- performance(predIDTAXgB,"mat")
perfMCCIDgB <- data.frame(perfMCCIDgB@y.values, perfMCCIDgB@x.values)
colnames(perfMCCIDgB) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDgB$`Taxonomic Rank` <- "Genus"

# Species
perfMCCIDsB <- performance(predIDTAXsB,"mat")
perfMCCIDsB <- data.frame(perfMCCIDsB@y.values, perfMCCIDsB@x.values)
colnames(perfMCCIDsB) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDsB$`Taxonomic Rank` <- "Species"

perfMCCIDB <- rbind(perfMCCIDfB, perfMCCIDgB, perfMCCIDsB)

########
# Bees 
# perfMCCIDB <- rbind(perfMCCIDgB, perfMCCIDsB)
########

# MCC curves for TOPBLASTHIT
blastPlotMCC1 <- ggplot(data = perfMCCIDB, aes(`Confidence Cutoff`,`MCC Score`, col=`Taxonomic Rank`)) +
  geom_line(aes(linetype=`Taxonomic Rank`), size = 1.1, alpha = 1) +
  scale_linetype_manual(values = lines) +
  ggtitle("TOPBLASTHIT") +
  labs(x = "Confidence Threshold") +
  ylim(0, 1.01) +
  scale_color_locuszoom() +
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "none")

#########
# RDP ROC
amphRDP <- amphFamRDP
amphRDP$outcomeG <- amphGenRDP$outcomeG
amphRDP$confG <- amphGenRDP$confG
amphRDP$outcomeS <- amphSpRDP$outcomeS
amphRDP$confS <- amphSpRDP$confS
amphRDP$conf <- amphRDP$confF

plotROClistR <- list(Family=roc(amphRDP$outcome, amphRDP$conf, ci=TRUE, auc=TRUE), 
                     Genus=roc(amphRDP$outcomeG, amphRDP$confG, ci=TRUE, auc=TRUE), 
                     Species=roc(amphRDP$outcomeS, amphRDP$confS, ci=TRUE, auc=TRUE))

# plotROClistR <- list(Family=roc(amphRDP$outcome, amphRDP$confF, ci=TRUE, auc=TRUE))

########
# For bees
# plotROClistR <- list(Genus=roc(amphRDP$outcomeG, amphRDP$confG, ci=TRUE, auc=TRUE), 
#                      Species=roc(amphRDP$outcomeS, amphRDP$confS, ci=TRUE, auc=TRUE))

# Extract coordinate values for family
coordsF <- coords(plotROClistR[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsF$threshold[1] <- 0
coordsF$threshold[nrow(coordsF)] <- 100
coordsF = coordsF[order(coordsF[,'tpr'],-coordsF[,'fpr']),]
coordsF = coordsF[!duplicated(coordsF$tpr),]

# Grouping and values for legend
coordsF$group <- as.character(paste("Family"," (AUC=", round(as.numeric(plotROClistR$Family$auc), 2), ")", sep=""))

########
# Aves & Diptera
plotROClistR <- list(Genus=roc(amphGen$outcome, amphGen$familyConf, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for genus
coordsG <- coords(plotROClistR[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsG$threshold[1] <- 0
coordsG$threshold[nrow(coordsG)] <- 100
coordsG = coordsG[order(coordsG[,'tpr'],-coordsG[,'fpr']),]
coordsG = coordsG[!duplicated(coordsG$tpr),]

coordsG$group <- as.character(paste("Genus"," (AUC=", round(as.numeric(plotROClistR$Genus$auc), 2), ")", sep=""))

########
# Aves & Diptera
plotROClistR <- list(Species=roc(amphGen$outcome, amphGen$familyConf, ci=TRUE, auc=TRUE))

#########

# Extract coordinate values for species
coordsS <- coords(plotROClistR[[1]], input = "threshold", ret=c("threshold","tpr", "fpr"), transpose = FALSE)
coordsS$threshold[1] <- 0
coordsS$threshold[nrow(coordsS)] <- 100
coordsS = coordsS[order(coordsS[,'tpr'],-coordsS[,'fpr']),]
coordsS = coordsS[!duplicated(coordsS$tpr),]

coordsS$group <- as.character(paste("Species"," (AUC=", round(as.numeric(plotROClistR$Species$auc), 2), ")", sep=""))

coordsAllR <- rbind(coordsF, coordsG, coordsS) 

#######
# For bees
# coordsAllR <- rbind(coordsG, coordsS) 

#######

# ROC curves for RDP
gpROC_RDP <- ggplot() + 
  geom_line(data = coordsAllR, aes(`fpr`, `tpr`,  col=`threshold`), size = 1.3, alpha = 1) +
  scale_color_gradientn(name = paste("Confidence", "Threshold", sep="\n"), colours = pal) +
  xlab("False Positive Rate") +
  ggtitle("RDP") +
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_grid(~ group, scales = "free_x", space = "free_x")

gpROC_RDP <- gpROC_RDP + annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + theme(panel.spacing = unit(1, "lines"))

##############
# MCC Curves - RDP
# Generate a prediction variable using percent identity as the cutoff, correct/incorrect as the outcome (1/0)
predIDTAXR <- prediction(amphRDP$confF, amphRDP$outcome)
predIDTAXgR <- prediction(amphGenRDP$confG, amphGenRDP$outcome)
predIDTAXsR <- prediction(amphSpRDP$confS, amphSpRDP$outcome)

########
# Aves and Diptera
predIDTAXgR <- prediction(amphGenRDP$confG, amphGenRDP$outcome)
predIDTAXsR <- prediction(amphSpRDP$confS, amphSpRDP$outcome)
#######

# Family
perfMCCIDfR <- performance(predIDTAXR,"mat")
perfMCCIDfR <- data.frame(perfMCCIDfR@y.values, perfMCCIDfR@x.values)
colnames(perfMCCIDfR) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDfR$`Taxonomic Rank` <- "Family"

# Genus
perfMCCIDgR <- performance(predIDTAXgR,"mat")
perfMCCIDgR <- data.frame(perfMCCIDgR@y.values, perfMCCIDgR@x.values)
colnames(perfMCCIDgR) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDgR$`Taxonomic Rank` <- "Genus"

# Species
perfMCCIDsR <- performance(predIDTAXsR,"mat")
perfMCCIDsR <- data.frame(perfMCCIDsR@y.values, perfMCCIDsR@x.values)
colnames(perfMCCIDsR) <- c("MCC Score", "Confidence Cutoff")
perfMCCIDsR$`Taxonomic Rank` <- "Species"

########
# For bees
# perfMCCIDR <- rbind(perfMCCIDgR, perfMCCIDsR)
########

perfMCCIDR <- rbind(perfMCCIDfR, perfMCCIDgR, perfMCCIDsR)

# MCC curves for RDP
RDPPlotMCC1 <- ggplot(data = perfMCCIDR, aes(`Confidence Cutoff`,`MCC Score`, col=`Taxonomic Rank`)) +
  geom_line(aes(linetype=`Taxonomic Rank`), size = 1.1, alpha = 1) +
  scale_linetype_manual(values = lines) +
  ggtitle("RDP") +
  labs(x = "Confidence Threshold") +
  ylim(0, 1.01) +
  scale_color_locuszoom() +
  theme(text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "none")

##########
taxGroup <- c("Actinopterygii", "Amphibia", "Anthophila", "Araneae", "Aves", "Diptera", "Gastropoda", "Hymenoptera", "Lepidoptera", "Mammalia")

# Aggregate ROC curves together
completeROC <- grid_arrange_shared_legend(gpROC_ID, gpROC_B, gpROC_RDP, ncol = 3, nrow = 1, 
                                          top=textGrob(taxGroup[10], gp = gpar(fontsize = 20)), position='bottom')

# Aggregate MCC curves together
completeMCC <- grid_arrange_shared_legend(idtaxaPlotMCC1, blastPlotMCC1, RDPPlotMCC1, ncol = 3, nrow = 1, position='bottom')

