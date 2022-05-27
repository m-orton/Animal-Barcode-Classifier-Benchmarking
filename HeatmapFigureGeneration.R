# ggplot approach
library(wesanderson)
library("ggplot2")
library(gridExtra)
library(grid)

##########
# Average % Incorrect Identifications
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal <- rev(pal)

colnames(ResultsHeatmapValues)[2] <- c("Avg % Incorrect Identifications")
topBlastHitEdit <- which(ResultsHeatmapValues$`Software Tool` == "BLAST")
ResultsHeatmapValues$`Software Tool`[topBlastHitEdit] <- "TOPBLASTHIT"
ResultsHeatmapValues$`fTaxonomic Group`= factor(ResultsHeatmapValues$`Taxonomic Group`, 
                                               levels=c('Cross-Taxon Average',
                                                        'Actinopterygii','Amphibia','Anthophila','Araneae',
                                                        'Aves','Diptera','Gastropoda','Hymenoptera',
                                                        'Lepidoptera','Mammalia'))

mine.heatmap1 <- ggplot(data = ResultsHeatmapValues, mapping = aes(x = `Taxonomic Rank`,
                                                                   y = `Software Tool`,
                                                                   fill = `Avg % Incorrect Identifications`)) +
  scale_fill_gradientn(colours=pal,
                       trans = "log10", 
                       limits = c(0.2,64)) + 
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        plot.margin=unit(c(-0.2,1,1,1),"cm")) +
  facet_wrap(~ `fTaxonomic Group`, nrow = 3) +
  coord_equal() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))

mine.heatmap1

##########
# AUC Heatmap
pal <- wes_palette("Zissou1", 15, type = "continuous")
pal <- rev(pal)

topBlastHitEdit <- which(AUCHeatmapValuesSep29$`Software Tool` == "BLAST")
AUCHeatmapValuesSep29$`Software Tool`[topBlastHitEdit] <- "TOPBLASTHIT"
AUCHeatmapValuesSep29$`fTaxonomic Group`= factor(AUCHeatmapValuesSep29$`Taxonomic Group`, 
                                                levels=c('Cross-Taxon Average',
                                                         'Actinopterygii','Amphibia','Anthophila','Araneae',
                                                         'Aves','Diptera','Gastropoda','Hymenoptera',
                                                         'Lepidoptera','Mammalia'))

mine.heatmap2 <- ggplot(data = AUCHeatmapValuesSep29, mapping = aes(x = `Taxonomic Rank`,
                                                                    y = `Software Tool`,
                                                                    fill = `AUC (TP/FP)`)) +
  scale_fill_gradientn(colours=pal) + 
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        plot.margin=unit(c(-0.2,1,1,1),"cm")) +
  facet_wrap(~ `fTaxonomic Group`, nrow = 3) +
  coord_equal() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))
