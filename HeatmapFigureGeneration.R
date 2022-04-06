# ggplot approach
library(wesanderson)
library("ggplot2")
library(gridExtra)
library(grid)

pal <- wes_palette("Zissou1", 10, type = "continuous")
pal <- rev(pal)
# pal <- brewer.pal(n = 11, name = 'Spectral')
# pal <- rev(pal)
colnames(ResultsHeatmapValues)[2] <- c("Avg % Incorrect Classifications")

# First 3 Groups Ordered
group1 <- c("Cross-Taxon Average", "Actinopterygii", "Amphibia", "Anthophila", "Araneae", "Aves")
group1 <- which(ResultsHeatmapValues$`Taxonomic Group` %in% group1)
group1 <- ResultsHeatmapValues[group1,]

# fill = `% Incorrect Classifications`
mine.heatmap1 <- ggplot(data = group1, mapping = aes(x = `Taxonomic Rank`,
                                                     y = `Software Tool`,
                                                     fill = `Avg % Incorrect Classifications`)) +
  scale_fill_gradientn(colours=pal,
                       trans = "log10", # can remove trans and limits, this log transforms the color gradient used
                       limits = c(0.2,60)) + #, 
                       # breaks = c(1,2,4,8,16,40,60)) +
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        plot.margin=unit(c(1,1,-0.2,1),"cm")) + 
  facet_grid(~ `Taxonomic Group`) +
  coord_equal() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))

# , scales = "free_x", space = "free_x"
mine.heatmap1

# Next 3
group2 <- c("Diptera", "Gastropoda", "Hymenoptera", "Lepidoptera", "Mammalia")
group2 <- which(ResultsHeatmapValuesSep29$`Taxonomic Group` %in% group2)
group2 <- ResultsHeatmapValuesSep29[group2,]

# fill = `% Incorrect Classifications`
mine.heatmap2 <- ggplot(data = group2, mapping = aes(x = `Taxonomic Rank`,
                                                     y = `Software Tool`,
                                                     fill = `Avg % Incorrect Classifications`)) +
  scale_fill_gradientn(colours=pal,
                       trans = "log10", 
                       limits = c(0.2,64)) + #, 
                       # breaks = c(0.5,1,2,4,8,16,32,64)) +
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(-0.2,1,1,1),"cm")) + 
  facet_grid(~ `Taxonomic Group`) +
  coord_equal()

# , scales = "free_x", space = "free_x"
mine.heatmap2

completeheat <- grid.arrange(mine.heatmap1, mine.heatmap2)

##########
# AUC Heatmap
ResultsHeatmapValuesSep29 <- AUCHeatmapValuesSep29

ct <- c("Cross-Taxon Average")
ct <- which(ResultsHeatmapValuesSep29$`Taxonomic Group` %in% ct)
ct <- ResultsHeatmapValuesSep29[ct,]

mine.heatmapct <- ggplot(data = ct, mapping = aes(x = `Taxonomic Rank`,
                                                  y = `Software Tool`,
                                                  fill = `AUC (TP/FP)`)) +
  scale_fill_gradientn(colours=pal) + #, 
  # breaks = c(1,2,4,8,16,40,60)) +
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=18),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        plot.margin=unit(c(1,1,-0.2,1),"cm")) + 
  facet_grid(~ `Taxonomic Group`) +
  coord_equal() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))

# , scales = "free_x", space = "free_x"
mine.heatmapct

# First 3 Groups Ordered
group1 <- c("Cross-Taxon Average", "Actinopterygii", "Amphibia", "Anthophila", "Araneae", "Aves")
group1 <- which(ResultsHeatmapValuesSep29$`Taxonomic Group` %in% group1)
group1 <- ResultsHeatmapValuesSep29[group1,]

# fill = `% Incorrect Classifications`
mine.heatmap1 <- ggplot(data = group1, mapping = aes(x = `Taxonomic Rank`,
                                                     y = `Software Tool`,
                                                     fill = `AUC (TP/FP)`)) +
  scale_fill_gradientn(colours=pal) + #, 
  # breaks = c(1,2,4,8,16,40,60)) +
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=16),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.direction = "horizontal",
        legend.position="top",
        plot.margin=unit(c(1,1,-0.2,1),"cm")) + 
  facet_grid(~ `Taxonomic Group`) +
  coord_equal() +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))

# , scales = "free_x", space = "free_x"
mine.heatmap1

# Next 3
group2 <- c("Diptera", "Gastropoda", "Hymenoptera", "Lepidoptera", "Mammalia")
group2 <- which(ResultsHeatmapValuesSep29$`Taxonomic Group` %in% group2)
group2 <- ResultsHeatmapValuesSep29[group2,]

# fill = `% Incorrect Classifications`
mine.heatmap2 <- ggplot(data = group2, mapping = aes(x = `Taxonomic Rank`,
                                                     y = `Software Tool`,
                                                     fill = `AUC (TP/FP)`)) +
  scale_fill_gradientn(colours=pal) + #, 
  # breaks = c(0.5,1,2,4,8,16,32,64)) +
  geom_tile(color="white", size=0.1) +
  theme_dark() +
  theme(text=element_text(size=18),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(-0.2,1,1,1),"cm")) + 
  facet_grid(~ `Taxonomic Group`) +
  coord_equal()

# , scales = "free_x", space = "free_x"
mine.heatmap2

completeheat <- grid.arrange(mine.heatmapct, mine.heatmap1, mine.heatmap2)

