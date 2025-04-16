library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(plyr)
library(units)
library(ggpubr)
library(MASS)
library(performance)
library(bruceR)
library(compositions)
library(ggpubr)
library(compositions)
library(betareg)
library(glmnetUtils)
library(glmnet)
library(car)
library(ggforce)
library(betapart)
library(VennDiagram)
library(gridSVG)
library(iNEXT)
library(ggtext)

g <- bruceR::RGB(34, 174, 0)
br <- bruceR::RGB(139, 71, 37)
bl <- bruceR::RGB(99, 184, 255)

################################################################################
########################## Swarm distribution analysis #########################
################################################################################
otu.df <- read.csv("Urban_BGB_OTU.csv")

# Calculate total number of reads and otu richness per sample, excluding metadata columns
otu.df$no_reads <- rowSums(otu.df[, -which(colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type"))])
otu.df$otu_richness <- specnumber(otu.df[, -which(colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type", "no_reads"))])

soil.df <- otu.df[otu.df$Type == "Soil", ]
water.df <- otu.df[otu.df$Type == "Water", ]

# Filter out unwanted or unclassified phyla and remove rows with NA values
taxo.df <- read.csv("Urban_BGB_taxo.csv")
taxo.target.df <- taxo.df[!(taxo.df$Phylum %in% c("unclassified_Root",
                                                  "Oomycota_4762",
                                                  "Discosea_555280",
                                                  "phylum_class_order_Collodictyonidae_190322",
                                                  "Rhodophyta_2763") | 
                                 is.na(taxo.df$Phylum) == T), ]

otu.target.df <- otu.df[, which(colnames(otu.df) %in% taxo.target.df$otu_id)]
otu.df$sp_richness <- specnumber(otu.target.df)

# Calculate the number of reads assigned or unassigned to target OTUs overall, and based on medium type
otu.df$no_reads_assigned <- rowSums(otu.df[, which(colnames(otu.df) %in% colnames(otu.target.df))])
otu.df$no_reads_unassigned <- rowSums(otu.df[, -which(colnames(otu.df) %in% colnames(otu.target.df) | colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type"))])

otu.df$fraction_assigned_otu <- otu.df$sp_richness/otu.df$otu_richness
otu.df$fraction_assigned_reads <- otu.df$no_reads_assigned/otu.df$no_reads

soil.assigned.df <- data.frame(rowSums(otu.df[otu.df$Type == "Soil", -which(colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type"))]), 
                               otu.df[otu.df$Type == "Soil", ]$Site_No, 
                               otu.df[otu.df$Type == "Soil", ]$no_reads_assigned, 
                               otu.df[otu.df$Type == "Soil", ]$sp_richness,
                               otu.df[otu.df$Type == "Soil", ]$otu_richness, 
                               "Soil")
colnames(soil.assigned.df) <- c("no_reads", "Site_No", "no_reads_assigned", "sp_richness", "otu_richness", "Type")
water.assigned.df <- data.frame(rowSums(otu.df[otu.df$Type == "Water", -which(colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type"))]), 
                                otu.df[otu.df$Type == "Water", ]$Site_No, 
                                otu.df[otu.df$Type == "Water", ]$no_reads_assigned, 
                                otu.df[otu.df$Type == "Water", ]$sp_richness,
                                otu.df[otu.df$Type == "Water", ]$otu_richness, 
                                "Water")
colnames(water.assigned.df) <- c("no_reads", "Site_No", "no_reads_assigned", "sp_richness", "otu_richness", "Type")
otu.assigned.df <- rbind(soil.assigned.df, water.assigned.df)

# Figure S1
ggplot(otu.assigned.df, aes(x = otu_richness, y = sp_richness, color = Type)) + 
  geom_point(cex = 2) + 
  geom_smooth(method = "lm") + 
  scale_color_manual(values=c("sienna4", "steelblue1")) + 
  labs(x = "OTU richness", y = "Species richness") + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) + 
  guides(color = guide_legend(title="Medium type"))

# Figure S2
ggplot(otu.assigned.df, aes(x = Type, y = otu_richness, color = Type)) + 
  geom_violin(trim = F, alpha = 0, size = 1) + 
  ylab(label = "OTU richness") + 
  stat_summary(fun.y = base::mean, geom = "point", aes(color = Type), size = 4) + 
  stat_summary(fun.data = mean_sd, geom = "errorbar", aes(color = Type), width  = 0.03, size = 1) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Type), alpha = 0.4) + 
  scale_color_manual(values=c("sienna4", "steelblue1"))  +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) + 
  guides(color = guide_legend(title="Medium type"))

################################################################################
################################## Ordinations #################################
################################################################################
soil.df <- otu.df[otu.df$Type == "Soil", ]
water.df <- otu.df[otu.df$Type == "Water", ]

# Perform Non-Metric Multidimensional Scaling (NMDS) on OTU data, excluding metadata columns, and extract NMDS site scores
otu.nmds <- metaMDS(otu.df[, -c(which(colnames(otu.df) == "Site_No"): 
                               which(colnames(otu.df) == "fraction_assigned_reads"))],
                   distance = "jaccard",
                   k = 2,
                   trymax = 1000, 
                   wascores = TRUE)
otu.scrs <- as.data.frame(scores(otu.nmds, display = "sites")) 
otu.scrs$Type <- otu.df$Type

# Test for significant differences in community composition between soil and water using PERMANOVA (adonis2)
adonis2(otu.df[, -c(which(colnames(otu.df) == "Site_No"): 
                      which(colnames(otu.df) == "fraction_assigned_reads"))] ~ otu.df$Type, method = "jaccard")

# Perform beta-dispersion analysis (homogeneity of group dispersions) using Jaccard distance and create a dataframe to store distances from the beta-dispersion analysis
otu.disp <- betadisper(vegdist(otu.df[, -c(which(colnames(otu.df) == "Site_No"): 
                                             which(colnames(otu.df) == "fraction_assigned_reads"))], 
                               method = "jaccard", binary = T), 
                       otu.df$Type, 
                       type = "centroid")
otu.disp.df <- data.frame(otu.disp$distances, otu.disp$group, "Within type")
colnames(otu.disp.df) <- c("Distance", "Type", "Scale")
otu.disp.df$Site_No <- otu.df$Site_No

# Linear mixed model to test the effect of 'Type' on beta-dispersion distances (random effect: Site_No)
summary(lmer(Distance ~ Type + (1 | Site_No), data = otu.disp.df))

# Separate beta-dispersion analysis for soil and water samples
terr.disp <- betadisper(vegdist(soil.df[, -c(which(colnames(soil.df) == "Site_No"): 
                                              which(colnames(soil.df) == "fraction_assigned_reads"))], 
                                method = "jaccard", binary = T), 
                        soil.df$Site_No, type = "centroid")
aqua.disp <-betadisper(vegdist(water.df[, -c(which(colnames(water.df) == "Site_No"): 
                                               which(colnames(water.df) == "fraction_assigned_reads"))], 
                               method = "jaccard", binary = T), 
                       water.df$Site_No, type = "centroid")

# Combine soil and water beta-dispersion results into a single dataframe for modeling and visualization purposes
terr.disp.df <- data.frame(terr.disp$distances, "Soil", "Within sites")
aqua.disp.df <- data.frame(aqua.disp$distances, "Water", "Within sites")
colnames(terr.disp.df) <- c("Distance", "Type", "Scale")
colnames(aqua.disp.df) <- c("Distance", "Type", "Scale")

bgb.disp.df <- rbind(terr.disp.df, aqua.disp.df)
bgb.disp.df$Site_No <- c(soil.df$Site_No, water.df$Site_No)

summary(lmer(Distance ~ Type + (1 | Site_No), data = bgb.disp.df))

# Add beta-dispersion results to the main dataframe and create 'scaletype' factor for plotting
otu.disp.df <- rbind(otu.disp.df, bgb.disp.df)
otu.disp.df$scaletype <- paste0(otu.disp.df$Type, otu.disp.df$Scale)
otu.disp.df$scaletype <- factor(otu.disp.df$scaletype, levels = c("WaterWithin type", "WaterWithin sites", "SoilWithin type", "SoilWithin sites"), ordered = T)

summary(lmer(Distance ~ Scale + (1 | Site_No), data = otu.disp.df[otu.disp.df$Type == "Soil", ]))
summary(lmer(Distance ~ Scale + (1 | Site_No), data = otu.disp.df[otu.disp.df$Type == "Water", ]))

# Figure 2
plot1 <- ggplot(otu.scrs, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, color = Type), size = 3) +
  annotate(geom = "text", label = paste0("stress = ", signif(otu.nmds$stress, 5)), x=0.2, y=0.3, size = 6) + 
  scale_color_manual(values = c("sienna4", "steelblue1")) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) + 
  guides(color = guide_legend(title="Medium type"))

plot2 <- ggplot(otu.disp.df, aes(x = scaletype, y = Distance, color = Type)) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Type), alpha = 0.2) + 
  geom_violin(trim = F, alpha = 0, size = 0.6) + 
  stat_summary(fun.y = median, geom = "point", aes(color = Type), position = position_dodge(width = .9), size = 4) +
  stat_summary(fun.ymin = function(x) median(x) - sd(x), fun.ymax = function(x) min(median(x) + sd(x), 1), 
               geom = "errorbar", aes(color = Type), position = position_dodge(width = .9), width  = 0.05, size = 1.1) +
  geom_signif(test="t.test", 
              comparisons = list(c("SoilWithin type", "SoilWithin sites"), 
                                 c("WaterWithin type", "WaterWithin sites"),
                                 c("SoilWithin type", "WaterWithin type"),
                                 c("SoilWithin sites", "WaterWithin sites")),
              map_signif_level=TRUE,  y_position = c(0.7,0.7,0.75,0.81), color = "black", size = 1, textsize = 6) +
  coord_cartesian(ylim = c(0.15, 0.85)) + 
  labs(y = "Distance to centroid") +
  scale_color_manual(values=c("sienna4", "steelblue1")) + 
  scale_x_discrete(labels= c("Medium\ntype", "Site", "Medium\ntype", "Site")) +
  theme_pubr() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(2, 'cm'))

ggarrange(plot1, plot2, common.legend = T, align = "hv") 

# We repeat the same steps as before but for soil and water samples individually 
# Figure S4
soil.nmds <- metaMDS(soil.df[, -c(which(colnames(soil.df) == "Site_No"): 
                                    which(colnames(soil.df) == "fraction_assigned_reads"))],
                     distance = "jaccard",
                     k = 2,
                     trymax = 1000, 
                     wascores = TRUE)
adonis2(soil.df[, -c(which(colnames(soil.df) == "Site_No"): 
                      which(colnames(soil.df) == "fraction_assigned_reads"))] ~ soil.df$Site_No, method = "jaccard")
soil.scrs <- as.data.frame(scores(soil.nmds, display = "sites")) 
soil.scrs$Site_No <- as.numeric(soil.df$Site_No)
soil.scrs <- soil.scrs[order(soil.scrs$Site_No), ]
soil.scrs$Site_No <- as.factor(soil.scrs$Site_No)

water.nmds <- metaMDS(water.df[, -c(which(colnames(water.df) == "Site_No"): 
                                      which(colnames(water.df) == "fraction_assigned_reads"))],
                      distance = "jaccard",
                      k = 2,
                      trymax = 1000, 
                      wascores = TRUE)
adonis2(water.df[, -c(which(colnames(water.df) == "Site_No"): 
                      which(colnames(water.df) == "fraction_assigned_reads"))] ~ water.df$Site_No, method = "jaccard")
water.scrs <- as.data.frame(scores(water.nmds, display = "sites")) 
water.scrs$Site_No <- as.numeric(water.df$Site_No)
water.scrs <- water.scrs[order(water.scrs$Site_No), ]
water.scrs$Site_No <- as.factor(water.scrs$Site_No)

plot1 <- ggplot(soil.scrs, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(NMDS1, NMDS2, color = Site_No), size = 2) +
  geom_mark_hull(aes(fill = Site_No, color = Site_No), expand = 0, radius = 0, alpha = 0.1) + 
  annotate(geom = "text", label = paste0("stress = ", signif(soil.nmds$stress, 4)), x=-0.13, y=0.15, size = 6) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.margin=unit(c(2,0.2,0.2,0.2), 'cm')) + 
  guides(colour = guide_legend(nrow = 3)) +
  scale_color_discrete(breaks=soil.scrs$Site_No) +
  scale_fill_discrete(breaks=soil.scrs$Site_No)

plot2 <- ggplot(water.scrs, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(NMDS1, NMDS2, color = Site_No), size = 2) +
  geom_mark_hull(aes(fill = Site_No, color = Site_No), expand = 0, radius = 0, alpha = 0.1) +  
  annotate(geom = "text", label = paste0("stress = ", signif(water.nmds$stress, 4)), x=0.35, y=0.67, size = 6) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.margin=unit(c(2,0.2,0.2,0.2), 'cm')) + 
  guides(colour = guide_legend(nrow = 3)) +
  scale_color_discrete(breaks=soil.scrs$Site_No) +
  scale_fill_discrete(breaks=soil.scrs$Site_No) 

ggarrange(plot1, plot2, nrow = 1, common.legend = T, legend = "bottom",
          labels = c("Soil", "Water"),
          font.label = list(size = 18))


################################################################################
######################### Community structure diversity ########################
################################################################################

otu.df <- read.csv("Urban_BGB_OTU.csv")
otu.df$otu_richness <- specnumber(otu.df[, -which(colnames(otu.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type"))])
soil.df <- otu.df[otu.df$Type == "Soil", ]
water.df <- otu.df[otu.df$Type == "Water", ]

# Aggregate the data by site and by type and calculate the mean and standard deviation of OTU richness for each site
# The data has first to be separated into soil and water samples so that we do not merge all six samples into a single value. 
terr.df <- aggregate(soil.df[, -which(colnames(soil.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type", "otu_richness"))], 
                     list(soil.df$Site_No), sum)
terr.df$otu_richness_mean <- aggregate(soil.df$otu_richness, list(soil.df$Site_No), mean)[,2]
terr.df$otu_richness_sd <- aggregate(soil.df$otu_richness, list(soil.df$Site_No), sd)[,2]
terr.df$Type <- "Soil"
colnames(terr.df)[1] <- "Site_No"

aqua.df <- aggregate(water.df[, -which(colnames(water.df) %in% c("Site_No", "Uniq_cd", "Rplct_N", "Type", "otu_richness"))], 
                     list(water.df$Site_No), sum)
aqua.df$otu_richness_mean <- aggregate(water.df$otu_richness, list(water.df$Site_No), mean)[,2]
aqua.df$otu_richness_sd <- aggregate(water.df$otu_richness, list(water.df$Site_No), sd)[,2]
aqua.df$Type <- "Water"
colnames(aqua.df)[1] <- "Site_No"

bgb.df <- rbind(terr.df, aqua.df)

# Calculate Jaccard distance, the number of shared, and distinct OTU per site
bgb.shared.df <- as.data.frame(unique(bgb.df$Site_No))
colnames(bgb.shared.df) <- "Site_No"
bgb.shared.df$jaccard <- NA
bgb.shared.df$nb.otu.shared <- NA
bgb.shared.df$nb.otu.terr <- NA
bgb.shared.df$nb.otu.aqua <- NA

# Vector that will store the name of each shared OTU
otu.shared <- c()
for (site in unique(bgb.df$Site_No )) {
  nb.otu.shared <- 0
  nb.otu.terr <- 0
  nb.otu.aqua <- 0
  
  site_ij <- bgb.df[bgb.df$Site_No == site, ]
  raw <- site_ij[site_ij$Site_No == i, 
          -c(which(colnames(site_ij) == "Site_No"), which(colnames(site_ij) == "otu_richness_mean"):which(colnames(site_ij) == "Type"))] # Keeping only the OTU table
  
  raw[raw > 1] <- 1
  raw <- raw[, -which(colSums(raw) == 0)]
  i <- 0 
  for (sp in raw) {
    i <- i + 1
    if (sp[1] > 0 &  sp[2] > 0) { # sp[1] is always soil, sp[2] is always water. If both are present, increase the number of shared OTUs. Otherwise, increse the number of terrestrial- or aquatic-specific OTUs. 
      otu.shared <- append(otu.shared, colnames(raw)[i])
      nb.otu.shared <- nb.otu.shared + 1
    }
    if (sp[1] > 0 & sp[2] == 0) {
      nb.otu.terr <- nb.otu.terr + 1
    }
    if (sp[1] == 0 & sp[2] > 0) {
      nb.otu.aqua <- nb.otu.aqua + 1
    }
  }
  jaccard <- vegdist(raw, method = "jaccard")
  bgb.shared.df[bgb.shared.df$Site_No == site, "jaccard"] <- as.numeric(jaccard)
  bgb.shared.df[bgb.shared.df$Site_No == site, "nb.otu.shared"] <- nb.otu.shared
  bgb.shared.df[bgb.shared.df$Site_No == site, "nb.otu.terr"] <- nb.otu.terr
  bgb.shared.df[bgb.shared.df$Site_No == site, "nb.otu.aqua"] <- nb.otu.aqua
}

# Filter dataset for shared OTUs based on the vector with the names of all shared OTUs
otu.shared.df <- bgb.df[, which(colnames(bgb.df) %in% unique(otu.shared))]
otu.shared.df$Site_No <- bgb.df$Site_No

# Ensure shared OTUs have the same presence/absence across both types for each site
for (site in unique(otu.shared.df$Site_No)) {
  site_ij <- otu.shared.df[otu.shared.df$Site_No == site, ]
  site_ij <- site_ij[, -which(colnames(site_ij) == "Site_No")]
  for (i in 1:ncol(site_ij)) {
    if (site_ij[1, i] == 0 | site_ij[2, i] == 0) {
      otu.shared.df[otu.shared.df$Site_No == site, i] <- 0
      otu.shared.df[otu.shared.df$Site_No == site, i] <- 0
    }
  }
}

# Aggregate shared OTUs data by site
otu.shared.df <- aggregate(otu.shared.df[, -which(colnames(otu.shared.df) == "Site_No")], by = list(otu.shared.df$Site_No), FUN = sum)
colnames(otu.shared.df)[1] <- "Site_No"
bgb.shared.df <- join(bgb.shared.df, otu.shared.df, by = "Site_No")

# Create alpha diversity dataframe for terrestrial and aquatic sites
terr.df <- bgb.df[bgb.df$Type == "Soil", ]
aqua.df <- bgb.df[bgb.df$Type == "Water", ]

terr.alpha <- data.frame(bgb.df[bgb.df$Type == "Soil", ]$otu_richness_mean, 
                         bgb.shared.df$jaccard, 
                         bgb.df[bgb.df$Type == "Soil", ]$Site_No)
terr.alpha <- terr.alpha[order(terr.alpha$bgb.df.bgb.df.Type.....Soil.....Site_No), ]
aqua.alpha <- data.frame(bgb.df[bgb.df$Type == "Water", ]$otu_richness_mean, 
                         bgb.df[bgb.df$Type == "Water", ]$Site_No)
aqua.alpha <- aqua.alpha[order(aqua.alpha$bgb.df.bgb.df.Type.....Water.....Site_No), ]

# Combine terrestrial and aquatic alpha diversity into a single dataframe for visualization
bgb.alpha <- cbind(terr.alpha, aqua.alpha$bgb.df.bgb.df.Type.....Water.....otu_richness_mean)
colnames(bgb.alpha) <- c("terr_otu",
                         "jaccard",
                         "Site_No", 
                         "aqua_otu")

alpha_plot <- ggplot(bgb.alpha, aes(x = (terr_otu), y = (aqua_otu))) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 1.3, alpha = 0.2) +
  geom_point(aes(size = 1 - jaccard, color = 1 - jaccard)) + 
  geom_smooth(method = "lm", color = "blue", alpha = 0.3, linewidth = 1.3) +
  labs(y = "Mean aquatic OTU richness", x = "Mean terrestrial OTU richness") +
  # coord_cartesian(xlim = c(83, 368), ylim = c(83, 368))
  scale_size_continuous(range = c(2,6)) +
  scale_color_gradient2(high = g, low = "white") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18), 
        legend.position = "bottom") + 
  guides(size = guide_legend(title = "Jaccard similarity"), 
         color = guide_legend(title = "Jaccard similarity"))

# For confidentially reasons and as many of the sampled ponds were on private ground, 
# we will only disclose data on pond location if asked specifically and if given a reason.

# Calculate pairwise Jaccard dissimilarities for beta diversity (terrestrial and aquatic)
terr.beta <- data.frame(as.vector(vegdist(terr.df[, -which(colnames(terr.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], 
                                          method = "jaccard", 
                                          binary = T)),
                        which(lower.tri(vegdist(terr.df[, -which(colnames(terr.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], # We only want the upper part of the matrix so that every combination of samples is only present once (and not twice if we take the entire matrix)
                                                method = "jaccard", 
                                                binary = T)), 
                              arr.ind = TRUE))
aqua.beta <- data.frame(as.vector(vegdist(aqua.df[, -which(colnames(aqua.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], 
                                          method = "jaccard", 
                                          binary = T)),
                        which(lower.tri(vegdist(aqua.df[, -which(colnames(aqua.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], 
                                                method = "jaccard", 
                                                binary = T)), 
                              arr.ind = TRUE))
colnames(terr.beta) <- c("terr_dist", "Site2", "Site1")
colnames(aqua.beta) <- c("aqua_dist", "Site2", "Site1")

# Merge terrestrial and aquatic beta distances into one dataframe
bgb.beta <- cbind(terr.beta, aqua.beta)
bgb.beta <- bgb.beta[, -c(5,6)]

# Linear regression between terrestrial and aquatic pairwise dissimilarities and confirmation with a Mantel test
summary(betareg::betareg(aqua_dist ~ scale(terr_dist), 
                         data = bgb.beta))
mantel(vegdist(terr.df[, -which(colnames(terr.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], 
               method = "jaccard", 
               binary = T), 
       vegdist(aqua.df[, -which(colnames(aqua.df) %in% c("Site_No", "Type", "otu_richness_mean", "otu_richness_sd"))], 
               method = "jaccard", 
               binary = T))

beta_plot <- ggplot(bgb.beta, aes(x = terr_dist, y = aqua_dist)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1.3, alpha = 0.2) +
  geom_point(alpha = 0.5, cex = 3) + 
  geom_smooth(method = "lm", color = "blue", alpha = 0.3, linewidth = 1.3) +
  labs(y = "Pairwise aquatic dissimilarity", x = "Pairwise terrestrial dissimilarity") +
  # coord_cartesian(xlim = c(83, 368), ylim = c(83, 368))
  scale_size_continuous(range = c(2,6)) +
  scale_color_gradient2(high = g, low = "white") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) 

# Figure 3
ggarrange(alpha_plot, beta_plot, common.legend = T,
          labels = c("A", "B"),
          font.label = list(size = 30))


################################################################################
################################# Venn diagram ################################# 
################################################################################
grid.newpage()

# Figure 4A

draw.pairwise.venn(area1 = mean(bgb.shared.df$nb.otu.terr) + mean(bgb.shared.df$nb.otu.shared) ,
                            area2 = mean(bgb.shared.df$nb.otu.aqua) + mean(bgb.shared.df$nb.otu.shared) ,
                            cross.area = mean(bgb.shared.df$nb.otu.shared) + sd(bgb.shared.df$nb.otu.shared),
                            col = c("sienna4", "steelblue1"),
                            alpha = 0,
                            lty = 2,
                            lwd = 2,
                            label.col = "transparent")

draw.pairwise.venn(area1 = mean(bgb.shared.df$nb.otu.terr) + mean(bgb.shared.df$nb.otu.shared) ,
                            area2 = mean(bgb.shared.df$nb.otu.aqua) + mean(bgb.shared.df$nb.otu.shared) ,
                            cross.area = mean(bgb.shared.df$nb.otu.shared) - sd(bgb.shared.df$nb.otu.shared),
                            col = c("sienna4", "steelblue1"),
                            alpha = 0,
                            lty = 2,
                            lwd = 2,
                            label.col = "transparent")

draw.pairwise.venn(area1 = mean(bgb.shared.df$nb.otu.terr) + mean(bgb.shared.df$nb.otu.shared),
                            area2 = mean(bgb.shared.df$nb.otu.aqua) + mean(bgb.shared.df$nb.otu.shared),
                            cross.area = mean(bgb.shared.df$nb.otu.shared),
                            col = c("sienna4", "steelblue1"),
                            fill = c("sienna4", "steelblue1"),
                            category = c("Soil", "Water"),
                            cat.fontfamily = rep("sans", 2),
                            cat.cex	= 1.2,
                            alpha = 0.2,
                            fontfamily = rep("sans", 3),
                            cex = 2)

# grid.export("Fig4_A.svg")

### The file is then modified in Inkscape so that all circles from the same medium type have the same size, but keeping their location regarding the cross-area. 


################################################################################
################################## Elastic Net #################################
################################################################################

covariates.df <- read.csv("Urban_BGB_covariates.csv")
bgb.df <- merge(bgb.df, covariates.df, by = "Site_No")

### For water samples only
aqua.df <- bgb.df[bgb.df$Type == "Water", ]

# Define the predictors to be used in the regression model
predictors <- c(
  "Grey50_frac",
  "Grey500_frac",
  "Green50_frac",
  "Green500_frac",
  "Blue50_frac",
  "Blue500_frac",
  "LowManagement50_frac",
  "LowManagement500_frac",
  "LowRandom50_frac",
  "LowRandom500_frac",
  "Distance_Water",
  "Distance_Pond",
  "Distance_Forest",
  "Connectivity_aquatic",
  "fish_PCR",
  "tot_volume_filtered",
  "population",
  "employment",
  "overwarming",
  "no2",
  "Walls")

# Create a dataframe of the predictors, scaling the numeric predictors
predictors.df <- aqua.df[, which(colnames(aqua.df) %in% predictors)]
predictors.df <- as.data.frame(scale(predictors.df))
predictors.df$fish_PCR <- aqua.df$fish_PCR # Add fish presence (TRUE-FALSE) to predictors
predictors.df$Walls <- aqua.df$Walls # Add walls presence (TRUE-FALSE) to predictors

# Initialize a dataframe to store cross-validation results for each alpha value
res <- data.frame(
  alpha = numeric(),
  lambda = numeric(),
  mse = numeric()
)

# Perform cross-validation for Elastic Net model across a range of alpha values
for (alpha in seq(0.1, 1, len = 10)) {
  cv.glmnet <- cv.glmnet(
    as.matrix(predictors.df), # Predictor matrix
    aqua.df$otu_richness_mean, # Response variable (OTU richness)
    alpha = alpha, # Elastic net mixing parameter
    nfold = 54, # LOO cross-validation
    family = "gaussian"
  )
  mse.min <- min(cv.glmnet$cvm) # Find the minimum MSE for the current alpha
  res <- rbind(res,
               data.frame(
                 alpha = alpha,
                 lambda = cv.glmnet$lambda.min,
                 mse = mse.min
               ))
  
}
# Fit the final Elastic Net model using the optimal lambda and alpha 
glnmnet <- glmnet(
  as.matrix(predictors.df),
  aqua.df$otu_richness_mean,
  alpha = 1,
  lambda = 7.985904,
  family = "gaussian"
)

# Fit a linear model using key predictors from the results of the Elastic Net
aqua.lm <- lm(aqua.df$otu_richness_mean ~ LowManagement50_frac + tot_volume_filtered,
                     data = predictors.df)

# Extract coefficients and standard errors of the linear model for latter purposes
aqua.lm.res <- data.frame(aqua.lm$coefficients[-1],  sqrt(diag(vcov(aqua.lm)))[-1], "water")

### Repeat the same steps for soil samples only
terr.df <- bgb.df[bgb.df$Type == "Soil",]

predictors <- c(
  "Grey50_frac",
  "Grey500_frac",
  "Green50_frac",
  "Green500_frac",
  "Blue50_frac",
  "Blue500_frac",
  "LowManagement50_frac",
  "LowManagement500_frac",
  "LowRandom50_frac",
  "LowRandom500_frac",
  "Distance_Water",
  "Distance_Pond",
  "Distance_Forest",
  "Connectivity_terrestrial",
  "fish_PCR",
  "population",
  "employment",
  "overwarming",
  "no2")
predictors.df <- terr.df[, which(colnames(terr.df) %in% predictors)]
predictors.df <- as.data.frame(scale(predictors.df))
predictors.df$fish_PCR <- terr.df$fish_PCR

res <- data.frame(
  alpha = numeric(),
  lambda = numeric(),
  mse = numeric()
)

for (alpha in seq(0.1, 1, len = 10)) {
  cv.glmnet <- cv.glmnet(
    as.matrix(predictors.df),
    terr.df$otu_richness_mean,
    alpha = alpha,
    nfold = 54,
    family = "gaussian"
  )
  mse.min <- min(cv.glmnet$cvm)
  res <- rbind(res,
               data.frame(
                 alpha = alpha,
                 lambda = cv.glmnet$lambda.min,
                 mse = mse.min
               ))
}

#Lowest AIC
glnmnet <- glmnet(
  as.matrix(predictors.df),
  terr.df$otu_richness_mean,
  alpha = 0.1,
  lambda = 69.590095,
  family = "gaussian"
)

terr.lm <- lm(terr.df$otu_richness_mean ~ Green500_frac,
                    data = predictors.df)
terr.lm.res <- data.frame(terr.lm$coefficients[-1], sqrt(diag(vcov(terr.lm)))[-1], "soil")

### Repeat the same steps for the jaccard index only
bgb.shared.df <- merge(bgb.shared.df, covariates.df, by = "Site_No")

predictors <- c(
  "Grey50_frac",
  "Grey500_frac",
  "Green50_frac",
  "Green500_frac",
  "Blue50_frac",
  "Blue500_frac",
  "LowManagement50_frac",
  "LowManagement500_frac",
  "LowRandom50_frac",
  "LowRandom500_frac",
  "Distance_Water",
  "Distance_Pond",
  "Distance_Forest",
  "Connectivity_aquatic",
  "Connectivity_terrestrial",
  "fish_PCR",
  "population",
  "employment",
  "overwarming",
  "no2",
  "meterstocc",
  "Walls",
  "tot_volume_filtered"
)

predictors.df <- bgb.shared.df[, which(colnames(bgb.shared.df) %in% predictors)]

predictors.df <- as.data.frame(scale(predictors.df))
predictors.df$fish_PCR <- bgb.shared.df$fish_PCR
predictors.df$Walls <- bgb.shared.df$Walls

# Transformation of the jaccard distance into the Jaccard similarity index (1 - Jaccard distance) and normalization using a Box-Cox transformation
bgb.shared.df$sim <- 1 - bgb.shared.df$jaccard
jaccard.lm <- lm(sim ~ 1, data = bgb.shared.df)
b <- boxcox(jaccard.lm)
lambda <- b$x[which.max(b$y)]
bgb.shared.df$sim <- (bgb.shared.df$sim ^ lambda - 1) / lambda

res <- data.frame(
  alpha = numeric(),
  lambda = numeric(),
  mse = numeric()
)

for (alpha in seq(0.1, 1, len = 10)) {
  cv.glmnet <- cv.glmnet(
    as.matrix(predictors.df),
    bgb.shared.df$sim,
    alpha = alpha,
    nfold = 54,
    family = "gaussian"
  )
  mse.min <- min(cv.glmnet$cvm)
  res <- rbind(res,
               data.frame(
                 alpha = alpha,
                 lambda = cv.glmnet$lambda.min,
                 mse = mse.min
               ))
}

#Lowest MSE
glnmnet <- glmnet(
  as.matrix(predictors.df),
  bgb.shared.df$sim,
  alpha = 0.4,
  lambda = 0.5043072,
  family = "gaussian"
)

jaccard.lm <- lm(bgb.shared.df$sim ~ Grey500_frac + population,
                       data = predictors.df)
jaccard.lm.res <- data.frame(jaccard.lm$coefficients[-1], sqrt(diag(vcov(jaccard.lm)))[-1], "jaccard")

# Figure 4B
# Regrouping the results from the linear models
aqua.lm.res$parameter <- rownames(aqua.lm.res)
terr.lm.res$parameter <- rownames(terr.lm.res)
jaccard.lm.res$parameter <- rownames(jaccard.lm.res)
colnames(aqua.lm.res) <- c("coef", "stderr", "model", "parameter")
colnames(terr.lm.res) <- c("coef", "stderr", "model", "parameter")
colnames(jaccard.lm.res) <- c("coef", "stderr", "model", "parameter")

# Add placeholder data for missing combinations of models and parameters
fig4.df.2 <- rbind(aqua.lm.res, terr.lm.res, jaccard.lm.res)
complete <- rbind(data.frame(coef = 0, 
                             stderr = 0, 
                             model = "soil", 
                             parameter = "tot_volume_filtered"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "soil", 
                             parameter = "LowManagement50_frac"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "soil", 
                             parameter = "Grey500_frac"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "soil", 
                             parameter = "population"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "water", 
                             parameter = "Green500_frac"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "water", 
                             parameter = "Grey500_frac"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "water", 
                             parameter = "population"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "jaccard", 
                             parameter = "LowManagement50_frac"),
                  data.frame(coef = 0, 
                             stderr = 0, 
                             model = "jaccard", 
                             parameter = "Green500_frac"))
fig4.df.2 <- rbind(fig4.df.2, complete)
fig4.df.2$model2 <- ifelse(fig4.df.2$model == "jaccard", "Similarity model", "Richness models")

fig4.df.2 <- data.table(fig4.df.2)

# Categorize models into "Similarity model" or "Richness models"
fig4.df.2[model == "water", x_min := -22]
fig4.df.2[model == "water", x_max := 22]
fig4.df.2[model == "soil", x_min := -22]
fig4.df.2[model == "soil", x_max := 22]
fig4.df.2[model == "jaccard", x_min := -0.61]
fig4.df.2[model == "jaccard", x_max := 0.61]

fig4.df.2$parameter <- c("Fraction of low-\nmaintained green\nspaces (50 m buffer)", "Volume filtered", "Fraction of green\nspaces (500 m buffer)", "Fraction of\nimpervious spaces\n(500 m buffer)", 
                         "Human population\ndensity", "Volume filtered", "Fraction of low-\nmaintained green\nspaces (50 m buffer)", "Fraction of\nimpervious spaces\n(500 m buffer)", 
                         "Human population\ndensity", "Fraction of green\nspaces (500 m buffer)", "Fraction of\nimpervious spaces\n(500 m buffer)", "Human population\ndensity", 
                         "Fraction of low-\nmaintained green\nspaces (50 m buffer)", "Fraction of green\nspaces (500 m buffer)")

fig4.2.plot <- ggplot(fig4.df.2, aes(x = coef, y = parameter, color = model)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr), width = 0.3, linewidth = 1.07, position = position_dodge(width=0.5)) +
  geom_point(size = 3, position = position_dodge(width=0.5)) +
  facet_grid(~ factor(model2, c("Richness models", "Similarity model")), scales = "free_x") +
  geom_blank(aes(x = x_min)) +
  geom_blank(aes(x = x_max)) +
  scale_color_manual(values = c(g, br, bl), labels = c("Jaccard similarity", "Terrestrial OTU richness", "Aquatic OTU richness")) +
  labs(x = "Coefficients", y = "Parameters") +
  theme_pubclean(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_text(size = 14), 
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 18)) + 
  guides(colour = guide_legend(title = ""))


#Figure 4C-F
fig4.df.3 <- bgb.shared.df
fig4.df.3$terr.otu <- terr.df$otu_richness_mean
fig4.df.3$aqua.otu <- aqua.df$otu_richness_mean

# We use a centered log.-ratio so that all variables can be displayed on the same scale. 
fig4.df.3$jaccard.clr <- matrix(clr(1 - fig4.df.3$jaccard), ncol = 1)
fig4.df.3$terr.otu.clr <- matrix(clr(fig4.df.3$terr.otu), ncol = 1)
fig4.df.3$aqua.otu.clr <- matrix(clr(fig4.df.3$aqua.otu), ncol = 1)

plot1 <- ggplot(data = fig4.df.3, aes(x = LowManagement50_frac)) + 
  geom_point(aes(y = terr.otu.clr, color = "Soil"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = terr.otu.clr, fill = "Soil"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = terr.otu.clr, color = "Soil"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = jaccard.clr, color = "Jaccard"), alpha = 1) +
  stat_smooth(geom="ribbon", aes(y = jaccard.clr, fill = "Jaccard"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = jaccard.clr, color = "Jaccard"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = aqua.otu.clr, color = "Water"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = aqua.otu.clr, fill = "Water"), method = "lm", alpha = 0.2) +
  stat_smooth(geom="line", aes(y = aqua.otu.clr, color = "Water"), method = "lm", alpha = 1, size = 1) +
  labs(x = "Fraction of low-maintained\ngreen spaces (50 m buffer)", y = "Centered log. ratio") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16)) +
  scale_color_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                     labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  scale_fill_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                    labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  guides(colour = guide_legend(title = ""), 
         fill = "none")

plot2 <- ggplot(data = fig4.df.3, aes(x = Green500_frac)) +
  geom_point(aes(y = aqua.otu.clr, color = "Water"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = aqua.otu.clr, fill = "Water"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = aqua.otu.clr, color = "Water"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = jaccard.clr, color = "Jaccard"), alpha = 1) +
  stat_smooth(geom="ribbon", aes(y = jaccard.clr, fill = "Jaccard"), method = "lm", alpha = 0.2) +
  stat_smooth(geom="line", aes(y = jaccard.clr, color = "Jaccard"), method = "lm", alpha = 1, size = 1) +
  geom_point(aes(y = terr.otu.clr, color = "Soil"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = terr.otu.clr, fill = "Soil"), method = "lm", alpha = 0.2) +
  stat_smooth(geom="line", aes(y = terr.otu.clr, color = "Soil"), method = "lm", alpha = 1, size = 1) +
  labs(x = "Fraction of green spaces\n (500 m buffer)", y = "Centered log. ratio") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_blank()) +
  scale_color_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                     labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  scale_fill_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                    labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  guides(colour = guide_legend(title = ""),
         fill = "none")

plot3 <- ggplot(data = fig4.df.3, aes(x = Grey500_frac)) + 
  geom_point(aes(y = aqua.otu.clr, color = "Water"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = aqua.otu.clr, fill = "Water"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = aqua.otu.clr, color = "Water"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = terr.otu.clr, color = "Soil"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = terr.otu.clr, fill = "Soil"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = terr.otu.clr, color = "Soil"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = jaccard.clr, color = "Jaccard"), alpha = 1) +
  stat_smooth(geom="ribbon", aes(y = jaccard.clr, fill = "Jaccard"), method = "lm", alpha = 0.2) +
  stat_smooth(geom="line", aes(y = jaccard.clr, color = "Jaccard"), method = "lm", alpha = 1, size = 1) +
  labs(x = "Fraction of grey spaces\n (500 m buffer)", y = "Centered log. ratio") +
  theme_pubr() + 
  theme(legend.text = element_markdown(size = 18), 
        legend.title = element_markdown(size = 18), 
        axis.text = element_markdown(size = 14), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_blank()) +
  scale_color_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                     labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  scale_fill_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                    labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  guides(colour = guide_legend(title = ""),
         fill = "none")

plot4 <- ggplot(data = fig4.df.3, aes(x = population)) + 
  geom_point(aes(y = aqua.otu.clr, color = "Water"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = aqua.otu.clr, fill = "Water"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = aqua.otu.clr, color = "Water"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = terr.otu.clr, color = "Soil"), alpha = 0.5) +
  stat_smooth(geom="ribbon", aes(y = terr.otu.clr, fill = "Soil"), method = "lm", alpha = 0.1) +
  stat_smooth(geom="line", aes(y = terr.otu.clr, color = "Soil"), method = "lm", alpha = 0.4, size = 1) +
  geom_point(aes(y = jaccard.clr, color = "Jaccard"), alpha = 1) +
  stat_smooth(geom="ribbon", aes(y = jaccard.clr, fill = "Jaccard"), method = "lm", alpha = 0.2) +
  stat_smooth(geom="line", aes(y = jaccard.clr, color = "Jaccard"), method = "lm", alpha = 1, size = 1) +
  labs(x = "Human population density<br>(inhabitants per 100 m<sup>2</sup>)", y = "Centered log. ratio") +
  theme_pubr() + 
  theme(legend.text = element_markdown(size = 18), 
        legend.title = element_markdown(size = 18), 
        axis.text = element_markdown(size = 14), 
        axis.title.x = element_markdown(size = 16), 
        axis.title.y = element_blank()) +
  scale_color_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                     labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  scale_fill_manual(values = c(Jaccard = g, Soil = br, Water = bl),
                    labels = c(Jaccard = "Jaccard similarity", Soil = "Terrestrial OTU richness", Water = "Aquatic OTU richness")) + 
  guides(colour = guide_legend(title = ""), 
         fill = "none")

fig4.3.plot <- ggarrange(plot1, plot2, plot3, plot4, 
                         nrow = 1, 
                         ncol = 4, 
                         align = "h",
                         labels = c("C", "D", "E", "F"),
                         common.legend = T) + 
  theme(plot.margin = unit(c(0.5, 0, 0, 0), 'cm')) +
  guides(colour = guide_legend(title = ""),
         fill = "none")

ggarrange(fig4.2.plot, fig4.3.plot, nrow = 2, ncol = 1, common.legend = T,
          heights  = c(1, 0.6),
          labels = c("B", ""),
          font.label = list(size = 30)) + 
  guides(colour = guide_legend(title = ""), 
         fill = "none")


################################################################################
############################ individual correlations ###########################
################################################################################
# Figure S6

plot_manag1 <- ggplot(data = aqua.df, aes(y = otu_richness_mean, x = LowManagement50_frac)) + 
  geom_point(cex = 3, color = bl) + 
  geom_smooth(method = "lm", alpha = 0.2, color = bl) +
  labs(x = "Fraction of low-maintained\ngreen spaces (50m buffer)", y = "Aquatic OTU richness") +
  scale_color_manual(values = c(bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_manag2 <- ggplot(data = terr.df, aes(y = otu_richness_mean, x = LowManagement50_frac)) + 
  geom_point(cex = 3, color = br) + 
  geom_smooth(method = "lm", alpha = 0.2, color = br) +
  labs(x = "Fraction of low-maintained\ngreen spaces (50m buffer)", y = "Terrestrial OTU richness") +
  scale_color_manual(values = c(br, bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_manag3 <- ggplot(data = bgb.shared.df, aes(y = 1- jaccard, x = LowManagement50_frac)) + 
  geom_point(cex = 3, color = g) + 
  geom_smooth(method = "lm", alpha = 0.2, color = g) +
  labs(x = "Fraction of low-maintained\ngreen spaces (50m buffer)", y = "Jaccard similarity") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_green1 <- ggplot(data = aqua.df, aes(y = otu_richness_mean, x = Green500_frac)) + 
  geom_point(cex = 3, color = bl) + 
  geom_smooth(method = "lm", alpha = 0.2, color = bl) +
  labs(x = "Fraction of green\nspaces (500m buffer)", y = "Aquatic OTU richness") +
  scale_color_manual(values = c(bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_green2 <- ggplot(data = terr.df, aes(y = otu_richness_mean, x = Green500_frac)) + 
  geom_point(cex = 3, color = br) + 
  geom_smooth(method = "lm", alpha = 0.2, color = br) +
  labs(x = "Fraction of green\nspaces (500m buffer)", y = "Terrestrial OTU richness") +
  scale_color_manual(values = c(br, bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_green3 <- ggplot(data = bgb.shared.df, aes(y = 1- jaccard, x = Green500_frac)) + 
  geom_point(cex = 3, color = g) + 
  geom_smooth(method = "lm", alpha = 0.2, color = g) +
  labs(x = "Fraction of green\nspaces (500m buffer)", y = "Jaccard similarity") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_grey1 <- ggplot(data = aqua.df, aes(y = otu_richness_mean, x = Grey500_frac)) + 
  geom_point(cex = 3, color = bl) + 
  geom_smooth(method = "lm", alpha = 0.2, color = bl) +
  labs(x = "Fraction of impervious\nsurface (500m buffer)", y = "Aquatic OTU richness") +
  scale_color_manual(values = c(bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_grey2 <- ggplot(data = terr.df, aes(y = otu_richness_mean, x = Grey500_frac)) + 
  geom_point(cex = 3, color = br) + 
  geom_smooth(method = "lm", alpha = 0.2, color = br) +
  labs(x = "Fraction of impervious\nsurface (500m buffer)", y = "Terrestrial OTU richness") +
  scale_color_manual(values = c(br, bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_grey3 <- ggplot(data = bgb.shared.df, aes(y = 1- jaccard, x = Grey500_frac)) + 
  geom_point(cex = 3, color = g) + 
  geom_smooth(method = "lm", alpha = 0.2, color = g) +
  labs(x = "Fraction of impervious\nsurface (500m buffer)", y = "Jaccard similarity") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_pop1 <- ggplot(data = aqua.df, aes(y = otu_richness_mean, x = population)) + 
  geom_point(cex = 3, color = bl) + 
  geom_smooth(method = "lm", alpha = 0.2, color = bl) +
  labs(x = "Population density", y = "Aquatic OTU richness") +
  scale_color_manual(values = c(bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_pop2 <- ggplot(data = terr.df, aes(y = otu_richness_mean, x = population)) + 
  geom_point(cex = 3, color = br) + 
  geom_smooth(method = "lm", alpha = 0.2, color = br) +
  labs(x = "Population density", y = "Terrestrial OTU richness") +
  scale_color_manual(values = c(br, bl)) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")
plot_pop3 <- ggplot(data = bgb.shared.df, aes(y = 1- jaccard, x = population)) + 
  geom_point(cex = 3, color = g) + 
  geom_smooth(method = "lm", alpha = 0.2, color = g) +
  labs(x = "Population density", y = "Jaccard similarity") +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_markdown(size = 18),
        legend.position="bottom")

plot_manag_tot <- ggarrange(plot_manag1 + rremove("xlab"), plot_manag2 + rremove("xlab"), plot_manag3 + rremove("xlab"),
                            nrow = 1, ncol = 3, labels = NULL)
plot_manag_tot <- annotate_figure(plot_manag_tot, bottom = text_grob("Fraction of low-maintained green spaces (50m buffer)", size = 18))
plot_green_tot <- ggarrange(plot_green1 + rremove("xlab"), plot_green2 + rremove("xlab"), plot_green3 + rremove("xlab"),
                            nrow = 1, ncol = 3, labels = NULL)
plot_green_tot <- annotate_figure(plot_green_tot, bottom = text_grob("Fraction of green spaces (500m buffer)", size = 18))
plot_grey_tot <- ggarrange(plot_grey1 + rremove("xlab"), plot_grey2 + rremove("xlab"), plot_grey3 + rremove("xlab"),
                           nrow = 1, ncol = 3, labels = NULL)
plot_grey_tot <- annotate_figure(plot_grey_tot, bottom = text_grob("Fraction of impervious surface (500m buffer)", size = 18))
plot_pop_tot <- ggarrange(plot_pop1 + rremove("xlab"), plot_pop2 + rremove("xlab"), plot_pop3 + rremove("xlab"),
                          nrow = 1, ncol = 3, labels = NULL)
plot_pop_tot <- annotate_figure(plot_pop_tot, bottom = text_grob(expression(paste("Human population density (inhabitants per 100 m"^2, ")")), size = 18))

ggarrange(plot_manag_tot, 
          plot_green_tot, 
          plot_grey_tot,
          plot_pop_tot, 
          nrow = 4, ncol = 1)

################################################################################
############################## Rarefaction curves ##############################
################################################################################
# Figure S5

# Initialize an empty list to store the data for water sites
water.accu.list <- list()
for (i in unique(water.df$Site_No)) {
  site <- water.df[water.df$Site_No == i, 
                   -c(which(colnames(water.df) == "Site_No"):which(colnames(water.df) == "otu_richness"))]
  site[site > 1] <- 1
  list <- list(t(as.data.frame(site)))
  water.accu.list <- append(water.accu.list, list)
}
water.accu.inext <- iNEXT(water.accu.list, q = 0, datatype = "incidence_raw") # Apply the iNEXT function to the water data list to estimate species richness
water.accu.df <- fortify(water.accu.inext, type = 1) # Convert the output of iNEXT into a data frame suitable for plotting
water.accu.df$Type <- "Water"

# Repeat the same process for soil samples
soil.accu.list <- list()
for (i in unique(soil.df$Site_No)) {
  site <- soil.df[soil.df$Site_No == i, 
                  -c(which(colnames(soil.df) == "Site_No"):which(colnames(soil.df) == "otu_richness"))]
  site[site > 1] <- 1
  list <- list(t(as.data.frame(site)))
  soil.accu.list <- append(soil.accu.list, list)
}
soil.accu.inext <- iNEXT(soil.accu.list, q = 0, datatype = "incidence_raw")
soil.accu.df <- fortify(soil.accu.inext, type = 1)
soil.accu.df$Type <- "Soil"

# Combine the data frames for water and soil into a single data frame
otu.accu.df <- rbind(water.accu.df, soil.accu.df)
otu.accu.point <- otu.accu.df[which(otu.accu.df$Method=="Observed"), ] # Extract the "observed" data points (i.e., the maximal number of samples per (= 3))
otu.accu.point_simple <- otu.accu.point # Create a simplified version of the observed data by averaging values for each type so that they can be displayed with different shapes latter on
otu.accu.point_simple[otu.accu.point_simple$Type == "Water", 7] <- mean(otu.accu.point_simple[otu.accu.point_simple$Type == "Water", 7])
otu.accu.point_simple[otu.accu.point_simple$Type == "Soil", 7] <- mean(otu.accu.point_simple[otu.accu.point_simple$Type == "Soil", 7])

# Simplify the water data by computing mean and standard deviation for each x value so that we can create a ribbon using these values latter on
water.accu.simple <- otu.accu.df[otu.accu.df$Type == "Water", ]
for (i in unique(water.accu.simple$x)) {
  mean <- mean(water.accu.simple[water.accu.simple$x == i, 7])
  sd <- sd(water.accu.simple[water.accu.simple$x == i, 7])
  water.accu.simple[water.accu.simple$x == i, 7] <- mean
  water.accu.simple[water.accu.simple$x == i, 8] <- mean - sd 
  water.accu.simple[water.accu.simple$x == i, 9] <- mean + sd
}
# Do the same for soil samples
soil.accu.simple <- otu.accu.df[otu.accu.df$Type == "Soil", ]
for (i in unique(soil.accu.simple$x)) {
  mean <- mean(soil.accu.simple[soil.accu.simple$x == i, 7])
  sd <- sd(soil.accu.simple[soil.accu.simple$x == i, 7])
  soil.accu.simple[soil.accu.simple$x == i, 7] <- mean
  soil.accu.simple[soil.accu.simple$x == i, 8] <- mean - sd
  soil.accu.simple[soil.accu.simple$x == i, 9] <- mean + sd
}
otu.accu.simple <- rbind(water.accu.simple, soil.accu.simple)
otu.accu.simple <- (unique(otu.accu.simple[, c(6:10)]))

# Separate the data into intra-site and extra-site categories based on x values so that rarefied and extrapolated data can be plotted with different line types 
otu.accu.df_intra <- otu.accu.df[otu.accu.df$x <= 3, ]
otu.accu.df_extra <- otu.accu.df[otu.accu.df$x >= 3, ]
otu.accu.df_intra$AssemblageType <- paste0(otu.accu.df_intra$Assemblage, otu.accu.df_intra$Type) # Do the same for water and soil samples so that they can be ploted with different shapes
otu.accu.df_extra$AssemblageType <- paste0(otu.accu.df_extra$Assemblage, otu.accu.df_extra$Type)
otu.accu.simple_intra <- otu.accu.simple[otu.accu.simple$x <= 3, ]
otu.accu.simple_extra <- otu.accu.simple[otu.accu.simple$x >= 3, ]

ggplot(otu.accu.df, aes(x = x, y = y)) +
  geom_line(aes(fill = Type, colour = AssemblageType), linetype = "solid", lwd = 0.8, data = otu.accu.df_intra, alpha = 0.4) + 
  geom_line(aes(fill = Type, colour = AssemblageType), linetype = "dashed", lwd = 0.8, data = otu.accu.df_extra, alpha = 0.4) +
  geom_point(data = otu.accu.point, aes(shape = Type, color = Type), size = 3, alpha = 0.4) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Type), alpha = 0.3, data = otu.accu.simple) + 
  geom_line(aes(fill = Type), linetype = "solid", lwd = 1.2, data = otu.accu.simple_intra) + 
  geom_line(aes(fill = Type), linetype = "dashed", lwd = 1.2, data = otu.accu.simple_extra) + 
  geom_point(data = otu.accu.point_simple, aes(shape = Type), size = 5) +
  xlab("Number of samples per site") + 
  ylab("OTU richness") + 
  scale_color_manual(values = c(rep(c("sienna4", "steelblue1"), 55))) +
  scale_fill_manual(values=c("sienna4", "steelblue1")) + 
  guides(color = FALSE, 
         fill = guide_legend(title = "Medium type"),
         shape = guide_legend(title = "Medium type")) +
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) 


################################################################################
########################## swarms composition analysis #########################
################################################################################

rownames(taxo.df) <- taxo.df$otu_id
LH.df <- read.csv("Urban_BGB_life_history.csv")
rownames(LH.df) <- LH.df$otu_id

# Function to classify taxa based on their name and other details
classify_taxa <- function(name, sp_name, reads) {
  # Initialize a vector with zeros for each taxonomic classification
  taxa <- c(unclassified = 0, Annelida = 0, Else = 0, Collembola = 0, Arachnida = 0,
            Crustacea = 0, Myriapoda = 0, Insecta = 0, unclassified_arthropods = 0,
            unclassified_frac = 0, Annelida_frac = 0, Else_frac = 0, Collembola_frac = 0,
            Arachnida_frac = 0, Crustacea_frac = 0, Myriapoda_frac = 0, Insecta_frac = 0,
            unclassified_arthropods_frac = 0)
  
  # Classify based on the name of the taxon
  if (is.na(name) || name == "unclassified") {
    taxa["unclassified"] <- 1
    taxa["unclassified_frac"] <- reads
  } else if (name == "Annelida") {
    taxa["Annelida"] <- 1
    taxa["Annelida_frac"] <- reads
  } else if (name %in% c("Chordata", "Rotifera", "Discosea", 
                         "Mollusca", "Gastrotricha", "Oomycota", "Cnidaria")) {
    taxa["Else"] <- 1
    taxa["Else_frac"] <- reads
  } else if (name == "Arthropoda") {
    # Further classify within the Arthropoda phylum
    if (is.na(sp_name) || sp_name == "unclassified") {
      taxa["unclassified_arthropods"] <- 1
      taxa["unclassified_arthropods_frac"] <- reads
    } else if (sp_name == "Collembola") {
      taxa["Collembola"] <- 1
      taxa["Collembola_frac"] <- reads
    } else if (sp_name == "Arachnida") {
      taxa["Arachnida"] <- 1
      taxa["Arachnida_frac"] <- reads
    } else if (sp_name %in% c("Ostracoda", "Branchiopoda", 
                              "Hexanauplia", "Malacostraca")) {
      taxa["Crustacea"] <- 1
      taxa["Crustacea_frac"] <- reads
    } else if (sp_name %in% c("Diplopoda", "Chilopoda", "Symphyla")) {
      taxa["Myriapoda"] <- 1
      taxa["Myriapoda_frac"] <- reads
    } else if (sp_name == "Insecta") {
      taxa["Insecta"] <- 1
      taxa["Insecta_frac"] <- reads
    }
  }
  return(taxa)
}

# Function to classify lifecycle stages
classify_lifecycle <- function(life, reads) {
  # Initialize a vector with zeros for each lifecycle classification
  taxa <- c(Aquatic = 0, Terrestrial = 0, Both = 0, na = 0, 
            Aquatic_frac = 0, Terrestrial_frac = 0, Both_frac = 0, na_frac = 0)
  # Classify based on lifecycle stage
  if (is.na(life)) {
    taxa["na"] <- 1
    taxa["na_frac"] <- reads
  } else if (life == "Aquatic") {
    taxa["Aquatic"] <- 1
    taxa["Aquatic_frac"] <- reads
  } else if (life == "Terrestrial") {
    taxa["Terrestrial"] <- 1
    taxa["Terrestrial_frac"] <- reads
  } else if (life == "Semi-terrestrial" | life == "Semi-aquatic") {
    taxa["Both"] <- 1
    taxa["Both_frac"] <- reads
  }
  return(taxa)
}

# Process soil samples
terr.compo <- do.call(rbind, lapply(1:nrow(terr.df), function(i) {
  site_no <- terr.df$Site_No[i]
  site <- terr.df[terr.df$Site_No == i, 
                  -c(which(colnames(terr.df) == "Site_No"), which(colnames(terr.df) == "otu_richness_mean"):which(colnames(terr.df) == "Type"))]
  site <- site[, colSums(site) != 0] # Remove columns with all zero values
  taxa_counts <- colSums(do.call(rbind, lapply(1:ncol(site), function(j) {  # Classify each species in the site data
    sp <- colnames(site)[j]
    name <- taxo.df[sp, "Phylum"]
    sp_name <- taxo.df[sp, "Class"]
    reads <- as.numeric(site[, j])
    classify_taxa(name, sp_name, reads)
  })))
  res <- data.frame(rowid = i, t(taxa_counts)) # Create a data frame with the counts for each taxa and return
  do.call(rbind, lapply(names(taxa_counts)[1:9], function(taxa) {
    data.frame(site_no, taxo = taxa, richness = as.numeric(res[taxa]), reads = as.numeric(res[paste0(taxa, "_frac")]))
  }))
}))

# Do the same for water samples
aqua.compo <- do.call(rbind, lapply(1:nrow(aqua.df), function(i) {
  site_no <- aqua.df$Site_No[i]
  site <- aqua.df[aqua.df$Site_No == i, 
                  -c(which(colnames(aqua.df) == "Site_No"), which(colnames(aqua.df) == "otu_richness_mean"):which(colnames(aqua.df) == "Type"))]
  site <- site[, colSums(site) != 0]
  taxa_counts <- colSums(do.call(rbind, lapply(1:ncol(site), function(j) {
    sp <- colnames(site)[j]
    name <- taxo.df[sp, "Phylum"]
    sp_name <- taxo.df[sp, "Class"]
    reads <- as.numeric(site[, j])
    classify_taxa(name, sp_name, reads)
  })))
  res <- data.frame(rowid = i, t(taxa_counts))
  do.call(rbind, lapply(names(taxa_counts)[1:9], function(taxa) {
    data.frame(site_no, taxo = taxa, richness = as.numeric(res[taxa]), reads = as.numeric(res[paste0(taxa, "_frac")]))
  }))
}))

# Do the same for shared OTUs
shared.compo <- do.call(rbind, lapply(1:nrow(bgb.shared.df), function(i) {
  site_no <- bgb.shared.df$Site_No[i]
  site <- bgb.shared.df[bgb.shared.df$Site_No == i, 
                        -c(which(colnames(bgb.shared.df) == "Site_No"):which(colnames(bgb.shared.df) == "nb.otu.aqua"))]
  site <- site[, colSums(site) != 0]
  taxa_counts <- colSums(do.call(rbind, lapply(1:ncol(site), function(j) {
    sp <- colnames(site)[j]
    name <- taxo.df[sp, "Phylum"]
    sp_name <- taxo.df[sp, "Class"]
    reads <- as.numeric(site[, j])
    classify_taxa(name, sp_name, reads)
  })))
  res <- data.frame(rowid = i, t(taxa_counts))
  do.call(rbind, lapply(names(taxa_counts)[1:9], function(taxa) {
    data.frame(site_no, taxo = taxa, richness = as.numeric(res[taxa]), reads = as.numeric(res[paste0(taxa, "_frac")]))
  }))
}))

# Do the same for shared OTUs, this time classifying their life cycle
shared.life <- do.call(rbind, lapply(1:nrow(bgb.shared.df), function(i) {
  site_no <- bgb.shared.df$Site_No[i]
  site <- bgb.shared.df[bgb.shared.df$Site_No == i, 
                        -c(which(colnames(bgb.shared.df) == "Site_No"):which(colnames(bgb.shared.df) == "nb.otu.aqua"))]
  site <- site[, colSums(site) != 0]
  taxa_counts <- colSums(do.call(rbind, lapply(1:ncol(site), function(j) {
    sp <- colnames(site)[j]
    life <- LH.df[sp, "Life_cycle"]
    reads <- as.numeric(site[, j])
    classify_lifecycle(life, reads)
  })))
  res <- data.frame(Site_No = i, t(taxa_counts))
  do.call(rbind, lapply(names(taxa_counts)[1:4], function(taxa) {
    data.frame(site_no, taxo = taxa, richness = as.numeric(res[taxa]), reads = as.numeric(res[paste0(taxa, "_frac")]))
  }))
}))

colnames(terr.compo) <- c("Site_No", "taxa", "richness", "reads")
colnames(aqua.compo) <- c("Site_No", "taxa", "richness", "reads")
colnames(shared.compo) <- c("Site_No", "taxa", "richness", "reads")
colnames(shared.life) <- c("Site_No", "LifeCycle", "richness", "reads")

terr.compo <- transform(terr.compo[order(as.numeric(terr.compo$Site_No)), ], Site_No = as.factor(Site_No))
aqua.compo <- transform(aqua.compo[order(as.numeric(aqua.compo$Site_No)), ], Site_No = as.factor(Site_No))
shared.compo <- transform(shared.compo[order(as.numeric(shared.compo$Site_No)), ], Site_No = as.factor(Site_No))
shared.life <- transform(shared.life[order(as.numeric(shared.life$Site_No)), ], Site_No = as.factor(Site_No))

# Figure S3

soil.plot_frac <- ggplot(terr.compo[terr.compo$taxa != "unclassified", ], aes(x = Site_No, y = richness, fill = taxa, color = taxa)) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "Site", y = "Fraction of OTUs") + 
  scale_fill_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                    labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  scale_color_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                     labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  theme_pubr() + 
  theme(axis.text.x = element_blank(), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.margin=unit(c(1,0.2,0.2,0.2), 'cm')) 

water.plot_frac <- ggplot(aqua.compo[aqua.compo$taxa != "unclassified", ], aes(x = Site_No, y = richness, fill = taxa, color = taxa)) + 
  geom_bar(position = "fill", stat = "identity", ) + 
  labs(x = "Site", y = "Fraction of OTUs") + 
  scale_fill_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                    labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  scale_color_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                     labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  theme_pubr() + 
  theme(axis.text.x = element_blank(), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.margin=unit(c(1,0.2,0.2,0.2), 'cm')) 

shared.plot_frac <- ggplot(shared.compo[shared.compo$taxa != "unclassified", ], aes(x = Site_No, y = richness, fill = taxa, color = taxa)) + 
  geom_bar(position = "fill", stat = "identity", ) + 
  labs(x = "Site", y = "Fraction of OTUs") + 
  scale_fill_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                    labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  scale_color_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                     labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  theme_pubr() + 
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.margin=unit(c(1,0.2,0.2,0.2), 'cm')) 

ggarrange(soil.plot_frac, water.plot_frac, shared.plot_frac, nrow = 3, common.legend = T,
          labels = c("Soil", "Water", "Shared"), 
          font.label = list(size = 18)) 

# Figure S7

terr.compo$Type <- "Soil"
aqua.compo$Type <- "Water"
shared.compo$Type <- "Shared"
shared.life$Type <- "Shared"

otu.compo <- rbind(terr.compo, aqua.compo, shared.compo)

plot1 <- ggplot(otu.compo[otu.compo$taxa != "unclassified", ], aes(x = Type, y = richness, fill = taxa, color = taxa)) + 
  geom_bar(position = "fill", stat = "identity", ) + 
  labs(x = "", y = "Fraction of assigned OTUs") + 
  scale_x_discrete(labels= c("Shared", "Soil", "Water")) +
  scale_color_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                     labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  scale_fill_manual(values = c("chocolate4", "firebrick", "goldenrod", "cornflowerblue", "lightgrey", "forestgreen", "sienna2", "gray40"), 
                    labels = c("Annelida", "Arachnida", "Collembola", "Crustacea", "Else", "Insecta", "Myriapoda", "Unclassified arthropods")) + 
  theme_pubr() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))

plot2 <- ggplot(shared.life[shared.life$LifeCycle != "na", ], aes(x = Type, y = richness, fill = LifeCycle, color = LifeCycle)) + 
  geom_bar(position = "fill", stat = "identity", ) + 
  labs(x = "", y = "Fraction of OTUs") + 
  scale_color_manual(values = c(bl, g, br)) + 
  scale_fill_manual(values = c(bl, g, br)) + 
  theme_pubr() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))

ggarrange(plot1, plot2, nrow = 1, common.legend = T, 
          labels = c("A", "B"),
          font.label = list(size = 30),
          hjust = 0, 
          widths = c(1, 0.5))
