


     ####################
     ### OTU Analyses ###
     ####################



# SET UP EVERYTHING

library(Biostrings)
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(scales)
library(vegan)
library(tidyverse)
library(RDPutils)
library(ggforce)
library(ggrepel)
library(microViz)
library(patchwork)
library(extrafont)
library(grid)
library(readr)

loadfonts(device = 'win')

setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/UNOISE3_alpha_5_output")

metaData <- read.csv('../0_metadata/Metadata_1719.csv', 
                     header = TRUE, 
                     check.names = FALSE, 
                     fileEncoding = "latin1")
freqTable <- read.table('zotutable_output_a5.txt', header = TRUE, sep = '\t', 
                           row.names = 1, check.names = FALSE, comment.char = '')
taxonomyTable <- read.table('taxonomy_output_a5.txt', header = TRUE, sep = '\t', 
                            row.names = 1, check.names = FALSE, 
                            comment.char = '', fill = TRUE, blank.lines.skip = FALSE)
sequenceTable <- readDNAStringSet('zotus_output_a5.fasta')

rownames(metaData) <- metaData$name


# The code below are some filtering as I need to combine data from 2017 and only TPC and CI from
# 2019. We can ignore this code later


if(FALSE){
row_sum <- rowSums(freqTable)
updated_freqTable <- freqTable[row_sum != 0, ]
CI_TPC_OTU <- rownames(taxonomyTable)
all_OTU <- rownames(updated_freqTable)
unmatched_rows <- updated_freqTable %>% 
  filter(!all_OTU %in% CI_TPC_OTU)

write.table(unmatched_rows, 'unmatched_rows.txt',
            append = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)
}

library(MetBrewer)
library(RColorBrewer)
library(cowplot)

rev(met.brewer("Derain"))

y32 <- c(met.brewer("Redon"), met.brewer("Monet"))
c1 <- as.vector(met.brewer("Cross"))
ycol <- c(rep(c1[1],18),rep(c1[5],5),rep(c1[9],18))
brewerPlus <- distinct_palette()


## site_col <- c("#ca0e12", "#f6bd21", "#2aa7de", "#25377f")

site_col <- c("#df8d71", "#d8b847", "#75884b", "#5b859e")

### If the phylum is unassigned, it is currently showing up as NA. Therefore, we
### need to change it into "Unassigned" before plotting

metaData2017 <- filter(metaData, yearCollected == 2017)

taxonomyTable <- taxonomyTable[, 2:7]


phyloseq <- phyloseq(otu_table(as.matrix(freqTable), taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomyTable)),
                       sample_data(metaData2017),
                       refseq(sequenceTable))

ps <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

merged_ps <- merge_samples(ps,"eventID")

#write.table(merged_ps@tax_table, 
#           'zotutable_filtered.txt', append = FALSE, sep = '\t', dec = '.', 
#            row.names = TRUE, col.names = NA)

### refilling the metadata which was lost in the previous step

sample_data(merged_ps)$monthCollected <- metaData$monthCollected[1:12]
sample_data(merged_ps)$name <- metaData$name[1:12]
sample_data(merged_ps)$eventID <- metaData$eventID[1:12]
sample_data(merged_ps)$country <- metaData$country[1:12]
sample_data(merged_ps)$locality <- metaData$locality[1:12]
sample_data(merged_ps)$continentOcean <- metaData$continentOcean[1:12]
sample_data(merged_ps)$impact <- metaData$impact[1:12]



     #################################################
     ## Bar Charts for Proportions by size-fraction ##
     #################################################



ps_fraction <- ps

ps_join(ps_fraction, metaData2017)

### If the phylum is unassigned, it is currently showing up as NA. Therefore, we
### need to change it into "Unassigned" before plotting

taxonomyTable_updated <- taxonomyTable
taxonomyTable_updated$phylum[is.na(taxonomyTable_updated$phylum)] <- "Unassigned"
tax_table(ps_fraction) <- as.matrix(taxonomyTable_updated)

desired_order <- c("COI_100_HK7", "COI_100_HK8", "COI_100_HK9", 
                   "COI_100_HK4", "COI_100_HK5", "COI_100_HK6", 
                   "COI_100_HK1", "COI_100_HK2", "COI_100_HK3", 
                   "COI_100_HK10", "COI_100_HK11", "COI_100_HK12", 
                   "COI_500_HK7", "COI_500_HK8", "COI_500_HK9", 
                   "COI_500_HK5", "COI_500_HK6", 
                   "COI_500_HK1", "COI_500_HK2", "COI_500_HK3", 
                   "COI_500_HK10", "COI_500_HK11", "COI_500_HK12",
                   "COI_Sessile_HK7", "COI_Sessile_HK8", "COI_Sessile_HK9", 
                   "COI_Sessile_HK4", "COI_Sessile_HK5", "COI_Sessile_HK6", 
                   "COI_Sessile_HK1", "COI_Sessile_HK2", "COI_Sessile_HK3", 
                   "COI_Sessile_HK10", "COI_Sessile_HK11", "COI_Sessile_HK12")


## Main figure

taxa_order = c("Annelida", "Arthropoda", "Bacillariophyta", 
              "Bryozoa", "Chordata", "Cnidaria", "Echinodermata", "Entoprocta", 
              "Mollusca", "Phaeophyceae", "Platyhelminthes", "Porifera", 
              "Rhodophyta", "Unassigned")

q <- ps_fraction %>% 
  tax_fix() %>% 
  comp_barplot(tax_level = "phylum", 
               n_taxa = 15, 
               merge_other = TRUE, 
               tax_order = taxa_order, 
               bar_width = 0.9, 
               x = "abbreviationLabel", 
               bar_outline_colour = NA,
               label = "abbreviationLabel",
               sample_order = desired_order, 
               palette = y32) +
  coord_flip() +
  facet_grid(cols = vars(sizeFraction), 
             scales = "fixed") +
  labs(x = "Individual ARMS by Location", 
       y = "Relative Abundance") +
  guides(fill = guide_legend(title = "Phylum")) +
  theme(
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.title = element_text(size = 12, face = "bold")
  )
q

#ggsave("Proportions_by_size_fraction.jpg", q, dpi = 600, width = 15, height = 5)



     #########################
     ## community structure ##
     #########################



# changing the data from abundance to presence absence

ps_pre_abs <- merged_ps
ps_pre_abs <- prune_taxa(taxa_sums(ps_pre_abs) > 0, ps_pre_abs)
ps_pre_abs@otu_table[1:12, ][ps_pre_abs@otu_table[1:12, ] > 0] <- 1

### The code below plots all phylum including the unassigned
### However, this might not be the best way to represent the data since
### the differences are not obvious (and according to documentation in ggplot,
### this could result in a deceiving plot)

taxonomy_pre_abs <- as(otu_table(ps_pre_abs), "matrix")
taxonomy_pre_abs <- as.data.frame(t(taxonomy_pre_abs))
OTUsum <- as.data.frame(colSums(taxonomy_pre_abs))
colnames(OTUsum) <- "OTU_Count"

# The figure below is not very informative. Not used currently
if(FALSE){
  ggplot(OTUsum, aes(x = rownames(OTUsum), y = OTU_Count)) +
  geom_col() +
  scale_x_discrete() +
  panel_border() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OTUs Count per site") +
  xlab("ARMS")
}



     ######################
     ## Ordination Plots ##
     ######################



## Based on presence-absence data ##
ps_pre_abs.ord <- ordinate(ps_pre_abs, "PCoA", "jaccard")

p1 <- plot_ordination(ps_pre_abs, ps_pre_abs.ord, color = "locality")
p1 <- p1 +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(shape = "Month Collected", color = "Location")
p1

ggsave("PCoA_jaccard_pre_abs.png", p, dpi = 600, width = 7, height = 5)


## Based on counts data ##
merged_ps.ord <- ordinate(merged_ps, "PCoA", "bray", autotransform = TRUE)
p2 <- plot_ordination(merged_ps, merged_ps.ord, color = "locality")
p2 <- p2 +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(shape = "Month Collected", color = "Location")
p2

ggsave("PCoA_bray_count.png", p2, dpi = 600, width = 7, height = 5)



## Based on counts data, but facet by phylum 
q1 <- plot_ordination(merged_ps, merged_ps.ord, type = "OTU", 
                     color = "phylum")

q1 <- q1 +
  geom_point(size = 2) + 
  scale_color_manual(values = c(y32)) +
  geom_text_repel(aes(label = order), size = 2) +
  theme(text = element_text(family = "Georgia")) +
  theme_minimal_grid() +
  facet_wrap(~ phylum) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10))

q1

layout <- c(
  area(1, 1),
  area(1, 2),
  area(2, 1, 3, 2)
)
plot(layout)

p1 + p2 + q1 + plot_layout(design = layout)

ggsave("PCoA_bray_count_facet_order.jpg", q1, dpi = 600, width = 14, height = 10)



     ########################################
     ## Ordination Plots by Size Fractions ##
     ########################################



## Based on count data ##

## 106 nm
ps_106 <- ps_filter(ps_fraction, sizeFraction == "106 nm")

ps_106.ord <- ordinate(ps_106, "PCoA", "bray", autotransform = TRUE)

p_106_bray <- plot_ordination(ps_106, ps_106.ord, color = "locality")

p_106_bray <- p_106_bray +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
#  scale_y_continuous(limits = c(-0.5, 0.5)) +
#  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Bray-Curtis Distance, 106 μm Fraction")
p_106_bray


## 500 nm
ps_500 <- ps_filter(ps_fraction, sizeFraction == "500 nm")

ps_500.ord <- ordinate(ps_500, "PCoA", "bray", autotransform = TRUE)

p_500_bray <- plot_ordination(ps_500, ps_500.ord, color = "locality")

p_500_bray <- p_500_bray +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse(aes(filter = locality != "Che Lei Pai")) +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Bray-Curtis Distance, 500 μm Fraction") +
  theme(legend.position = "none")
p_500_bray


## Sessile
ps_sessile <- ps_filter(ps_fraction, sizeFraction == "Sessile")

ps_sessile.ord <- ordinate(ps_sessile, "PCoA", "bray", autotransform = TRUE)

p_sessile_bray <- plot_ordination(ps_sessile, ps_sessile.ord, color = "locality")

p_sessile_bray <- p_sessile_bray +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Bray-Curtis Distance, Sessile Fraction")
p_sessile_bray


####################################
## Based on presence-absence data ##
####################################


ps_count <- tax_transform(ps, trans = "binary")
ps_pa <- ps_get(ps_count)

## 106 nm
ps_106_pa <- ps_filter(ps_pa, sizeFraction == "106 nm")

ps_106_pa.ord <- ordinate(ps_106_pa, "PCoA", "jaccard")

p_106_jaccard <- plot_ordination(ps_106_pa, ps_106_pa.ord, color = "locality")

p_106_jaccard <- p_106_jaccard +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Jaccard Distance, 106 μm Fraction")
p_106_jaccard


## 500 nm
ps_500_pa <- ps_filter(ps_pa, sizeFraction == "500 nm")

ps_500_pa.ord <- ordinate(ps_500_pa, "PCoA", "jaccard")

p_500_jaccard <- plot_ordination(ps_500_pa, ps_500_pa.ord, color = "locality")

p_500_jaccard <- p_500_jaccard +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse(aes(filter = locality != "Che Lei Pai")) +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Jaccard Distance, 500 μm Fraction") +
  theme(legend.position = "none")
p_500_jaccard


## Sessile
ps_sessile_pa <- ps_filter(ps_pa, sizeFraction == "Sessile")

ps_sessile_pa.ord <- ordinate(ps_sessile_pa, "PCoA", "jaccard")

p_sessile_jaccard <- plot_ordination(ps_sessile_pa, ps_sessile_pa.ord, color = "locality")

p_sessile_jaccard <- p_sessile_jaccard +
  geom_point(size = 3) + 
  scale_color_manual(values = site_col) +
  theme_light() +
  ggforce::geom_mark_ellipse() +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(color = "Location") +
  ggtitle("Jaccard Distance, Sessile Fraction")
p_sessile_jaccard



t1 <- p_500_bray + p_500_jaccard
t2 <- p_106_bray + p_106_jaccard
t3 <- p_sessile_bray +p_sessile_jaccard

t2 / t1 / t3 + plot_layout(guides = "collect")


bray_fraction <- p_106_bray / p_500_bray / p_sessile_bray
bray_fraction


jaccard_fraction <- p_106_jaccard / p_500_jaccard / p_sessile_jaccard 
jaccard_fraction


pcoa_fraction <- (bray_fraction | jaccard_fraction) + plot_layout(guides = "collect")
pcoa_fraction

ggsave("PCoA_size_fraction.jpg", pcoa_fraction, width = 10, height = 12, dpi = 600)



     ##################################################
     ## Counts by Proportions (Including Unassigned) ##
     ##################################################



### If the phylum is unassigned, it is currently showing up as NA. Therefore, we
### need to change it into "Unassigned" before plotting

taxonomyTable_updated <- taxonomyTable
taxonomyTable_updated$phylum[is.na(taxonomyTable_updated$phylum)] <- "Unassigned"
ps_pre_abs_phylum <- ps_pre_abs
tax_table(ps_pre_abs_phylum) <- as.matrix(taxonomyTable_updated)

Phyfac <- factor(tax_table(ps_pre_abs_phylum)[, "phylum"])
OTUtab_unique <-  apply(otu_table(ps_pre_abs_phylum), MARGIN = 1, function(x) {
  tapply(x, INDEX = Phyfac, FUN = sum, na.rm = F, simplify = TRUE)
})

OTUtab_unique <- as.data.frame(as.table(OTUtab_unique))

### Create another dataframe which only contains the assigned phylum for plotting

OTUtab_unique_assigned <- OTUtab_unique[OTUtab_unique$Var1 != "Unassigned", ]




if(FALSE){
### OTU proportions per phylum

ARMS_labels <- c("2017 OCT - CI", "2017 OCT - CI", "2017 OCT - CI", 
                 "2017 OCT - TPC", "2017 OCT - TPC", "2017 OCT - TPC", 
                 "2019 JAN - CI", "2019 JAN - CI", "2019 JAN - CI", 
                 "2019 JAN - TPC", "2019 JAN - TPC", "2019 JAN - TPC", 
                 "2019 AUG - CI", "2019 AUG - CI", "2019 AUG - CI", 
                 "2019 AUG - TPC", "2019 AUG - TPC", "2019 AUG - TPC")


f <- ggplot(OTUtab_unique, aes(x=Var2, y=Freq, fill=Var1, order=as.factor(Var1)))+
  geom_col(position = "fill")+
  scale_x_discrete()+
  panel_border() +
  scale_fill_manual(values=y32) +
  scale_x_discrete(label = ARMS_labels) +
  scale_y_continuous(expand=expansion(c(0,0)), labels = label_number(accuracy = 0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = ycol))+
  ylab("ZOTU Proportion per Phylum")+
  xlab("ARMS by sites and dates")+
  labs(fill = "Phylum")

f
}


### OTU proportions per phylum for 2017

ARMS_labels <- c("CI-1", "CI-2", "CI-3",
                 "CLP-1", "CLP-2", "CLP-3",
                 "PI-1", "PI-2", "PI-3",  
                 "TPC-1", "TPC-2", "TPC-3")

desired_order <- c("HKARMS_07", "HKARMS_08", "HKARMS_09", 
                   "HKARMS_04", "HKARMS_05", "HKARMS_06", 
                   "HKARMS_01", "HKARMS_02", "HKARMS_03", 
                   "HKARMS_10", "HKARMS_11", "HKARMS_12")

OTUtab_unique_2017 <- OTUtab_unique
OTUtab_unique_2017$Var2 <- factor(OTUtab_unique_2017$Var2, levels = desired_order)


f <- ggplot(OTUtab_unique_2017, aes(x=Var2, y=Freq, fill=Var1, order=as.factor(Var1)))+
  geom_col(position = "fill")+
  scale_x_discrete()+
  panel_border() +
  scale_fill_manual(values = y32[1:14]) +
  scale_x_discrete(label = ARMS_labels) +
  scale_y_continuous(expand=expansion(c(0,0)), labels = label_number(accuracy = 0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none"
        ) +
  ylab("Proportion per Phylum") +
  xlab("ARMS by sites") +
#  ggtitle("ZOTU Proportions per Site") +
  labs(fill = "Phylum")

f

ggsave("OTU_proportion.png", f, dpi = 600, device = png, width = 8, height = 6)




### Plot counts for assigned data

ggplot(OTUtab_unique_assigned, aes(x=Var2, y=Freq, fill=Var1, order=as.factor(Var1)))+
  geom_col()+
  scale_x_discrete()+
  panel_border() +
  scale_fill_manual(values=y32) +
  scale_x_discrete(label = ARMS_labels) +
  scale_y_continuous(expand=expansion(c(0,0)), labels = label_number(accuracy = 0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = ycol))+
  ylab("OTU Proportion per Phylum")+
  xlab("ARMS by sites and dates")+
  labs(fill = "Phylum")



     ##########################
     ## Excluding Unassigned ##
     ##########################



Phyfac <- factor(tax_table(ps_pre_abs)[, "phylum"])
OTUtab_unique <-  apply(otu_table(ps_pre_abs), MARGIN = 1, function(x) {
  tapply(x, INDEX = Phyfac, FUN = sum, na.rm = F, simplify = TRUE)
})

OTUtab_unique <- as.data.frame(as.table(OTUtab_unique))

### OTU proportions per phylum

ggplot(OTUtab_unique, aes(x=Var2, y=Freq, fill=Var1, order=as.factor(Var1)))+
  geom_col(position = "fill")+
  scale_x_discrete()+
  panel_border() +
  scale_fill_manual(values=y32) +
  scale_x_discrete(label = ARMS_labels) +
  scale_y_continuous(expand=expansion(c(0,0)), labels = label_number(accuracy = 0.1)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = ycol))+
  ylab("Proportion of OTUs per Phylum")+
  xlab("ARMS")+
  labs(fill = "Phylum")



     ###################
     ### ZOTU Counts ###
     ###################



OTUsum_location_year <- rownames_to_column(OTUsum, var = "Sample")
OTUsum_location_year <- OTUsum_location_year %>% 
  mutate(locality = sample_data(ps_pre_abs)$locality,
         yearCollected = sample_data(ps_pre_abs)$yearCollected,
         monthCollected = sample_data(ps_pre_abs)$monthCollected)



# Create the plot with improved aesthetics
c <- ggplot(OTUsum_location_year, 
            aes(x = as.numeric(factor(locality)), y = OTU_Count)) +
  geom_point(aes(fill = ifelse(Sample == "HKARMS_07", "grey", locality)),
             size = 4, shape = 21, color = "white") +
  
# Color scheme  
  scale_fill_manual(values = c("#df8d71",  "#d8b847", "grey", "#75884b","#5b859e")) +
  
# Polynomial smooth line for data excluding 'HKARMS_07'
  geom_smooth(data = OTUsum_location_year[OTUsum_location_year$Sample != 'HKARMS_07', ],
              aes(x = as.numeric(factor(locality)), y = OTU_Count),
              method = "lm", formula = y ~ poly(x, 2),
              se = TRUE, color = "dimgrey", linetype = "dashed", size = 1.5, alpha = 0.1, 
              level = 0.95) +

  # Customizing x-axis labels
  scale_x_continuous(breaks = 1:length(unique(OTUsum_location_year$locality)), 
                     labels = c("Center Island", "Che Lei Pai", "Port Island", "Tung Ping Chau"), 
                     expand = expansion(mult = c(0.1, 0.1))) +

  
  
# Theme and aesthetics  
  theme_light() +  
  ylab("ZOTU Count")+
  xlab("Location")+
  theme(
#    text = element_text(family = "Calibri"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
#    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
c

ggsave("OTU_counts.png", c, dpi = 600, width = 8, height = 6)




### Stitch the figures together

q <- q + ggtitle("(a)")
f <- f + ggtitle("(b)")
c <- c + ggtitle("(c)")

community_structure_figure <- q / (f | c) + plot_layout(guides = "collect")
community_structure_figure

ggsave("community_structure.png", community_structure_figure, 
       dpi = 600, width = 12, height = 8)



     ##############################################
     ## ZOTU counts, but with size fraction data ##
     ##############################################



ps_fraction_pre_abs <- tax_transform(ps_fraction, trans = "binary")

pre_abs_zotutable <- as(otu_table(ps_fraction_pre_abs), "matrix")

pre_abs_zotutable <- as.data.frame(t(pre_abs_zotutable))

pre_abs_zotutable <- as.data.frame(rowSums(pre_abs_zotutable))

colnames(pre_abs_zotutable) <- "OTU_Count"


zotu_fraction_sum <- rownames_to_column(pre_abs_zotutable, var = "Sample")
zotu_fraction_sum <- zotu_fraction_sum %>% 
  mutate(locality = sample_data(ps_fraction_pre_abs)$locality)

###
{zotu_fraction_sum <- zotu_fraction_sum %>%
  filter(Sample != "COI_500_HK7")}
###


k <- ggplot(zotu_fraction_sum, aes(x = locality, y = OTU_Count, fill = factor(locality), color = factor(locality))) +
  # Transparent violin plot with border
  geom_violin(fill = NA, linewidth = 1, alpha = 0.3) +
  # Adding the data points
  geom_point(size = 4.5, shape = 21, color = "white", stroke = 1.5) +
  # Adding a smooth line (linear model)
  geom_smooth(aes(group = 1), method = "loess", color = "blue", se = FALSE, size = 1.2, linetype = "dashed") +
  # Custom color palette
  scale_fill_manual(values = site_col[1:4]) +
  scale_color_manual(values = site_col[1:4]) +
  # Light theme with customizations
  theme_light() +  
  # Custom labels
  ylab("ZOTU Count") +
  xlab("Location") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("size_fraction_zotu_count_removed.jpeg", k, 
       dpi = 600, width = 12, height = 8)

k



     ###############################
     ## 500 fraction outlier test ##
     ###############################



library(outliers)


setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/UNOISE3_alpha_5_output")
freqTable <- read.table('zotutable_a5.txt', header = TRUE, sep = '\t', 
                        row.names = 1, check.names = FALSE, comment.char = '')
taxonomyTable <- read.table('sintaxLineageCRABS.txt', header = TRUE, sep = '\t', 
                            row.names = 1, check.names = FALSE, 
                            comment.char = '', fill = TRUE, blank.lines.skip = FALSE)
physeq <- phyloseq(otu_table(freqTable, taxa_are_rows = TRUE), 
                   sample_data(metaData),
                   tax_table(as.matrix(taxonomyTable)))

physeq

ps_500_unfiltered <- ps_filter(physeq, sizeFraction == "500 nm", yearCollected == 2017)
ps_500_unfiltered <- tax_transform(ps_500_unfiltered, trans = "binary")
ps_500_count <- ps_get(ps_500_unfiltered)

zotu_500 <- as(otu_table(ps_500_count), "matrix")

zotu_500 <- as.data.frame(t(zotu_500))

zotu_500_sum <- as.data.frame(rowSums(zotu_500))

colnames(zotu_500_sum) <- "OTU_Count"



### boxplot test
zotu_counts <- zotu_500_sum[[1]]  # Extract the first column

boxplot(zotu_counts, main = "Boxplot of ZOTU counts")

outliers <- boxplot.stats(zotu_counts)$out  # This extracts the outliers
print(outliers)



# Perform Grubbs' test for outliers
grubbs.test(zotu_counts)

outlier_value <- grubbs.test(zotu_counts)$alternative  # Shows which value is an outlier
print(outlier_value)



# Perform Dixon's test for outliers
outlier_test <- dixon.test(zotu_500_sum$OTU_Count)
outlier_test



# ggplot
ggplot(zotu_500_sum, aes(x = "", y = OTU_Count)) +
  geom_boxplot() +
  geom_point()

ggplot(zotu_500_sum, aes(x = "", y = OTU_Count)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 4, 
               coef = 1) +
  labs(y = "ZOTU Counts", title = "Boxplot of ZOTU Counts with Custom Outlier Detection") +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal()



     ########################
     ## Summary statistics ##
     ########################



mean_ZOTU_count <- OTUsum_location_year[-7, ] %>% 
  group_by(locality) %>% 
  summarize(mean_ZOTU = mean(OTU_Count))
mean_ZOTU_count



     ###############
     ## PERMANOVA ##
     ###############


# Bray-Curtis distance
ARMS_dist <- phyloseq::distance(merged_ps, method = "bray")
sample_info <- data.frame(sample_data(merged_ps))
locality_adonis2 <- adonis2(ARMS_dist ~ locality,  
                            data = sample_info,
                            permutations = 999, 
                            method = "bray")
locality_adonis2


# Jaccard Distance
ARMS_dist <- phyloseq::distance(ps_pre_abs, method = "jaccard")
sample_info <- data.frame(sample_data(ps_pre_abs))
locality_adonis2 <- adonis2(ARMS_dist ~ locality,  
                            data = sample_info,
                            permutations = 999, 
                            method = "jaccard")
locality_adonis2


# Bray-Curtis distance using different size fractions
ARMS_dist <- phyloseq::distance(ps_sessile, method = "bray")
sample_info <- data.frame(sample_data(ps_sessile))
locality_adonis2 <- adonis2(ARMS_dist ~ locality,  
                            data = sample_info,
                            permutations = 999, 
                            method = "bray")
locality_adonis2


# Jaccard distance using different size fractions
ARMS_dist <- phyloseq::distance(ps_sessile_pa, method = "jaccard")
sample_info <- data.frame(sample_data(ps_sessile_pa))
locality_adonis2 <- adonis2(ARMS_dist ~ locality,  
                            data = sample_info,
                            permutations = 999, 
                            method = "jaccard")
locality_adonis2



     ###########
     ## ANOVA ##
     ###########



library(car)
library(outliers)
library(report)
library(multcomp)

## Test for normality, one assumption of ANOVA
shapiro.test(OTUsum_location_year$OTU_Count)

## Test for equality of variances, another assumption of ANOVA
leveneTest(ZOTU_count$OTU_Count ~ ZOTU_count$locality)

## Test for outlier in Center Island, which suggests it is NOT an outlier
CI_counts <- OTUsum_location_year[7:9, 1:2]
CI_outlier_test <- dixon.test(CI_counts$OTU_Count)
CI_outlier_test

## Looks like there are significant differences between sites if removing the 
## "outlier" in Center Island

ZOTU_count <- OTUsum_location_year[-7, ]
ZOTU_count$locality <- as.factor(ZOTU_count$locality)
aov_count <- aov(OTU_Count ~ locality, data = ZOTU_count)
summary(aov_count)
report(aov_count)


## Post-hoc with Tukey HSD test
posthoc_count <- TukeyHSD(aov_count)
posthoc_count
plot(posthoc_count)



## Post-hoc Tukey HSD test
posthoc_count <- glht(aov_count,
                      linfct = mcp(locality = "Tukey")
                      )
summary(posthoc_count)
plot(posthoc_count)


## Post-hoc with Dunnett's test
#ZOTU_count$locality <- relevel(ZOTU_count$locality, ref = "Center Island")
#levels(ZOTU_count$locality)
aov_count <- aov(OTU_Count ~ locality, data = ZOTU_count)
posthoc_count <- glht(aov_count,
                      linfct = mcp(locality = "Dunnett")
)
summary(posthoc_count)
plot(posthoc_count)



     ##################
     ## Venn Diagram ##
     ##################



library(MiscMetabar)

venn_diagram <- ggvenn_pq(ps_pre_abs, fact = "locality") +
  ggplot2::scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(title = "Shared number of ZOTUs")
venn_diagram

ggsave("venn_diagram.jpeg", venn_diagram, 
       dpi = 600, width = 12, height = 8)



     ############################
     ## Figure 5 with MicroViz ##
     ############################



setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/UNOISE3_alpha_5_output")

taxonomyTable_order <- read.table('taxonomy_order_and_below_a5.txt', header = TRUE, sep = '\t', 
                            row.names = 1, check.names = FALSE, 
                            comment.char = '', fill = TRUE, blank.lines.skip = FALSE)

merged_ps_order <- phyloseq(otu_table(merged_ps),
                            sample_data(merged_ps), 
                            tax_table(as.matrix(taxonomyTable_order)))

ps_microviz <- merged_ps %>% 
  tax_fix() %>% 
  phyloseq_validate()

tax_agg(ps_microviz, "order") %>% 
  ps_get() %>% 
  tax_table()

site_color_heatmap <- c(
  "Centre Island" = "#df8d71", 
  "Che Lei Pai" = "#d8b847", 
  "Port Island" = "#75884b", 
  "Tung Ping Chau" = "#5b859e"
)


ps_microviz <- merged_ps_order %>% 
  tax_fix() %>% 
  phyloseq_validate()

taxa_heatmap <- ps_microviz %>%
  tax_transform(trans = "hellinger") %>%
  tax_filter(use_counts = TRUE) %>%
  comp_heatmap(
    colors = heat_palette(palette = "Rocket", rev = TRUE), grid_col = NA,
    sample_side = "top", name = "Abd.",
    show_row_names = FALSE, 
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))),
    sample_anno = sampleAnnotation(
      Location = anno_sample("locality"),
      col = list(Location = site_color_heatmap), border = FALSE))
taxa_heatmap

png("taxa_heatmap.png", width = 8, height = 10, units = "in", res = 600)
taxa_heatmap
dev.off()


## adding extra columns to help with making the legend

ps_microviz <- ps_microviz %>%
  ps_mutate(
    CI = if_else(locality == "Centre Island", true = 1L, false = 0L),
    CLP = if_else(locality == "Che Lei Pai", true = 1L, false = 0L),
    PI = if_else(locality == "Port Island", true = 1L, false = 0L),
    TPC = if_else(locality == "Tung Ping Chau", true = 1L, false = 0L)
  ) 

##
## Correlation heatmap between taxa and site, but with p-value ##
##

## compute correlations, with p values, and store in a dataframe
correlations_df <- ps_microviz %>% 
  tax_transform(trans = "log2", zero_replace = "halfmin") %>% 
  tax_model(
    rank = "order", variables = list("CI", "CLP", 
                                     "PI", "TPC"), 
    type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE
  ) %>% 
  tax_models2stats(rank = "order")

# get nice looking ordering of correlation estimates using hclust
taxa_hclust <- correlations_df %>% 
  dplyr::select(term, taxon, estimate) %>% 
  tidyr::pivot_wider(
    id_cols = taxon, names_from = term, values_from = estimate
  ) %>% 
  tibble::column_to_rownames("taxon") %>% 
  as.matrix() %>% 
  stats::dist(method = "euclidean") %>% 
  hclust(method = "ward.D2") 

taxa_order <- taxa_hclust$labels[taxa_hclust$order]

cor_p_value_map <- correlations_df %>% 
  mutate(p.FDR = p.adjust(p.value, method = "fdr")) %>%
  ggplot(aes(x = term, y = taxon)) +
  geom_raster(aes(fill = estimate)) +
  geom_point(
    data = function (x) x %>% filter(p.value < 0.05),
    shape = "asterisk", size = 2
  ) +
  geom_point(
    data = function (x) x %>% filter(p.FDR < 0.05),
    shape = "circle", size = 3
  ) +
  scale_y_discrete(limits = taxa_order) +
  scale_fill_distiller(palette = "RdBu") +
  theme_minimal() +
  labs(
    x = "Location", y = "Order-Level Taxa", fill = "Spearman's\nRank\nCorrelation",
    caption = paste(
      "Asterisk indicates p < 0.05, not FDR adjusted",
      "Filled circle indicates FDR corrected p < 0.05", sep = "\n"
    ))

cor_p_value_map

ggsave("correlation_p_value_map.jpeg", cor_p_value_map, device = jpeg, 
       width = 8, height = 10, dpi = 600)

## Correlation heatmap between taxa and environmental conditions, but with p-value ##

## compute correlations, with p values, and store in a dataframe
correlations_df_env <- ps_microviz %>% 
  tax_transform(trans = "log2", zero_replace = "halfmin") %>% 
  tax_model(
    rank = "order", variables = list("CI", "CLP", 
                                     "PI", "TPC"), 
    type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE
  ) %>% 
  tax_models2stats(rank = "order")

correlations_df_env



     ####################################
     ## PCoA and environmental factors ##
     ####################################



##### Load the environmental data #####
setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/analyses/pcoa_envfit")
env_data <- read.csv('env_data_mean_size_fraction.csv', 
                     header = TRUE, 
                     check.names = FALSE, 
                     row.names = 1)


##### Perform ordination with PCoA #####
ZOTU_hel <- decostand(otu_table(ps), method = "hellinger")

ps_hel <- phyloseq(otu_table(ZOTU_hel, taxa_are_rows = TRUE), 
                   sample_data(ps))

ps_ZOTU.ord <- ordinate(ps_hel, "PCoA", "bray")

ps_ZOTU <- plot_ordination(ps, ps_ZOTU.ord, color = "locality")


##### envfit to test environmental data #####
ZOTU_all_fraction <- as.data.frame(t(otu_table(ps)))
bray_dist <- vegdist(decostand(ZOTU_all_fraction, "hellinger"), method = "bray")
ZOTU_pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

envfit_result <- envfit(ZOTU_pcoa, env_data, perm = 999)
arrows <- as.data.frame(scores(envfit_result, display = "vectors"))
arrows$variable <- rownames(arrows)


##### Plot the ordination with envfit #####
ps_ZOTU <- ps_ZOTU +
  geom_point(size = 3) + 
  geom_text_repel(aes(label = name), size = 3) +
  scale_color_manual(values = site_col) +
  theme_minimal_grid() +
  labs(color = "Location")
  #  ggforce::geom_mark_ellipse() +
  #  scale_y_continuous(limits = c(-0.5, 0.5)) +
  #  scale_x_continuous(limits = c(-0.5, 0.5)) +

ps_ZOTU <- ps_ZOTU + 
  geom_segment(data = arrows, 
               aes(x = 0, y = 0, 
                   xend = Dim1 * 0.5, 
                   yend = Dim2 * 0.5),  # Adjust scaling
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = arrows, 
                  aes(x = Dim1 * 0.5, y = Dim2 * 0.5, label = variable), 
                  color = "blue", 
                  size = 4, 
                  box.padding = 0.5,    # Increase padding between labels
                  point.padding = 0.5,  # Increase padding between label and arrow
                  segment.size = 0.2)
ps_ZOTU

ggsave("env_data_on_PCoA.jpeg", ps_ZOTU, device = jpeg, width = 12, height = 9, dpi = 600)



     #######################################
     ## Environmental Data Panel Plotting ##
     #######################################



library(tidyverse)
library(report)

setwd("C:/Baker_Lab/ARMS_data/COI_Tolo_Harbour/analyses/pcoa_envfit")
library(broom)

env_data_all <- read.csv('env_data_all.csv', 
                     header = TRUE, 
                     check.names = FALSE, 
                     sep = ",")

env_data_long <- env_data_all %>% 
  pivot_longer(cols = -c(Dates, Station, `Sample No`, Depth, `Water Control Zone`, Location), 
               names_to = "Parameter", 
               values_to = "Value")

env_data_long <- env_data_all %>%
  mutate(across(-c(Dates, Station, `Sample No`, Depth, `Water Control Zone`, Location), 
                as.character)) %>%  # Convert all pivot columns to character
  pivot_longer(cols = -c(Dates, Station, `Sample No`, Depth, `Water Control Zone`, Location), 
               names_to = "Parameter", 
               values_to = "Value")



site_col <- c(
  "CI" = "#df8d71", 
  "CLP" = "#d8b847", 
  "PI" = "#75884b", 
  "TPC" = "#5b859e"
  )

env_data_long$Station <- factor(env_data_long$Station, 
                                levels = c("MM5", "MM17", "TM7TM8", "TM4"),
                                labels = c("MM5(TPC)", "MM17(PI)", "TM7TM8(CLP)", "TM4(CI)"))

env_data_long$Station <- factor(env_data_long$Location)

### ANOVA and p-value
# Perform ANOVA for each parameter

env_data_long <- env_data_long %>%
  mutate(Dates = as.Date(Dates, format = "%m/%d/%Y"))  # Convert "Dates" to Date class

anova_results <- env_data_long %>%
  # Convert Value column to numeric, coercing invalid values to NA
  mutate(Value = as.numeric(Value)) %>%
  # Remove rows with NA in Value
  filter(!is.na(Value)) %>%
  group_by(Parameter) %>%
  do(tidy(aov(Value ~ Location, data = .))) %>%
  filter(term == "Location") %>%
  select(Parameter, p.value)

# Print ANOVA results to check
print(anova_results)

env_data_long <- env_data_long %>%
  left_join(anova_results, by = "Parameter")

env_data_long$Value <- as.numeric(env_data_long$Value)


# Plotting
e <- ggplot(env_data_long, aes(x = Dates, y = Value, color = Location, group = Location)) +
  geom_point(size = 1) +
  geom_line(linewidth = 1) +
  theme_minimal() +  # Use a clean, minimalistic theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),) +
  labs(
    x = NULL,  # Add axis labels
    y = NULL,
    color = "Station"
  ) +
  scale_color_manual(values = site_col) + 
  scale_y_continuous(n.breaks = 5) +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 2)


e <- e + 
  geom_text(data = subset(env_data_long, Parameter == "5dOD (mg/L)"),
                   aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
                   size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "Chlo (μg/L)"),
            aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
            size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "DO (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.199"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "FC (cfu/100mL)"),
            aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
            size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "OP (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.686"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "PhP (μg/L)"),
            aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
            size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "Sal (psu)"),
            aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
            size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "TIN (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.401"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "TKN (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.202"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "TN (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.452"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "Turb (NTU)"),
            aes(x = max(Dates), y = Inf, label = "p < 0.01"), 
            size = 4, color = "red", hjust = 0.7, vjust = 1.1, inherit.aes = FALSE) +
  geom_text(data = subset(env_data_long, Parameter == "VSS (mg/L)"),
            aes(x = max(Dates), y = Inf, label = "p = 0.638"), 
            size = 4, hjust = 0.8, vjust = 1.1, inherit.aes = FALSE)
  
e

ggsave("env_data_panel_plot.jpeg", e, device = jpeg, width = 9, height = 12, dpi = 600)



### ANOVA of each single parameter, followed up with post-hoc Tukey's Test

aov_env_data <- aov(`VSS (mg/L)` ~ Location, data = env_data_all)
summary(aov_env_data)
report(aov_env_data)


TukeyHSD(aov_env_data)
plot(TukeyHSD(aov_env_data))

