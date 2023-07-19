library("tidyverse")
library ("phyloseq")
library("vegan")

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq.RData")

sampledata_meta_2$Sample = row.names(sampledata_meta_2)
physeq_for_beta = psmelt(physeq_for_beta)

physeq_for_beta_Species = physeq_for_beta

physeq_for_beta_Species = physeq_for_beta_Species %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq_for_beta_Species = pivot_wider(physeq_for_beta_Species, id_cols = Genus, names_from = "Sample", values_from = "Abundance")
physeq_for_beta_Species = as.data.frame(physeq_for_beta_Species)
row.names(physeq_for_beta_Species) = physeq_for_beta_Species[,1]
physeq_for_beta_Species = physeq_for_beta_Species[,-1]
physeq_for_beta_Species = t(physeq_for_beta_Species)

vdist = vegdist(physeq_for_beta_Species, method="bray")

vdist_NMDS = metaMDS(physeq_for_beta_Species, distance = "bray", k = 2, maxit = 999,  trymax = 500, wascores = TRUE)
NMDS = as.data.frame(scores(vdist_NMDS, display="site"))
NMDS$Sample = rownames(NMDS)

NMDS = left_join(NMDS, sampledata_meta_2, by ="Sample")

mycolors = c("#FF1202", "#401CE6", "#11a40e")

# Figure 8d. NMDS ordination plot of Vitis communities (ITS)
Vitis = c("V. amurensis", "V. coignetiae", "V. vinifera cv. ‘Cabernet Dorsa’ ")

NMDS$Species = factor(NMDS$Species, levels = Vitis)
NMDS$Plant = factor(NMDS$Plant, levels = c(V_amur, V_coig, V_vinifera))

ggplot(data = NMDS, aes(x = NMDS1, y = NMDS2, color = Species, shape = Plant)) + 
	geom_point(size=5) +
	scale_colour_manual(values = mycolors) +
	scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8)) +
	theme_bw() +
	guides(shape = guide_legend(order = 1), 
		color = guide_legend(order = 2)) +
	theme(legend.text = element_text(size = 15, face = "italic"),
		legend.title = element_text(size = 20)
	)
ggsave("./Figure 8d. NMDS ordination plot of Vitis communities (ITS).png", width = 11, height = 10)
#


#####################


physeq_for_beta_Plant = physeq_for_beta

physeq_for_beta_Plant = physeq_for_beta_Plant %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq_for_beta_Plant = left_join(physeq_for_beta_Plant, sampledata_meta_2, by="Sample")
physeq_for_beta_Plant = filter(physeq_for_beta_Plant, Plant != "V. vinifera-USA" & Plant != "V. vinifera-Germany")
physeq_for_beta_Plant = pivot_wider(physeq_for_beta_Plant, id_cols = Genus, names_from = "Sample", values_from = "Abundance")
physeq_for_beta_Plant = as.data.frame(physeq_for_beta_Plant)
row.names(physeq_for_beta_Plant) = physeq_for_beta_Plant[,1]
physeq_for_beta_Plant = physeq_for_beta_Plant[,-1]
physeq_for_beta_Plant = t(physeq_for_beta_Plant)

vdist2 = vegdist(physeq_for_beta_Plant, method="bray")

vdist_NMDS2 = metaMDS(physeq_for_beta_Plant, distance = "bray", k = 2, maxit = 999,  trymax = 500, wascores = TRUE)
NMDS2 = as.data.frame(scores(vdist_NMDS2, display="site"))
NMDS2$Sample = rownames(NMDS2)

NMDS2 = left_join(NMDS2, sampledata_meta_2, by ="Sample")

# Supporting information 8. NMDS ordination plot of V. amurensis and V. coignetiae plant samples (ITS)
NMDS2$Species = factor(NMDS2$Species, levels = Vitis)
NMDS2$Plant = factor(NMDS2$Plant, levels = c(V_amur, V_coig, V_vinifera))

ggplot(data = NMDS2, aes(x = NMDS1, y = NMDS2, color = Species, shape = Plant)) + 
	geom_point(size=5) +
	scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8)) +
	scale_colour_manual(values = mycolors) +
	theme_bw() +
		guides(shape = guide_legend(order = 1), 
		color = guide_legend(order = 2)) +
	theme(legend.text = element_text(size = 15, face = "italic"),
		legend.title = element_text(size = 20)
	)
ggsave("./Supporting information 8. NMDS ordination plot of V. amurensis and V. coignetiae plant samples (ITS).png", width = 11, height = 10)
#


# PERMANOVA test of V. amurensis and V. coignetiae plant's microbial communities (ITS)
beta_Plant = adonis2(physeq_for_beta_Plant ~ Plant, data = NMDS2, permutations = 999, method = "bray")
print("PERMANOVA test of V. amurensis and V. coignetiae plant's microbial communities (ITS)")
beta_Plant
write.table(beta_Plant, sep = "\t", quote = FALSE, file="./PERMANOVA test of V. amurensis and V. coignetiae plant's microbial communities (ITS).txt")

# PERMANOVA test of V. amurensis, V. coignetiae and V. vinifera cv. ‘Cabernet Dorsa’ microbial communities (ITS)
beta_Species = adonis2(physeq_for_beta_Species ~ Species, data = NMDS, permutations = 999, method = "bray")
print("PERMANOVA test of V. amurensis, V. coignetiae and V. vinifera cv. ‘Cabernet Dorsa’ microbial communities (ITS)")
beta_Species
write.table(beta_Species, sep = "\t", quote = FALSE, file="./PERMANOVA test of V. amurensis, V. coignetiae and V. vinifera cv. ‘Cabernet Dorsa’ microbial communities (ITS).txt")
