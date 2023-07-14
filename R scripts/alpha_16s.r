library("tidyverse")
library ("phyloseq")
library("vegan")


load("physeq.RData")


sampledata_meta_2[sampledata_meta_2[] == "V. vinifera cv. ‘Cabernet Dorsa’ "] = "V. vinifera-USA"
sampledata_meta_2[sampledata_meta_2[] == "V. vinifera cv. ‘Syrah’"] = "V. vinifera-Germany"

sampledata_meta_2$Sample = rownames(sampledata_meta_2)

# Figure 7c. Shannon′s alpha diversity boxplot (16s)
alpha_diversity = estimate_richness(physeq_for_alfa, split = TRUE, measures = NULL)
alpha_diversity$Sample = rownames(alpha_diversity)
alpha_diversity = left_join(sampledata_meta_2, alpha_diversity, by = "Sample")

Vitis = c("V. amurensis", "V. coignetiae", "V. vinifera-USA", "V. vinifera-Germany")

ggplot(data = alpha_diversity, aes(x = Species, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	scale_x_discrete(limits = Vitis) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("../Figure 7c. Shannon′s alpha diversity boxplot (16s).png", width = 7, height = 10)
#


# Pairwise Wilcoxon rank sum test of Vitis microbial communities (16s)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Species, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test of Vitis microbial communities (16s)")
Shannon