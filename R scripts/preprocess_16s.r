library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

asv_meta = read_qza("./dada2/FeatureTable[Frequency]_16s.qza") # Table of metagenome data reads for each ASV by samples
asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)

tax_meta = read_qza("feature-classifier_classify-sklearn/FeatureData[Taxonomy]_16s.qza") # Taxonomy table for metagenome data
tax_meta = parse_taxonomy(tax_meta$data, trim_extra=FALSE)
tax_meta[is.na(tax_meta)] <- "uncultured"

# Swap values ​​from Swap_unidentified in tax_meta to previous taxon level
Swap_unidentified = c("uncultured", "unidentified", "metagenome", "bacteriap25", "Unknown")
for(Swap in Swap_unidentified) {
    for(j in 1:7) {
        for (i in 1:nrow(tax_meta)) {
            if(grepl(Swap, tax_meta[i,j])) {
                tax_meta[i,j] = tax_meta[i,j-1] }
                else { tax_meta[i,j] = tax_meta[i,j]
            }       
        }
    }
}
#

# Filtering metagenome data from non-significant taxa
tax_meta = filter(tax_meta, 
Phylum != "d__Bacteria" & 
Phylum != "Unassigned" &
Kingdom != "d__Eukaryota" & 
Phylum != "k__Fungi" & 
Kingdom != "d__Archaea" &
Kingdom != "k__Viridiplantae" &
Genus != "g__Mitochondria" &
Genus != "g__Chloroplast")
#

tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")

sampledata_meta = read.table(file="./sampledata_16s_meta.txt", row.names = 1, header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome

V_amur = c("Gh","M","S-Va","P-1","P-2","P-3","P-4","P-5","P-6","Kh-1","Kh-2")
V_coig = c("S-1", "S-2", "S-3")
V_vinifera = c("V. vinifera-USA","V. vinifera-Germany")

sampledata_meta_2 = sampledata_meta

# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_meta)
tax_meta = tax_table(as.matrix(tax_meta))
asv_meta = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq = phyloseq(asv_meta, tax_meta, sampledata_meta)
#

physeq_read_counts = physeq
physeq_read_counts = merge_samples(physeq_read_counts, "Plant", fun = mean)

physeq_for_div = physeq

physeq = tax_glom(physeq, taxrank="Genus", NArm=TRUE) # Merge ASVs into genus level taxa

physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100) # Convert read counts to relative abundance (%)

physeq = merge_samples(physeq, "Plant", fun = mean) # Merge samples by plant
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)
physeq_for_info = physeq

physeq_main = filter_taxa(physeq, function(x) max(x) > 0.1, TRUE) # Genus level taxa with relative abundance > 0.1%

physeq_other = filter_taxa(physeq, function(x) max(x) < 0.1, TRUE) # Genus level taxa with relative abundance < 0.1%
physeq_other = merge_taxa(physeq_other, taxa_names(physeq_other))
taxa_names(physeq_other) = "Other"

# phyloseq object where genus level taxa with a relative abundance of less than 0.1% are placed in the "other" category
physeq = merge_phyloseq(physeq_main, physeq_other)
othertaxa = rbind(rank_names(tax_table(physeq)), "Other")
colnames(othertaxa) = othertaxa[1,]
othertaxa = othertaxa[-1,]
othertaxa = t(as.matrix(othertaxa))
row.names(othertaxa) = "Other"
tax_table(physeq)["Other",] <- tax_table(othertaxa)
#

# phyloseq object for beta diversity analysis
physeq_for_div = tax_glom(physeq_for_div, taxrank="Genus", NArm=TRUE)
physeq_for_beta = prune_taxa(taxa_names(physeq_main), physeq_for_div)
total = median(sample_sums(physeq_for_beta))
standf = function(x, t=total) round(t * (x / sum(x))) # for transformation for the even sample depth
physeq_for_beta = transform_sample_counts(physeq_for_beta, standf)
#

physeq_for_alfa = physeq_for_div # phyloseq object for alpha diversity analysis

####

tax_sow = read.table(file="./tax_16s_sow.txt", row.names = 1, header = TRUE, sep="\t") # Taxonomy table for cultivation-dependent approach data

sampledata_sow = read.table(file="./sampledata_16s_sow.txt", row.names = 1, header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the cultivation-dependent approach data

asv_sow = read.table(file="./asv_16s_sow.txt", header = TRUE, sep="\t") # Table of cultivation-dependent approach data reads for each strain by samples
asv_sow = asv_sow %>%
	group_by(., Species, Sample) %>%
	summarise(., Abundance = sum(Abundance))
asv_sow = pivot_wider(asv_sow, id_cols = "Species", names_from = "Sample", values_from = "Abundance")
asv_sow[is.na(asv_sow)] <- 0
asv_sow = as.data.frame(asv_sow)
row.names(asv_sow) = asv_sow[,1]
asv_sow = asv_sow[,-1]
# colnames(asv_sow) <- paste("sample", gsub("-", "_", colnames(asv_sow)), sep = "_")
asv_sow[] = apply(asv_sow[], 2, as.numeric)

# Creating a phyloseq object for cultivation-dependent approach data
sampledata_sow = sample_data(sampledata_sow)
tax_sow = tax_table(as.matrix(tax_sow))
asv_sow = otu_table(asv_sow, taxa_are_rows = TRUE)
physeq_sow = phyloseq(asv_sow, tax_sow, sampledata_sow)
physeq_sow = merge_samples(physeq_sow, "Plant", fun = sum) # Merge samples by plant
#


## End


physeq_info_before_filtration_V_amur = filter_taxa(prune_samples(V_amur, physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_V_coig = filter_taxa(prune_samples(V_coig, physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_V_USA = filter_taxa(prune_samples("V. vinifera-USA", physeq_for_info), function(x) max(x) > 0, TRUE)
physeq_info_before_filtration_V_Germany = filter_taxa(prune_samples("V. vinifera-Germany", physeq_for_info), function(x) max(x) > 0, TRUE)

physeq_info_after_filtration_V_amur = filter_taxa(prune_samples(V_amur, physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_V_coig = filter_taxa(prune_samples(V_coig, physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_V_USA = filter_taxa(prune_samples("V. vinifera-USA", physeq_main), function(x) max(x) > 0, TRUE)
physeq_info_after_filtration_V_Germany = filter_taxa(prune_samples("V. vinifera-Germany", physeq_main), function(x) max(x) > 0, TRUE)



print(paste("Taxa of genus level before filtration in V. amurensis samples:", length(taxa_names(physeq_info_before_filtration_V_amur))))
print(paste("Taxa of genus level before filtration in V. coignetiae samples:", length(taxa_names(physeq_info_before_filtration_V_coig))))
print(paste("Taxa of genus level before filtration in V. vinifera-USA samples:", length(taxa_names(physeq_info_after_filtration_V_USA))))
print(paste("Taxa of genus level before filtration in V. vinifera-Germany samples:", length(taxa_names(physeq_info_after_filtration_V_Germany))))

print(paste("Taxa of genus level with relative abundance > 0.1% in V. amurensis samples:", length(taxa_names(physeq_info_after_filtration_V_amur))))
print(paste("Taxa of genus level with relative abundance > 0.1% in V. coignetiae samples:", length(taxa_names(physeq_info_after_filtration_V_coig))))
print(paste("Taxa of genus level with relative abundance > 0.1% V. vinifera-USA samples:", length(taxa_names(physeq_info_after_filtration_V_USA))))
print(paste("Taxa of genus level with relative abundance > 0.1% V. vinifera-Germany samples:", length(taxa_names(physeq_info_after_filtration_V_Germany))))




save(physeq, physeq_for_div, V_amur, V_coig, V_vinifera, physeq_for_beta, physeq_for_alfa, physeq_main, sampledata_meta_2, sampledata_sow, physeq_read_counts, physeq_sow,  file = "physeq.RData")