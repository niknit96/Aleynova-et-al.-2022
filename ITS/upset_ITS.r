library ("RColorBrewer")
library("tidyverse")
library ("phyloseq")
library("ComplexHeatmap")
library("circlize")


dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq.RData")


physeq_main = psmelt(physeq_main)

pre_upset_data = pivot_wider(physeq_main, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data[, names(pre_upset_data) %in% V_amur]
V_amur = V_amur[V_amur %in% colnames(upset_data)]

upset_data[] = ifelse(upset_data[] > 0, 1, 0) # Transform data to upset format
upset_data = upset_data[!apply(upset_data == 0, 1, all),]

# A table that shows which plants contain taxa from Upset plots (for V. amurensis)
upset_data_for_suppl = upset_data[, V_amur]
upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)
for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]
upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}
Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="./Supporting information 9. Intersections in ITS metagenome data of V. amurensis plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 4b. Upset diagram for V. amurensis plants (ITS)
png("./Figure 4b. Upset diagram for V. amurensis plants (ITS).png",  width = 30, height = 15, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 35),
	set_order = V_amur, 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()


#######


upset_data = pre_upset_data[, names(pre_upset_data) %in% V_coig]
V_coig = V_coig[V_coig %in% colnames(upset_data)]

upset_data[] = ifelse(upset_data[] > 0, 1, 0)
upset_data = upset_data[!apply(upset_data == 0, 1, all),]

# A table that shows which plants contain taxa from Upset plots (for V. coignetiae)
upset_data_for_suppl = upset_data[, V_coig]

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)

for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
# upset_data_for_suppl = sapply(upset_data_for_suppl, paste)
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]


upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}

Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="./Supporting information 11. Intersections in ITS metagenome data of V. coignetiae plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

png("./Figure 5b. Upset diagram for V. coignetiae plants (ITS).png",  width = 10, height = 7, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 35),
	column_title = "", 
	set_order = V_coig, 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()


##########


upset_data = pre_upset_data

upset_data['V. amurensis'] = rowMeans(upset_data[,names(upset_data) %in% V_amur])
upset_data['V. coignetiae'] = rowMeans(upset_data[,names(upset_data) %in% V_coig])

Vitis = c("V. amurensis", "V. coignetiae", "V. vinifera-USA")

upset_data = upset_data[, Vitis]

upset_data[] = ifelse(upset_data[] > 0, 1, 0)

# Supporting information 14. Intersections in ITS metagenome data of V. amurensis, V. coignetiae and V. vinifera plant samples
upset_data_for_suppl = upset_data[, Vitis]

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)

for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}

upset_data_for_suppl = t(upset_data_for_suppl)[-1,]

upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}

Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="./Supporting information 14. Intersections in ITS metagenome data of V. amurensis, V. coignetiae and V. vinifera plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 8b. Upset diagram for Vitis plants (ITS)
png("./Figure 8b. Upset diagram for Vitis plants (ITS).png",  width = 10, height = 7, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_title = "", 
	set_order = Vitis, 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()
#

#######################


physeq_sow = transform_sample_counts(physeq_sow, function(x) x / sum(x) * 100)

physeq = psmelt(physeq_sow) %>%
	filter(Abundance > 0)

pre_upset_data = pivot_wider(physeq, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")

pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data[, names(pre_upset_data) %in% V_amur]
V_amur = V_amur[V_amur %in% colnames(upset_data)]

upset_data[is.na(upset_data)] = 0
upset_data[] = ifelse(upset_data[] > 0, 1, 0)
upset_data = upset_data[!apply(upset_data == 0, 1, all),]


# Supporting information 10. Intersections in ITS data of cultivate-dependent method of V. amurensis plant samples
upset_data_for_suppl = upset_data[, V_amur]

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)

for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}

upset_data_for_suppl = t(upset_data_for_suppl)[-1,]

upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}

Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="./Supporting information 10. Intersections in ITS data of cultivate-dependent method of V. amurensis plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 4d. Upset diagram for V. amurensis plants (ITS, culture-depend)
png("./Figure 4d. Upset diagram for V. amurensis plants (ITS, culture-depend).png",  width = 12, height = 8, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 35),
	column_title = "", 
	set_order = V_amur, 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()

#######


upset_data = pre_upset_data[, names(pre_upset_data) %in% V_coig]
V_coig = V_coig[V_coig %in% colnames(upset_data)]

upset_data[is.na(upset_data)] = 0
upset_data[] = ifelse(upset_data[] > 0, 1, 0)
upset_data = upset_data[!apply(upset_data == 0, 1, all),]

# A table that shows which plants contain taxa from Upset plots (for V. coignetiae)
upset_data_for_suppl = upset_data[, V_coig]

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)

for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}

upset_data_for_suppl = t(upset_data_for_suppl)[-1,]

upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}

Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Taxa of genus level", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="./Supporting information 12. Intersections in ITS data of cultivate-dependent method of V. coignetiae plant samples.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 5d. Upset diagram for V. coignetiae plants (ITS, culture-depend)
png("./Figure 5d. Upset diagram for V. coignetiae plants (ITS, culture-depend).png",  width = 10, height = 7, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 35),
	column_title = "", 
	set_order = V_coig, 
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 18)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 15)),
		annotation_name_gp = gpar(fontsize = 15),
		numbers_gp = gpar(fontsize = 15)))
draw(upset_plot)
dev.off()
#