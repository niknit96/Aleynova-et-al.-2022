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

# Supporting information 9. Heatmap for V. amurensis plants (ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Supporting information 9. Heatmap for V. amurensis plants (ITS).png",  width = 20, height = 15, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "Top 10 most represented taxa of genus level for each V. amurensis plant",
	column_title_gp = gpar(fontsize = 20, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = V_amur,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 20),
	column_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_rot = 0,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 10))
})

draw(Hmap)

dev.off()
#

##########

upset_data = pre_upset_data[, names(pre_upset_data) %in% V_coig]
V_coig = V_coig[V_coig %in% colnames(upset_data)]

# Supporting information 11. Heatmap for V. coignetiae plants (ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Supporting information 11. Heatmap for V. coignetiae plants (ITS).png",  width = 10, height = 9, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "Top 10 most represented taxa of genus level \nfor each V. coignetiae plant",
	column_title_gp = gpar(fontsize = 15, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = V_coig,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 13),
	column_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_rot = 0,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 15), 
		title_gp = gpar(fontsize = 15)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 10))
})

draw(Hmap)

dev.off()
#


##########



upset_data = pre_upset_data

upset_data['V. amurensis'] = rowMeans(upset_data[,names(upset_data) %in% V_amur])
upset_data['V. coignetiae'] = rowMeans(upset_data[,names(upset_data) %in% V_coig])

Vitis = c("V. amurensis", "V. coignetiae", "V. vinifera-USA")

upset_data = upset_data[, Vitis]

# Supporting information 14. Heatmap for Vitis plants (ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Supporting information 14. Heatmap for Vitis plants (ITS).png",  width = 10, height = 9, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "Top 10 most represented taxa of genus level \nfor Vitis plants",
	column_title_gp = gpar(fontsize = 20, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = Vitis,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 13),
	column_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_rot = 75,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 15), 
		title_gp = gpar(fontsize = 15)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 10))
})

draw(Hmap)

dev.off()
#

##########

upset_data = pre_upset_data

upset_data['V. amurensis'] = rowMeans(upset_data[,names(upset_data) %in% V_amur])
upset_data['V. coignetiae'] = rowMeans(upset_data[,names(upset_data) %in% V_coig])

upset_data = upset_data[, c("V. amurensis", "V. coignetiae")]

# Figure 6. V. amurensis and V. coignetiae heatmap (ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 20) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Figure 6. V. amurensis and V. coignetiae heatmap (ITS).png",  width = 20, height = 15, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	cluster_columns = FALSE,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 25),
	column_names_gp = gpar(fontsize = 40, fontface = "italic"),
	column_names_rot = 0,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 20))
})

draw(Hmap)

dev.off()
#

##########


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

# Supporting information 10. Heatmap for V. amurensis plants (culture-depend, ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 1000000) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Supporting information 10. Heatmap for V. amurensis plants (culture-depend, ITS).png",  width = 15, height = 10, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "Taxa of genus level for each V. amurensis plant",
	column_title_gp = gpar(fontsize = 20, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = V_amur,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 20),
	column_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_rot = 0,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 10))
})

draw(Hmap)

dev.off()
#


#######


upset_data = pre_upset_data[, names(pre_upset_data) %in% V_coig]
V_coig = V_coig[V_coig %in% colnames(upset_data)]

upset_data[is.na(upset_data)] = 0

# Supporting information 12. Heatmap for V. coignetiae plants (culture-depend, ITS)
Heatmap = upset_data
Heatmap = filter(Heatmap, rowSums(Heatmap) > 0)

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 1000000) %>%
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Supporting information 12. Heatmap for V. coignetiae plants (culture-depend, ITS).png",  width = 8, height = 5, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "Taxa of genus level \nfor each V. coignetiae plant",
	column_title_gp = gpar(fontsize = 20, fontface = "italic"),
	cluster_columns = FALSE,
	column_order = V_coig,
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 13),
	column_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_rot = 0,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 15), 
		title_gp = gpar(fontsize = 15)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 10))
})

draw(Hmap)

dev.off()
#
