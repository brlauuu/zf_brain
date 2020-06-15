library(shinycssloaders)
library(Seurat)
library(DT)
library(dplyr)
library(ggplot2)
library(URD)

seurat.version <<- packageVersion("Seurat")
urd.version <<- packageVersion("URD")

print(paste0("Seurat version being used: ", seurat.version))

f <- list.files("data")
path.to.load <- c("<select>", f[grepl(".*rds", f)])
clusters <- c("<select>")
genes <- c("<select>")

features.urd <- c("<select>")

units <- c("mm", "cm", "in")
device <- c("png", "eps", "ps", "tex", "pdf", "jpeg", "tiff", "bmp", "svg")

object <- NULL
object.urd <- NULL

loadObject <- function(path) {
	tryCatch({
		if (seurat.version > "2" && seurat.version < "3") {
			object <<- readRDS(paste0("data/", path))
		} else if (seurat.version > "3") {
			object <<- UpdateSeuratObject(readRDS(paste0("data/", path)))
		}
		print(paste0("Loaded Seurat object from path: ", "data/", path))
	}, error=function(cond) {
		print(paste("Invalid Seurat file loaded:", "data/", path))
		print("Here's the original error message:")
		print(cond)
		object <<- NULL
	})
}

loadObjectURD <- function(path) {
	object.urd <<- readRDS(paste0("data/", path))
	print(paste0("Loaded URD object from path: ", "data/", path))
}

listClusters <- function (path) {
	tryCatch({
		return(
			c(
				"all",
				meta.data %>%
					filter(STAGE == clustering.data$stage[clustering.data$path==path]) %>%
					select(CLUSTER.NAME)
			)
		)
	}, error=function(cond) {
		print(paste("Invalid Seurat file loaded:", "data/", path))
		print("Here's the original error message:")
		print(cond)
		return("<select>")
	})
}

listGenes <- function() {
	tryCatch({
		if (seurat.version > "2" && seurat.version < "3") {
			return(c("<select>", sort(rownames(object@raw.data))))
		} else if (seurat.version > "3") {
			return(c("<select>", sort(rownames(object))))
		}
	}, error=function(cond) {
		print(paste("Invalid Seurat file loaded:", "data/", path))
		print("Here's the original error message:")
		print(cond)
		return("<select>")
	})
}

listFeaturesURD <- function() {
	tryCatch({
		return(
			c(
				"pseudotime",
				colnames(object.urd@meta),
				rownames(object.urd@logupx.data)
			)
		)
	}, error=function(cond) {
		print(paste("Invalid URD file loaded:", "data/", path))
		print("Here's the original error message:")
		print(cond)
		return("<select>")
	})
}

plotDimPlot <- function(object) {
	if (seurat.version > "2" && seurat.version < "3") {
		DimPlot(
			object,
			reduction.use = "tsne",
			do.label = T,
			cols.use = colors.to.use
		)
	} else if (seurat.version > "3") {
		DimPlot(
			object,
			reduction = "tsne",
			label = T,
			cols = colors.to.use
		)
	}
}

plotViolinPlot <- function(object, genes) {
	if (seurat.version > "2" && seurat.version < "3") {
		VlnPlot(
			object,
			features = genes,
			cols.use = colors.to.use
		)
	} else if (seurat.version > "3") {
		VlnPlot(
			object,
			features = genes,
			cols = colors.to.use
		)
	}
}

plotFeaturePlot <- function(object, overlay, genes) {
	if (seurat.version > "2" && seurat.version < "3") {
		FeaturePlot(
			object,
			features = genes,
			reduction.use = "tsne",
			no.legend = F,
			cols.use = c("lightgray", "blue"),
			overlay = ((overlay == 1) && (length(genes) == 2))
		)
	} else if (seurat.version > "3") {
		FeaturePlot(
			object,
			features = genes,
			reduction = "tsne",
			blend = ((overlay == 1) && (length(genes) == 2))
		)
	}
}

meta.data <<- read.csv("data/ZFBrainAtlasMaster.csv")
meta.data <<- meta.data %>%
	select("STAGE", "CLUSTER", "ASSIGNED.CELL.TYPE.STATE", "ENRICHED.MARKERS")
meta.data <<-  meta.data %>%
	mutate(CLUSTER.NAME = paste(CLUSTER, ASSIGNED.CELL.TYPE.STATE, sep=" - "))

clustering <- c(
	"res.4.5",
	"res.4",
	"res.5",
	"res.5",
	"res.4.5",
	"res.5",
	"res.6",
	"res.6",
	"res.5",
	"res.5.5",
	"res.6",
	"res.5"
)

path <- c(
	"zf6s_cc_filt.cluster.rds",
	"zf10s_cc_filt.cluster.rds",
	"zf14s_cc_filt.cluster.rds",
	"zf18s_cc_filt.cluster.rds",
	"zf20hpf_cc_filt.cluster.rds",
	"zf24hpf_cc_filt.cluster.rds",
	"zf36hpf_cc_filt.cluster.rds",
	"zf2dpf_cc_filt.cluster.rds",
	"zf3dpf_cc_filt.cluster.rds",
	"zf5dpf_cc_filt.cluster.rds",
	"zf8dpf_cc_filt.cluster4.rds",
	"zf15dpf_PCAALL.rds"
)

stage <- c(
	"12 hpf",
	"14 hpf",
	"16 hpf",
	"18 hpf",
	"20 hpf",
	"24 hpf",
	"36 hpf",
	"2 dpf",
	"3 dpf",
	"5 dpf",
	"8 dpf",
	"15 dpf"
)

clustering.data <<- data.frame(stage, path, clustering)

colors.to.use <<- c(
	"pink1",
	"gold",
	"cyan3",
	"wheat3",
	"palegreen3",
	"skyblue1",
	"goldenrod",
	"sandybrown",
	"lightseagreen",
	"khaki3",
	"green3",
	"aquamarine3",
	"plum",
	"yellowgreen",
	"darkorange",
	"cornsilk3",
	"turquoise2",
	"orchid",
	"mediumpurple2",
	"chartreuse3",
	"seagreen3",
	"pink4",
	"slateblue1",
	"lemonchiffon2",
	"darkslategray3",
	"darkolivegreen3",
	"coral",
	"palevioletred",
	"darkgoldenrod2",
	"darksalmon",
	"sienna1",
	"steelblue1",
	"slategray2",
	"chocolate1",
	"lightgreen",
	"cadetblue3",
	"salmon",
	"yellow3",
	"thistle",
	"burlywood",
	"khaki",
	"mistyrose2",
	"azure3",
	"tomato",
	"darkseagreen3",
	"pink2",
	"lightgoldenrod",
	"peachpuff2",
	"grey70",
	"rosybrown3",
	"deepskyblue",
	"lavenderblush3",
	"lemonchiffon",
	"powderblue",
	"springgreen3",
	"bisque3",
	"lightcyan3",
	"gold3",
	"lightsalmon",
	"darkkhaki",
	"darkorchid2",
	"gray40",
	"hotpink",
	"lightgrey",
	"cyan2",
	"slategray3",
	"azure2",
	"wheat2",
	"mediumpurple",
	"darkgoldenrod3",
	"grey80",
	"mistyrose",
	"lavenderblush2",
	"seagreen2",
	"khaki2",
	"yellow2",
	"peachpuff",
	"rosybrown4",
	"springgreen2",
	"bisque2",
	"cornsilk2",
	"steelblue2",
	"sienna2",
	"aquamarine2",
	"palegreen2",
	"pink3",
	"darkolivegreen2",
	"steelblue3",
	"goldenrod3",
	"darkorchid",
	"brown",
	"indianred",
	"lemonchiffon4",
	"moccasin",
	"grey",
	"honeydew2",
	"coral2",
	"chocolate4",
	"darkturquoise",
	"firebrick2",
	"blueviolet"
)
