library(shinycssloaders)
library(Seurat)

stages <- c("<select>", list.files("data"))
clusters <- c("<select>")
genes <- c("<select>")

object <- NULL

clustering <- NULL

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

files <- c(
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
	"zf15dpf_PCAALL"
)

clustering.data <- data.frame(files, clustering)
