library(shinycssloaders)
library(Seurat)

stages <- c("<select>", list.files("data"))
clusters <- c("<select>")
genes <- c("<select>")

object <- NULL