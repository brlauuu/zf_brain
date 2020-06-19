if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    shiny,
    shinycssloaders,
    devtools,
    Seurat,
    DT,
    dplyr,
    ggplot2,
    install = TRUE
)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("destiny")) BiocManager::install("destiny")

if (!require("URD")) install_github(repo = "farrellja/URD")