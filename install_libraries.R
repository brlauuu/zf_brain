if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    shiny,
    shinycssloaders,
    devtools,
    Seurat,
    URD,
    DT,
    dplyr,
    ggplot2,
    install = TRUE
)