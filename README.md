# github.com/brlauuu/seurat_exploration

Simple R Shiny app for exploring different [Seurat](https://satijalab.org/seurat/) objects.

In order to use the app, add your Seurat objects into `data` directory and start the Shiny app.

Your `Seurat` files must be built under version `v2.3.4`. 

Requirements:

* Seurat version v2.3.4 install using instructions from [https://satijalab.org/seurat/install.html](https://satijalab.org/seurat/install.html)
* Some of the requirements that you need to install separately:
    * `mclust`
    * `Hmisc`
    * `latticeExtra`
    * devtools::install_github("cran/SDMTools")
    * devtools::install_github("cran/multtest")
