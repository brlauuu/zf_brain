
# Emergence of neuronal diversity during vertebrate brain development

Bushra Raj, Ph.D.; Jeffrey A. Farrell; Jialin Liu; Jakob El Kholtei; Adam N. Carte; Joaquin Navajas Acedo; Lucia Y. Du; Aaron McKenna; Đorđe Relić; Jessica M. Leslie; Alexander F. Schier

Neurogenesis comprises many highly regulated processes including proliferation, differentiation, and maturation. However, the transcriptional landscapes underlying brain development are poorly characterized. We describe a developmental single-cell catalog of ∼220,000 zebrafish brain cells encompassing 12 stages from embryo to larva. We characterize known and novel gene markers for ∼800 clusters and provide an overview of the diversification of neurons and progenitors across these time points. We also introduce an optimized GESTALT lineage recorder that enables higher expression and recovery of Cas9-edited barcodes to query lineage segregation. Cell type characterization indicates that most embryonic neural progenitor states are transitory and transcriptionally distinct from neural progenitors of post-embryonic stages. Reconstruction of cell specification trajectories reveals that late-stage retinal neural progenitors transcriptionally overlap cell states observed in the embryo. The zebrafish brain development atlas provides a resource to define and manipulate specific subsets of neurons and to uncover the molecular mechanisms underlying vertebrate neurogenesis.

Co-corresponding authors. Bushra Raj (bushranraj@gmail.com); Alexander F. Schier (alex.schier@unibas.ch)

https://doi.org/10.1016/j.neuron.2020.09.023

## Exploring the dataset using R Shiny App

Thank you for using this ShinyApp for exploring data from the work cited above.
In order to use the app, you have to follow the following steps:

* Download (and decompress!) `Seurat` and `URD` files from [GEO: GSE158142](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE158142) to `data` directory
* Install [`Seurat`](https://github.com/satijalab/seurat),
[`URD`](https://github.com/farrellja/URD) and
[`Shiny`](https://shiny.rstudio.com/) libraries.
* Start the ShinyApp by opening `ui.R` from RStudio and clicking `RunApp`.

Seurat objects were built using `Seurat v2.3.4`. We recommend using this
version of Seurat, however, we did implement conversion of objects to the new
Seurat version as well. In case you want to use `v2.3.4` of Seurat, see
instructions on how to set it up here
[https://satijalab.org/seurat/install.html](https://satijalab.org/seurat/install.html).

### Using ShinyApp

When you start the app you can see 3 tabs. One for exploring the `Seurat`
generated objects, one for exploring `URD` generated objects and this page.

Once you downloaded the right files to the `data` directory, the files are
going to appear in the "Select file" dropdown menu. Depending on the library
you are using, you need to select the right files.

Once you loaded a file, plots will be automatically generated. After that, you
can select metadata of interest that you would like to explore.

For both libraries, on the bottom of the sidebar, you can find download buttons
for each of the plots.

### `Seurat`

* Select all, one or several clusters to plot on the tSNE coordinates of.
* Select one, two (combined tSNE plot) or several (sequential violin plot)
genes to plot.

### `URD`

* Select one or two features to plot on the tree. Note: discrete and float
features cannot be combined.

## Error troubleshooting

In case of unexpected behavior, please consult first the log in the R console.
If that did not help, open an [new issue](https://github.com/brlauuu/zf_brain/issues) or contact me at dorde.relic@unibas.ch.
