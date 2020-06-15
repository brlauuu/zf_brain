# Emergence of neuronal diversity during vertebrate brain development

Raj, B., Farrell, J.A., McKenna, A., Leslie, J.L. and Schier, A.F., 2019. Emergence of neuronal diversity during vertebrate brain development. bioRxiv, p.839860.

https://doi.org/10.1101/839860

Co-corresponding authors. Email: bushranraj@gmail.com(B.R.); schier@mcb.harvard.edu(A.F.S.)

Neurogenesis  in  the  vertebrate  brain  comprises  many  steps  ranging  from  the  proliferation  of progenitors to the differentiation and maturation of neurons. Although these processes are highly regulated,  the  landscape  of  transcriptional  changes  and  progenitoridentities  underlying  brain development are poorly characterized. Here, we describe the first developmental single-cell RNA-seq catalog of more than 200,000 zebrafish brain cells encompassing 12 stages from 12 hours post-fertilization to 15 days post-fertilization. We characterize known and novel gene markers for more than 800 clusters across these timepoints. Our results capture the temporal dynamics of multiple neurogenic waves from embryo to larva that expand neuronal diversity from ~20cell types at 12hpf to ~100 cell types at 15 dpf. We find that most embryonic neural progenitor states are transient and transcriptionally distinct from long-lasting neural progenitors of post-embryonic stages.  Furthermore,  we  reconstruct  cell  specification  trajectories for  the  retina  and hypothalamus, and identify gene expression cascades and novel markers. Our analysis reveal that  late-stage  retinal  neural  progenitors  transcriptionally  overlap  cell  states  observed  in  the embryo, while hypothalamic neural progenitors become progressively distinct with developmental time.  These  data  provide  the  first  comprehensive  single-cell  transcriptomic  time  course  for vertebrate  brain  development  and  suggest  distinct  neurogenic  regulatory  paradigms  between different stages and tissues. 

## Exploring the dataset using R Shiny App

Thank you for using this ShinyApp for exploring data from the work cited above.
In order to use the app, you have to follow the following steps:

* Download prepared `Seurat` and `URD` fils from [LINK] to `data` directory
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
for each of th plots.

### `Seurat`

* Select all, one or several clusters to plot on the tSNE coordinates of.
* Select one, two (combined tSNE plot) or several (sequential violin plot)
genes to plot.

### `URD`

* Select one or two features to plot on the tree. Note: discrete and float
features cannot be combined.

## Error troubleshooting

In case of unexpected behavior, please consult first the log in the R console.
If that did not help, do not hesitate to contact me at dorde.relic@unibas.ch.
