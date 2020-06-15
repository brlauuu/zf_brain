
server <- function(input, output, session) {

	##############
	# Seurat
	##############
	getClusters = reactive({
		return(listClusters(input$path))
	})

	getGenes = reactive({
		return(listGenes())
	})

	observe({
		print(paste0("Selected ", input$path))
		if (input$path == "<select>") {
			return()
		}
		loadObject(input$path)
		updateSelectInput(
			session = session,
			inputId = "cluster",
			choices = getClusters(),
			selected = "all")
		updateSelectInput(
			session = session,
			inputId = "gene",
			choices = getGenes())
	})

	output$umap <- renderPlot({
		validate(
			need(
				input$path != "<select>" && input$cluster != "<select>",
				"Please select file to be laoded and cluster to be plotted."
			)
		)
		if ("all" %in% input$cluster) {
			plotDimPlot(object)
		} else {
			selected.clusters <- meta.data %>%
				filter(CLUSTER.NAME %in% input$cluster) %>%
				select(CLUSTER) %>%
				unlist()
			plotDimPlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% selected.clusters)
				)
			)
		}
	})

	output$violin <- renderPlot({
		validate(
			need(
				input$gene != "<select>",
				"Please select gene to be plotted over given clusters."
			)
		)
		if ("all" %in% input$cluster) {
			plotViolinPlot(object, input$gene)
		} else {
			selected.clusters <- meta.data %>%
				filter(CLUSTER.NAME %in% input$cluster) %>%
				select(CLUSTER) %>%
				unlist()
			plotViolinPlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% selected.clusters)
				),
				input$gene
			)
		}
	})

	output$umap_gene <- renderPlot({
		validate(
			need(
				input$gene != "<select>" && length(input$gene) < 3,
				"Please select 1 or 2 gene(s) to be plotted over given clusters."
			)
		)
		if ("all" %in% input$cluster) {
			plotFeaturePlot(object, input$overlay, input$gene)
		} else {
			selected.clusters <- meta.data %>%
				filter(CLUSTER.NAME %in% input$cluster) %>%
				select(CLUSTER) %>%
				unlist()
			plotFeaturePlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% selected.clusters)
				),
				input$overlay,
				input$gene
			)
		}
 	})

	output$metaTable <- DT::renderDataTable({
		validate(
			need(
				input$cluster != "<select>",
				""
			)
		)
		if ("all" %in% input$cluster) {
			subset(
				meta.data,
				subset=(
					STAGE==clustering.data$stage[clustering.data$path==input$path]
				)
			) %>% select(STAGE, CLUSTER.NAME, ENRICHED.MARKERS)
		} else {
			selected.clusters <- meta.data %>%
				filter(CLUSTER.NAME %in% input$cluster) %>%
				select(CLUSTER) %>%
				unlist()
			subset(
				meta.data,
				subset=(
					STAGE==clustering.data$stage[clustering.data$path==input$path]
					& CLUSTER %in% selected.clusters
				)
			)  %>% select(STAGE, CLUSTER.NAME, ENRICHED.MARKERS)
		}
	})

	output$downloadMainTsne <- downloadHandler(
		filename = function() { paste0("tSNE_snapshot.", input$format) },
		content = function(file) {
			if ("all" %in% input$cluster) {
				ggsave(
					file,
					plot = plotDimPlot(object),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			} else {
				selected.clusters <- meta.data %>%
					filter(CLUSTER.NAME %in% input$cluster) %>%
					select(CLUSTER) %>%
					unlist()
				ggsave(
					file,
					plot = plotDimPlot(
						SubsetData(
							object,
							cells.use = (object@ident %in% selected.clusters)
						)
					),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			}
		}
	)

	output$downloadGeneOverlay <- downloadHandler(
		filename = function() { paste0("cluster_expression_snapshot.", input$format) },
		content = function(file) {
			if ("all" %in% input$cluster) {
				ggsave(
					file,
					plot = plotFeaturePlot(object, input$overlay, input$gene),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			} else {
				selected.clusters <- meta.data %>%
					filter(CLUSTER.NAME %in% input$cluster) %>%
					select(CLUSTER) %>%
					unlist()
				ggsave(
					file,
					plot = plotFeaturePlot(
						SubsetData(
							object,
							cells.use = (object@ident %in% selected.clusters)
						),
						input$overlay,
						input$gene
					),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			}
		}
	)

	output$downloadViolin <- downloadHandler(
		filename = function() { paste0("violoin_plot_snapshot.", input$format) },
		content = function(file) {
			if ("all" %in% input$cluster) {
				ggsave(
					file,
					plot = plotViolinPlot(object, input$gene),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			} else {
				selected.clusters <- meta.data %>%
					filter(CLUSTER.NAME %in% input$cluster) %>%
					select(CLUSTER) %>%
					unlist()
				ggsave(
					file,
					plot = plotViolinPlot(
						SubsetData(
							object,
							cells.use = (object@ident %in% selected.clusters)
						),
						input$gene
					),
					width = as.double(input$width),
					height = as.double(input$height),
					units = input$units,
					device = input$format
				)
			}
		}
	)

	##############
	# URD
	##############
	getFeaturesURD = reactive({
		return(listFeaturesURD())
	})

	observe({
		print(paste0("Selected ", input$path.urd))
		if (input$path.urd == "<select>") {
			return()
		}
		loadObjectURD(input$path.urd)
		updateSelectInput(
			session = session,
			inputId = "feature.urd",
			choices = getFeaturesURD())
	})

	output$tree <- renderPlot({
		validate(
			need(
				input$path.urd != "<select>",
				"Please select file to be laoded."
			)
		)
		validate(
			need(
				length(input$feature.urd) < 3 && length(input$feature.urd) > 0,
				"Select at least 1 or 2 features to be plotted on the tree."
			)
		)
		if (length(input$feature.urd) < 2) {
			plotTree(
				object.urd,
				input$feature.urd
			)
		} else {
			plotTreeDual(
				object.urd,
				input$feature.urd[[1]],
				input$feature.urd[[2]]
			)

		}
	})

	output$downloadTree <- downloadHandler(
		filename = function() { paste0("tree_plot_snapshot.", input$format.urd) },
		content = function(file) {
			if (length(input$feature.urd) < 2) {
				ggsave(
					file,
					plot = plotTree(
						object.urd,
						input$feature.urd
					),
					width = as.double(input$width.urd),
					height = as.double(input$height.urd),
					units = input$units.urd,
					device = input$format.urd
				)
			} else {
				ggsave(
					file,
					plot = plotTreeDual(
						object.urd,
						input$feature.urd[[1]],
						input$feature.urd[[2]]
					),
					width = as.double(input$width.urd),
					height = as.double(input$height.urd),
					units = input$units.urd,
					device = input$format.urd
				)
			}
		}
	)
}