
server <- function(input, output, session) {
	get_clusters = reactive({
		return(listClusters(input$path))
	})

	get_genes = reactive({
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
			choices = get_clusters(),
			selected = "all")
		updateSelectInput(
			session = session,
			inputId = "gene",
			choices = get_genes())
	})

	output$umap <- renderPlot({
		validate(
			need(
				input$path != "<select>" & input$cluster != "<select>",
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
		filename = function() { paste("tSNE_snapshot", '.png', sep='') },
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
		filename = function() { paste("cluster_expression_snapshot", '.png', sep='')},
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
		filename = function() { paste("violoin_plot_snapshot", '.png', sep='')},
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
}