
server <- function(input, output, session) {
	get_clusters = reactive({
		return(listClusters())
	})

	get_genes = reactive({
		return()
	})
	observe({
		input$path
		print(paste0("Selected ", input$path))
		validate(
			need(
				input$path != "<select>",
				"Please select file to be laoded."
			)
		)
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
				input$path != "<select>",
				"Please select file to be laoded."
			)
		)
		validate(
			need(
				input$cluster != "<select>",
				"Please select cluster(s) to be plotted."
			)
		)
		if ("all" %in% input$cluster) {
			plotDimPlot(object)
		} else {
			plotDimPlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% input$cluster)
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
			plotViolinPlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% input$cluster)
				),
				input$gene
			)
		}
	})

	output$umap_gene <- renderPlot({
		validate(
			need(
				input$gene != "<select>",
				"Please select gene to be plotted over given clusters."
			)
		)
		validate(
			need(
				length(input$gene) < 3,
				"Please select 1 or 2 genes."
			)
		)
		# FeaturePlot(object, feature = input$gene)
		if ("all" %in% input$cluster) {
			plotFeaturePlot(object, input$overlay, input$gene)
		} else {
			plotFeaturePlot(
				SubsetData(
					object,
					cells.use = (object@ident %in% input$cluster)
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
				"Please select cluster(s)."
			)
		)
		if ("all" %in% input$cluster) {
			subset(
				meta.data,
				subset=(
					STAGE==clustering.data$stage[clustering.data$path==input$path]
				)
			)
		} else {
			subset(
				meta.data,
				subset=(
					STAGE==clustering.data$stage[clustering.data$path==input$path]
					& CLUSTER %in% input$cluster
				)
			)
		}
	})
}