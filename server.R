
server <- function(input, output, session) {
	get_clusters = reactive({
		# object <<- UpdateSeuratObject(readRDS(paste0("data/", input$stage)))
		object <<- readRDS(paste0("data/", input$stage))
		# c("<select>", "all", levels(unique(object@active.ident)))
		c("<select>", "all", levels(unique(object@ident)))
	})
	
	get_genes = reactive({
		# c("<select>", rownames(object))
		c("<select>", rownames(object@raw.data))
	})
	observe({
		input$stage
		print(paste0("Selected ", input$stage))
		validate(need(input$stage != "<select>", "PLease select stage to be laoded."))
		# object <<- UpdateSeuratObject(readRDS(paste0("data/", input$stage)))
		object <<- readRDS(paste0("data/", input$stage))
		updateSelectInput(
			session = session, 
			inputId = "cluster",
			choices = get_clusters())
		updateSelectInput(
			session = session,
			inputId = "gene",
			choices = get_genes())
	})
	
	output$umap <- renderPlot({
		validate(need(input$stage != "<select>", "PLease select stage to be laoded."))
		validate(need(input$cluster != "<select>", "PLease select cluster(s) to be plotted."))
		if (input$cluster == "all") {
			# DimPlot(object)
			DimPlot(object, reduction.use = "tsne")
		} else {
			# DimPlot(object[,object@active.ident == input$cluster])
			DimPlot(SubsetData(object, cells.use = object@ident == input$cluster), reduction.use = "tsne")
		}

	})
	
	output$violin <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		VlnPlot(object, features = input$gene)
	})
	
	output$umap_gene <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		# FeaturePlot(object, feature = input$gene)
		FeaturePlot(object, feature = input$gene, reduction.use = "tsne", no.legend = F)
	})
}