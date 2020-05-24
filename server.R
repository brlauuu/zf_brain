
server <- function(input, output, session) {
	get_clusters = reactive({
		object <<- readRDS(paste0("data/", input$stage))
		c("all", levels(unique(object$seurat_clusters)))
	})
	
	get_genes = reactive({
		c("<select>", rownames(object))
	})
	observe({
		input$stage
		print(paste0("Selected ", input$stage))
		validate(need(input$stage != "<select>", "PLease select stage to be laoded."))
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
		if (input$cluster == "all") {
			DimPlot(object)
		} else {
			DimPlot(object[,object$seurat_clusters == input$cluster])
		}

	})
	
	output$violin <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		VlnPlot(object, features = input$gene)
	})
	
	output$umap_gene <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		FeaturePlot(object, feature = input$gene)
	})
}