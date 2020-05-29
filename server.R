
server <- function(input, output, session) {
	get_clusters = reactive({
		# object <<- UpdateSeuratObject(readRDS(paste0("data/", input$path)))
		object <<- readRDS(paste0("data/", input$path))
		# c("<select>", "all", levels(unique(object@active.ident)))
		c("<select>", "all", levels(unique(object@ident)))
	})

	get_genes = reactive({
		# c("<select>", rownames(object))
		c("<select>", sort(rownames(object@raw.data)))
	})
	observe({
		input$path
		print(paste0("Selected ", input$path))
		validate(need(input$path != "<select>", "PLease select file to be laoded."))
		# object <<- UpdateSeuratObject(readRDS(paste0("data/", input$path)))
		object <<- readRDS(paste0("data/", input$path))
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
		validate(need(input$path != "<select>", "PLease select file to be laoded."))
		validate(need(input$cluster != "<select>", "PLease select cluster(s) to be plotted."))
		if ("all" %in% input$cluster) {
			# DimPlot(object)
			DimPlot(
				object, 
				reduction.use = "tsne", 
				do.label = T,
				cols.use = colors.to.use
			)
		} else {
			# DimPlot(object[,object@active.ident == input$cluster])
			DimPlot(
				SubsetData(object, cells.use = (object@ident %in% input$cluster)), 
				reduction.use = "tsne", 
				do.label = T,
				cols.use = colors.to.use
			)
		}

	})
	
	output$violin <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		if ("all" %in% input$cluster) {
			VlnPlot(
				object,
				features = input$gene,
				cols.use = colors.to.use
			)
		} else {
			VlnPlot(
				SubsetData(object, cells.use = (object@ident %in% input$cluster)),
				features = input$gene,
				cols.use = colors.to.use
			)
		}
		
	})
	
	output$umap_gene <- renderPlot({
		validate(need(input$gene != "<select>", "PLease select gene to be plotted over given clusters."))
		validate(need(length(input$gene) < 3, "PLease select 1 or 2 genes."))
		# FeaturePlot(object, feature = input$gene)
		if ("all" %in% input$cluster) {
			FeaturePlot(
				object,
				features = input$gene,
				reduction.use = "tsne",
				no.legend = F,
				cols.use = c("lightgray", "blue"),
				overlay = ((input$overlay == 1) && (length(input$gene) == 2))
			)
		} else {
			FeaturePlot(
				SubsetData(object, cells.use = (object@ident %in% input$cluster)),
				features = input$gene,
				reduction.use = "tsne",
				no.legend = F,
				cols.use = c("lightgray", "blue"),
				overlay = ((input$overlay == 1) && (length(input$gene) == 2))
			)
		}
 	})
	
	output$metaTable <- DT::renderDataTable({
		validate(need(input$cluster != "<select>", "PLease select cluster(s)."))
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