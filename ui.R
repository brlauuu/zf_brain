library(markdown)

navbarPage("Emergence of neuronal diversity during vertebrate brain development",
		   tabPanel("[INSERT SUBTITLE]",
		   		 sidebarLayout(
		   		 	sidebarPanel(
		   		 		helpText(paste0("Seurat v", seurat.version)),
		   		 		selectInput(
		   		 			"path",
		   		 			"Select file",
		   		 			choices=path.to.load),
		   		 		selectInput(
		   		 			"cluster",
		   		 			"Select cluster(s) -- n >= 1",
		   		 			choices=clusters,
		   		 			multiple = T),
		   		 		selectInput(
		   		 			"gene",
		   		 			"Select gene(s) -- 0 <= n <= 2",
		   		 			choices=genes,
		   		 			multiple = T),
		   		 		radioButtons("overlay",
							label = "2 gene plot:",
							choices = list("Overlay" = 1,"Separate" = 2)
		   		 		),
						downloadButton(
							"downloadMainTsne", 
							label = "Clusters tSNE"
						),
						downloadButton(
						   	"downloadGeneOverlay", 
						   	label = "Gene Expr tSNE"
						),
						downloadButton(
					   		"downloadViolin", 
					   		label = "Violin Plot"
						),
						width = "3"
		   		 	),
		   		 	mainPanel(
		   		 		fluidRow(
		   		 			withSpinner(DT::dataTableOutput("metaTable")),
							withSpinner(plotOutput("umap")),
							withSpinner(plotOutput("umap_gene")),
							withSpinner(plotOutput("violin"))
		   		 		)
		   		 	)
		   		 )
		   ),
   		   tabPanel("About",
   		   		 includeMarkdown("README.md")
   		   )
)

