library(markdown)

navbarPage("[INSERT TITLE]!",
		   tabPanel("[INSERT SUBTITLE]",
		   		 sidebarLayout(
		   		 	sidebarPanel(
		   		 		paste0("Seurat version used: ", seurat.version),
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

