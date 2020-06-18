library(markdown)

shinyUI(
	navbarPage("Emergence of neuronal diversity during vertebrate brain development",
		   tabPanel("Seurat",
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
						selectInput(
							"format",
							"Downlaod format",
							choices=device),
						selectInput(
							"units",
							"Size units",
							choices=units),
						textInput("height", "height", "200"),
						textInput("width", "width", "400"),
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
		   tabPanel("URD",
		   		 sidebarLayout(
		   		 	sidebarPanel(
		   		 		helpText(paste0("URD v", urd.version)),
		   		 		selectInput(
		   		 			"path.urd",
		   		 			"Select file",
		   		 			choices=path.to.load),
		   		 		selectInput(
		   		 			"feature.urd",
		   		 			"Select feature(s)",
		   		 			choices=features.urd,
		   		 			multiple = T),
		   		 		selectInput(
		   		 			"format.urd",
		   		 			"Downlaod format",
		   		 			choices=device),
		   		 		selectInput(
		   		 			"units.urd",
		   		 			"Size units",
		   		 			choices=units),
		   		 		textInput("height.urd", "height", "200"),
		   		 		textInput("width.urd", "width", "400"),
		   		 		downloadButton(
		   		 			"downloadTree",
		   		 			label = "Tree"
		   		 		),
		   		 		width = "3"
		   		 	),
		   		 	mainPanel(
		   		 		fluidRow(
		   		 			withSpinner(plotOutput("tree")),
		   		 		)
		   		 	)
		   		 )
		   ),
   		   tabPanel("About",
   		   		 includeMarkdown("README.md")
   		   )
	)
)

