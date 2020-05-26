library(markdown)

navbarPage("[INSERT TITLE]!",
		   tabPanel("[INSERT SUBTITLE]",
		   		 sidebarLayout(
		   		 	sidebarPanel(
		   		 		selectInput(
		   		 			"stage",
		   		 			"Select stage",
		   		 			choices=stages),
		   		 		selectInput(
		   		 			"cluster",
		   		 			"Select clusters",
		   		 			choices=clusters),
		   		 		selectInput(
		   		 			"gene",
		   		 			"Select genes",
		   		 			choices=genes),
		   		 		width = "3"
		   		 	),
		   		 	mainPanel(
		   		 		fluidRow(
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

