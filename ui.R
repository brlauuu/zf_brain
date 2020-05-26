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
		   		 			# column(4,
		   		 				   withSpinner(plotOutput("umap")),
		   		 			# ),
		   		 			# column(4,
		   		 				   withSpinner(plotOutput("violin")),
		   		 			# ),
		   		 			# column(4,
		   		 				   withSpinner(plotOutput("umap_gene"))
		   		 			# )
		   		 		)
		   		 		
		   		 	)
		   		 )
		   ),
   		   tabPanel("About",
   		   		 verbatimTextOutput("About")
   		   )
)

