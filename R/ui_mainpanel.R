main_panel <- mainPanel(
  tabsetPanel(
    id = 'panel' ,
    
    tabPanel("About", column(8, includeMarkdown("README.md"))),
    
    
    tabPanel("Primers", dataTableOutput("table_primers")),
    
    tabPanel("Primer sets", dataTableOutput("table_primer_sets")),
    
    tabPanel("Amplification - overview",
             withSpinner(uiOutput("plot_matches_all"))),
    
    tabPanel("Amplification - details", 
             withSpinner(uiOutput("plot_matches_one"))),
    
    tabPanel("Test your primer/probe" ,
             withSpinner(uiOutput('primer_matches_user')),
             uiOutput('primer_matches_user_graph')),
    
    tabPanel("Test your primer set" , 
             withSpinner(uiOutput('primer_set_matches_user')),
             uiOutput('primer_set_matches_user_graph'))

    # tabPanel("Amplification - one set", dataTableOutput("test"))
  )
)