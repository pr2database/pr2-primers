# User interface ----------------------------------------------------------

ui <- fluidPage(
  
  # Script to close the windows after some inactivity
  
  tags$script(inactivity), 
  
  # Title
  title = "Primer database",
  titlePanel(div(img(src='pr2_logo.png', width="80"),"The PR2 primer database")),
  
  
  # --- Side bar
  sidebarLayout(side_bar, main_panel)
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Stop the application of the session is closed (after 10 min)
  session$onSessionEnded(stopApp)

  # 1 - Primer table
  
    source("R/panels/panel_1_primers.R", local = TRUE)

  # 2 - Primer set table
     
    source("R/panels/panel_2_primer_sets.R", local = TRUE)
  

  # 3 - Plot matches
  
   source("R/panels/panel_3_summary.R", local = TRUE) 
  
  # 4 - Plot matches for one set of primer and one taxonomy group
    
    source("R/panels/panel_4_one_set.R", local = TRUE)
  
  # 5 - Compute matches for user provided primer set
  
    source("R/panels/panel_5_custom_set.R", local = TRUE)
     
  # 6 - Compute matches for user provided primer
  
    source("R/panels/panel_6_custom_primer.R", local = TRUE)
  
  # Utils - Dynamic taxonomy boxes
  
    source("R/panels/dynamic_taxo_boxes.R", local = TRUE)
  
    
   
  # Test
  output$test <- renderDataTable(plot_matches_sg(input$primer_set_id))    
}

# Run the shiny app -------------------------------------------------------


shinyApp(ui, server)
