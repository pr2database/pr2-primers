# --- Create the graph and count number of taxa    

data_4 <- eventReactive({input$button_plot
  input$primer_set_id},{
    # --- Determine at which level do the plot depending on what has been selected
    taxo <- taxo_selected(input$kingdom, input$supergroup, input$division, input$class) 
    plot_matches_detailed_taxa(input$primer_set_id, taxo$level, taxo$name)
  })

output$plot_matches_one <- renderUI({
  # req(plot_matches_detailed_taxa_data())
  
  tagList( 
    p("Overall statistics"),
    p(),
    renderTable(data_4()$summary_kingdom, width = 600, colnames = FALSE),
    p(strong("Top panel:"), "Location of mismatches for forward and reverse primer.", 
      strong("Bottom panel left"), "number of mismatches.",
      strong("Bottom panel right:"), "amplicon size"),
    renderPlot({
      data_4()$gg
    }, width = 1200, height=function(){300 + 90 + 20*data_4()$n_taxa}, res = 96)
  )
}) 