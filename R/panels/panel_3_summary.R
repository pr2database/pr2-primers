output$plot_matches_all <- renderUI({
  tagList(
    p(strong("Left panel:"), "% of sequence amplified.", 
      strong("Center panel:"), "number of mismatches.", 
      strong("Right panel:"), " Amplicon size"),
    renderPlot({plot_matches(input$kingdom_3, input$type)}, width = 1200, height = 1000, res = 96)
  )  
})