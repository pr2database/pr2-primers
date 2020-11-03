# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # 1 - Primer table
  output$table_primers <- renderDataTable({
    primers[, input$show_vars_primers, drop = FALSE]
  }, escape = FALSE,                                # To keep the HTML (doi links)
  options = list(pageLength = 25))
  
  output$table_primers_download <- downloadHandler(
          filename = function() {str_c("primers_", Sys.Date(), ".tsv")},
          content = function(path) {export(primers, file=path)}
  )      


  # 2 - Primer set table
  output$table_primer_sets <- renderDataTable({
        select(primer_sets[, input$show_vars_primer_sets, drop = FALSE], - primer_set_label_long)
      }, escape = FALSE,                                     # To keep the HTML (doi links)
      options = list(pageLength = 25))
  
  output$table_primer_sets_download <- downloadHandler(
          filename = function() {str_c("primer_sets_", Sys.Date(), ".tsv")},
          content = function(path) {export(primer_sets, file=path)}
  )      


  # 3 - Plot matches
  output$plot_matches_all <- renderUI({
        tagList(
          p(strong("Left panel:"), "% of sequence amplified.", 
             strong("Center panel:"), "number of mismatches.", 
             strong("Right panel:"), " Amplicon size"),
          renderPlot({plot_matches(input$kingdom_3, input$type)}, width = 1200, height = 1000, res = 96)
        )  
  })
  
  # 4 - Plot matches for one set of primer and one taxonomy group
    
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
          p(strong("Top panel:"), "Location of mismatches for fwd and reverse primer.", 
            strong("Bottom panel left"), "number of mismatches.",
            strong("Bottom panel right:"), "amplicon size"),
          renderPlot({
           data_4()$gg
          }, width = 1200, height=function(){300 + 90 + 20*data_4()$n_taxa}, res = 96)
        )
     }) 
     
  # 5 - Compute matches for user provided primers
     
  primer_match.df <- eventReactive(input$button_compute, ({
    fwd_valid <- primer_check(input$fwd_seq)
    shinyFeedback::feedback("fwd_seq", show = !fwd_valid, "Invalid fwd primer", color="red")
    
    rev_valid <- primer_check(input$rev_seq)
    shinyFeedback::feedback("rev_seq", show = !rev_valid, "Invalid rev primer", color="red")
    
    req(fwd_valid, rev_valid, cancelOutput = TRUE)

    primer_match(input$fwd_seq, input$rev_seq, as.integer(input$fwd_mis), as.integer(input$rev_mis) )
    # primer_match()
    })
  )
  
  primer_match_stats <- eventReactive(input$button_compute,({ primer_match.df() %>%
      group_by(kingdom) %>% 
      summarise(pct_fwd = sum(!is.na(fwd_pos))/n()*100,
                pct_rev = sum(!is.na(rev_pos))/n()*100,
                pct_amplified = sum(!is.na(ampli_size))/n()*100,
                mean_amplicon_size = mean(ampli_size, na.rm = TRUE)) %>%
      tidyr::pivot_longer(cols = pct_fwd:mean_amplicon_size, names_to = "Parameter", values_to = "Values") %>%
      mutate(Parameter = str_replace_all(Parameter,c("_" = " ",
                                                     "pct" = "% sequences",
                                                     "fwd" = "matching forward primer",
                                                     "rev" = "matching reverse primer")))
  }) 
  )


  output$primer_matches_user <- renderUI({
    req(primer_match.df())
    tagList(
      p("Overall statistics"),
      p(),
      renderTable(primer_match_stats(), width = 600, colnames = FALSE),
      downloadHandler(
          filename = function() {str_c("primer_match_pr2_", Sys.Date(), ".tsv")},
          content = function(path) {export(primer_match.df(), file=path)},
          outputArgs = list(label = "Download results"),
      )
    )  
  })  
  
    plot_matches_simple_taxa_data <- eventReactive(
      {input$button_plot_user
       input$button_compute },
      {
      # --- Determine at which level do the plot depending on what has been selected
      taxo <- taxo_selected(input$kingdom, input$supergroup, input$division, input$class)   
      plot_matches_simple_taxa(primer_match.df(), taxo$level, taxo$name)
    })
  

  output$primer_matches_user_graph <- renderUI({
    req(primer_match.df())
    tagList(
       p(""),
       renderPlot({
         plot_matches_simple_taxa_data()$gg
        }, width = 1200, height=function(){90 + 30*plot_matches_simple_taxa_data()$n_taxa}, res = 96)
    )  
  })  
        
   
     

  # Utils - Dynamic taxonomy boxes
  
    # --- Supergroups
    supergroups <- reactive({
      req(input$kingdom)
      filter(pr2_taxo, kingdom == input$kingdom) %>% 
        pull(supergroup) %>%  
        unique()  
    })
    
    observeEvent(supergroups(), {
      updateSelectInput(session, "supergroup", choices = c("All", supergroups()))
    })  
  
    # --- Divisions
    divisions <- reactive({
      req(input$supergroup)
      filter(pr2_taxo, supergroup == input$supergroup) %>% 
        pull(division) %>%  
        unique()  
    })
    
    observeEvent(divisions(), {
      updateSelectInput(session, "division", choices = c("All", divisions()))
    })

    # --- Classes    
    classes <- reactive({
      req(input$division)
      filter(pr2_taxo, division == input$division) %>% 
        pull(class) %>%  
        unique()
    })
    
    observeEvent(classes(), {
      updateSelectInput(session, "class", choices = c("All", classes()))
    })
    
   
  # Test
  output$test <- renderDataTable(plot_matches_sg(input$primer_set_id))    
}

# Run the shiny app -------------------------------------------------------


# shinyApp(ui, server)