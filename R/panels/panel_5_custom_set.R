primer_set_match.df <- eventReactive(input$button_compute_primer_set, ({
  fwd_valid <- primer_check(input$fwd_seq)
  shinyFeedback::feedback("fwd_seq", show = !fwd_valid, "Invalid fwd primer", color="red")
  
  rev_valid <- primer_check(input$rev_seq)
  shinyFeedback::feedback("rev_seq", show = !rev_valid, "Invalid rev primer", color="red")
  
  req(fwd_valid, rev_valid, cancelOutput = TRUE)
  
  primer_set_match(input$fwd_seq, input$rev_seq, as.integer(input$fwd_mis), as.integer(input$rev_mis) )
  # primer_set_match()
})
)

primer_set_match_stats <- eventReactive(input$button_compute_primer_set,({ primer_set_match.df() %>%
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


output$primer_set_matches_user <- renderUI({
  req(primer_set_match.df())
  tagList(
    p("Overall statistics"),
    p(),
    renderTable(primer_set_match_stats(), width = 600, colnames = FALSE),
    downloadHandler(
      filename = function() {str_c("primer_set_match_pr2_", Sys.Date(), ".tsv")},
      content = function(path) {export(primer_set_match.df(), file=path)},
      outputArgs = list(label = "Download results"),
    )
  )  
})  

plot_primer_set_matches_simple_taxa_data <- eventReactive(
  {input$button_plot_primer_set_user
    input$button_compute_primer_set },
  {
    # --- Determine at which level do the plot depending on what has been selected
    taxo <- taxo_selected(input$kingdom, input$supergroup, input$division, input$class)   
    plot_primer_set_matches_simple_taxa(primer_set_match.df(), taxo$level, taxo$name)
  })


output$primer_set_matches_user_graph <- renderUI({
  req(primer_set_match.df())
  tagList(
    p(""),
    renderPlot({
      plot_primer_set_matches_simple_taxa_data()$gg
    }, width = 1200, height=function(){90 + 30*plot_primer_set_matches_simple_taxa_data()$n_taxa}, res = 96)
  )  
})  
