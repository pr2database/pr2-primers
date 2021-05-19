print("Load 6")

primer_match.df <- eventReactive(input$button_compute_primer, ({
  primer_valid <- primer_check(input$primer_seq)
  shinyFeedback::feedback("primer_seq", show = !primer_valid, "Invalid primer", color="red")
  
  req(primer_valid, cancelOutput = TRUE)

  primer_match(input$primer_seq, as.integer(input$primer_mis), input$primer_type)
})
)

primer_match_stats <- eventReactive(input$button_compute_primer,({ primer_match.df() %>%
    group_by(kingdom) %>%
    summarise(pct_fwd = sum(!is.na(fwd_pos))/n()*100) %>%
    tidyr::pivot_longer(cols = pct_fwd, names_to = "Parameter", values_to = "Values") %>%
    mutate(Parameter = str_replace_all(Parameter,c("_" = " ",
                                                   "pct" = "% sequences",
                                                   "fwd" = "matching primer/probe")))
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

plot_primer_matches_simple_taxa_data <- eventReactive(
  {input$button_plot_primer_user
    input$button_compute_primer },
  {
    # --- Determine at which level do the plot depending on what has been selected
    taxo <- taxo_selected(input$kingdom, input$supergroup, input$division, input$class)
    plot_primer_matches_simple_taxa(primer_match.df(), taxo$level, taxo$name)
  })


output$primer_matches_user_graph <- renderUI({
  req(primer_match.df())
  tagList(
    p(""),
    renderPlot({
      plot_primer_matches_simple_taxa_data()$gg
    }, width = 1200, height=function(){90 + 30*plot_primer_matches_simple_taxa_data()$n_taxa}, res = 96)
  )
})
