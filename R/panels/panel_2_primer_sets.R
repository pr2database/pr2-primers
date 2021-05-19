output$table_primer_sets <- renderDataTable({
  select(primer_sets[, input$show_vars_primer_sets, drop = FALSE], - primer_set_label_long)
}, escape = FALSE,                                     # To keep the HTML (doi links)
options = list(pageLength = 25))

output$table_primer_sets_download <- downloadHandler(
  filename = function() {str_c("primer_sets_", Sys.Date(), ".tsv")},
  content = function(path) {export(select(primer_sets, -doi_html, -metabarcoding_doi_html), file=path)}
) 