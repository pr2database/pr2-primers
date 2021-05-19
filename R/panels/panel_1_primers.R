output$table_primers <- renderDataTable({
  primers[, input$show_vars_primers, drop = FALSE]
}, escape = FALSE,                                # To keep the HTML (doi links)
options = list(pageLength = 25))

output$table_primers_download <- downloadHandler(
  filename = function() {str_c("primers_", Sys.Date(), ".tsv")},
  content = function(path) {export(select(primers, -doi_html), file=path)}
)      
