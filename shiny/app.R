source("read_data.R")
source("graphics.R")

# User interface ----------------------------------------------------------

ui <- fluidPage(
  title = "Primer databases",
  titlePanel("Primer database for protist rRNA genes"),
  fluidRow(
    column(5, includeMarkdown("readme.md"))
           ),
  sidebarLayout(
    sidebarPanel(width = 3,
      conditionalPanel(
        'input.dataset === "Primers"',
        includeMarkdown("primers.md"),
        checkboxGroupInput("show_vars_primers", "Columns to show:",
                           names(primers), selected = names(primers))
      ),
      conditionalPanel(
        'input.dataset === "Primer sets"',
        includeMarkdown("primer_sets.md"),
        checkboxGroupInput("show_vars_primer_sets", "Columns to show:",
                           names(primer_sets), selected = names(primer_sets))
      ),conditionalPanel(
        'input.dataset === "Amplification - eukaryotes"',
          radioButtons("specificity", "Primer type:",
               c("General" = "general",
                 "Specific" = "specific")),
      ),conditionalPanel(
        'input.dataset === "Amplification - by group"',
          selectInput("primer_set_id", "Primer set",
               setNames(primer_sets_tested$primer_set_id, primer_sets_tested$primer_set_label_long),
               selected = 4),
          radioButtons("taxo_level", "Taxonomic level:",
               c("Supergroup" = "supergroup",
                 "Class" = "class")),
          selectInput("division", "If choose Class, select Division",
               setNames(pr2_taxo$division, str_c(pr2_taxo$supergroup, pr2_taxo$division, sep = "-")),
               selected = "Chlorophyta"),          
      ),
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset' ,
        tabPanel("Primers", dataTableOutput("table_primers")),
        tabPanel("Primer sets", dataTableOutput("table_primer_sets")),
        tabPanel("Amplification - eukaryotes", plotOutput("plot_matches")),
        tabPanel("Amplification - by group", plotOutput("plot_matches_group"))
        # tabPanel("Amplification - one set", dataTableOutput("test"))
      )
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output) {

  # Primer table
  output$table_primers <- renderDataTable({
    primers[, input$show_vars_primers, drop = FALSE]
  }, escape = FALSE,                                # To keep the HTML (doi links)
  options = list(pageLength = 25)) 

  # Primer set table
  output$table_primer_sets <- renderDataTable({
    select(primer_sets[, input$show_vars_primer_sets, drop = FALSE], - primer_set_label_long)
  }, escape = FALSE,                                     # To keep the HTML (doi links)
  options = list(pageLength = 25))

  # Plot matches
  output$plot_matches <- renderPlot({
    plot_matches(input$specificity)
  }, width = 1200, height = 600, res = 96)

  # Plot matches at group level
  output$plot_matches_group <- renderPlot({
    plot_matches_group(input$primer_set_id, input$taxo_level, input$division)
  }, width = 1200, height = 600, res = 96)

  # Test
  output$test <- renderDataTable(plot_matches_sg(input$primer_set_id))    
}

# Run the shiny app -------------------------------------------------------


shinyApp(ui, server)