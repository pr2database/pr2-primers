library(shiny)
library(dplyr)
library(stringr)
# library(DT)


# Read the data -----------------------------------------------------------


primers<- readRDS("primers.rds") %>% 
  filter(str_detect(gene, "rRNA"))
  
primer_sets<- readRDS("primer_sets.rds") %>% 
  filter(str_detect(gene, "rRNA"))

pr2_match_class <- readRDS("pr2_match_18S_rRNA_mismatches_2_summary_class.rds")


# Functions ---------------------------------------------------------------

# Function to plot mismatches at supergroup level for one primer set


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
        checkboxGroupInput("show_vars_table1", "Columns to show:",
                           names(primers), selected = names(primers))
      ),
      conditionalPanel(
        'input.dataset === "Primer sets"',
        includeMarkdown("primer_sets.md"),
        checkboxGroupInput("show_vars_table2", "Columns to show:",
                           names(primer_sets), selected = names(primer_sets))
      ),
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset' ,
        tabPanel("Primers", dataTableOutput("mytable1")),
        tabPanel("Primer sets", dataTableOutput("mytable2")),
        tabPanel("Primer sets", dataTableOutput("mytable2"))
      )
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output) {

  # choose columns to display
  output$mytable1 <- renderDataTable({
    primers[, input$show_vars_table1, drop = FALSE]
  }, escape = FALSE,                                # To keep the HTML (doi links)
  options = list(pageLength = 10)) 

  # choose columns to display
  output$mytable2 <- renderDataTable({
    primer_sets[, input$show_vars_table2, drop = FALSE]
  }, escape = FALSE,                                     # To keep the HTML (doi links)
  options = list(pageLength = 10))


}

# Run the shiny app -------------------------------------------------------


shinyApp(ui, server)