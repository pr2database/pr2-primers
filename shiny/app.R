library(shiny)
# library(DT)


primers<- readRDS("primers.rds")
primer_sets<- readRDS("primer_sets.rds") 


ui <- fluidPage(
  title = "Primer database",
  titlePanel("Primer database"),
  sidebarLayout(
    sidebarPanel(width = 2,
      conditionalPanel(
        'input.dataset === "Primers"',
        checkboxGroupInput("show_vars_table1", "Columns to show:",
                           names(primers), selected = names(primers))
      ),
      conditionalPanel(
        'input.dataset === "Primer sets"',
        checkboxGroupInput("show_vars_table2", "Columns to show:",
                           names(primer_sets), selected = names(primer_sets))
      ),
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset' ,
        tabPanel("Primers", dataTableOutput("mytable1")),
        tabPanel("Primer sets", dataTableOutput("mytable2"))
      )
    )
  )
)

server <- function(input, output) {

  # choose columns to display
  output$mytable1 <- renderDataTable({
    primers[, input$show_vars_table1, drop = FALSE]
  })

  # choose columns to display
  output$mytable2 <- renderDataTable({
    primer_sets[, input$show_vars_table2, drop = FALSE]
  })


}

shinyApp(ui, server)