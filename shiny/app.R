library(shiny)
library(dplyr)
library(stringr)
# library(DT)


primers<- readRDS("primers.rds") %>% 
  filter(str_detect(gene, "rRNA"))
  
primer_sets<- readRDS("primer_sets.rds") %>% 
  filter(str_detect(gene, "rRNA"))


ui <- fluidPage(
  title = "Primer databases",
  titlePanel("Primer database for protis rRNA genes"),
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