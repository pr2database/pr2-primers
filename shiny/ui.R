
# User interface ----------------------------------------------------------

ui <- fluidPage(
  title = "Primer database",
  titlePanel(div(img(src='pr2_logo.png', width="80"),"The PR2 primer database")),

  # fluidRow(
  #   # --- Header
  #   img(src='pr2_logo.png', align = "left"),
  #   column(8, includeMarkdown("readme.md"))
  #          ),
  
  # --- Side bar
  sidebarLayout(
    
    # --- Side bar Panel
    sidebarPanel(width = 3,
                 
      # --- About
      
      conditionalPanel(
        'input.panel == "About"',
        includeMarkdown("about.md")
      ),
      
      # --- Panel primers
      conditionalPanel(
        'input.panel == "Primers"',
        includeMarkdown("primers.md"),
        p(),
        downloadButton("table_primers_download", "Download primers"),
        p(),
        checkboxGroupInput("show_vars_primers", "Columns to show:",
                           names(primers), selected = c("primer_id", "gene", "organelle", "direction", 
                                                        "name", "sequence", "length", "start_yeast", "specificity", 
                                                        "reference"))
      ),
      
      # --- Panel primer sets
      conditionalPanel(
        'input.panel == "Primer sets"',
        includeMarkdown("primer_sets.md"),
        p(),
        downloadButton("table_primer_sets_download", "Download primer sets"),
        p(),
        checkboxGroupInput("show_vars_primer_sets", "Columns to show:",
                           names(primer_sets), selected = c("primer_set_id", "primer_set_name", "gene", "gene_region", 
                                                            "specificity", "fwd_name", "rev_name", "amplicon_size", "reference", 
                                                            "primer_set_label_long"))
      ),
      
      # --- General amplification properties
      conditionalPanel(
        'input.panel == "Amplification - overview"',
          h3("Precomputed results for primer sets"),
          p("Against PR2 sequence database"),
          radioButtons("type", "Primer type:",
               c("General" = "general",
                 "Specific" = "specific")),
      ),
      
      # --- Amplification for one primer set and one taxonomic level
      conditionalPanel(
        'input.panel == "Amplification - details"',
          h3("Precomputed results for primer sets"),
          p("Against PR2 sequence database"),
          selectInput("primer_set_id", "Primer set",
               setNames(primer_sets_tested$primer_set_id, primer_sets_tested$primer_set_label_long),
               selected = 4)
      ),
      
      # --- Amplification for own primer set
      conditionalPanel(
        'input.panel == "Test your primer set"',
        
         useShinyFeedback(), # include shinyFeedback
         h3("Test your primer set"),
         p("Primer set is tested against the PR2 database."),
         p("Use only UIPAC characters (", strong("ACGTRYSWKMBDHVN"),")."), 
         p(strong("Length of primers:"), "between 15 and 30 bp."),
        
         textInput("fwd_seq", "Primer forward (5' -> 3')", value = "GCCAGCAVCYGCGGTAAY"),
         radioButtons("fwd_mis", "Max mismatches", inline = TRUE,  choices = c(0, 1, 2), selected = 0),
         textInput("rev_seq", "Primer reverse (5' -> 3')", value = "CCGTCAATTHCTTYAART"),
         radioButtons("rev_mis", "Max mismatches", inline = TRUE,  choices = c(0, 1, 2), selected = 0),

         actionButton("button_compute", "Run"),
        
      ),
        
      # --- Dynamic boxes for taxonomy - See https://mastering-shiny.org/action-dynamic.html
      conditionalPanel(
        'input.panel == "Amplification - details" || input.panel == "Test your primer set"', 
          hr(),
          selectInput("supergroup", "Supergroup",
               choices = c("All", unique(pr2_taxo$supergroup))),
          selectInput("division", "Division",
               choices = NULL),
          selectInput("class", "Class",
               choices = NULL)
      ),
      conditionalPanel(
        'input.panel == "Amplification - details"',  
          actionButton("button_plot", "Update plot")
      ),
      conditionalPanel(
        'input.panel == "Test your primer set"',  
          actionButton("button_plot_user", "Update plot")
      ),
        

    ),
    mainPanel(
      tabsetPanel(
        id = 'panel' ,
        
        tabPanel("About", column(8, includeMarkdown("readme.md"))),

        
        tabPanel("Primers", dataTableOutput("table_primers")),
        
        tabPanel("Primer sets", dataTableOutput("table_primer_sets")),
        
        tabPanel("Amplification - overview",
                 withSpinner(uiOutput("plot_matches_all"))),
        
        tabPanel("Amplification - details", 
                 withSpinner(uiOutput("plot_matches_one"))),
        
        tabPanel("Test your primer set" , 
                 withSpinner(uiOutput('primer_matches_user')),
                 uiOutput('primer_matches_user_graph'))
        
       # tabPanel("Amplification - one set", dataTableOutput("test"))
      )
    )
  )
)
