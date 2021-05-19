# --- Side bar Panel
side_bar <- sidebarPanel(width = 3,
             
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
                                                               "name", "sequence", "length", "start_yeast", "type", "specificity", 
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
               selectInput("kingdom_3", "Kingdom",
                           choices = c("Eukaryota", "Bacteria", "Archaea")),
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
               
               actionButton("button_compute_primer_set", "Run"),
               
             ),
             
             # --- Amplification for own primer/primer
             conditionalPanel(
               'input.panel == "Test your primer/probe"',

               useShinyFeedback(), # include shinyFeedback
               h3("Test your primer/probe"),
               p("Primer/Probe is tested against the PR2 database."),
               p("Use only UIPAC characters (", strong("ACGTRYSWKMBDHVN"),")."),
               p(strong("Length of primers:"), "between 15 and 30 bp."),

               textInput("primer_seq", "Primer/probe (5' -> 3')", value = "GCCAGCAVCYGCGGTAAY"),
               radioButtons("primer_mis", "Max mismatches", inline = TRUE,  choices = c(0, 1, 2), selected = 0),
               radioButtons("primer_type", "Type", inline = TRUE,  choices = c("primer fwd","primer rev/probe"), selected = "primer fwd"),

               actionButton("button_compute_primer", "Run"),

             ),
             
             # --- Dynamic boxes for taxonomy - See https://mastering-shiny.org/action-dynamic.html
             conditionalPanel(
               'input.panel == "Amplification - details" ||
                input.panel == "Test your primer set" ||
                input.panel == "Test your primer/probe"',
               hr(),
               selectInput("kingdom", "Kingdom",
                           choices = c("Eukaryota", "Bacteria", "Archaea")),
               selectInput("supergroup", "Supergroup",
                           choices = NULL),
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
               actionButton("button_plot_primer_set_user", "Update plot")
             ),
             conditionalPanel(
               'input.panel == "Test your primer/probe"',
               actionButton("button_plot_primer_user", "Update plot")
             ),
             
             
)
