library(shiny)
library(markdown) # To display text boxes in md

library(dplyr)
library(stringr)
library(forcats)

library(ggplot2)
library(patchwork)
# library(DT)


# Read the data -----------------------------------------------------------

# --- Primers

primers<- readRDS("primers.rds") %>% 
  filter(str_detect(gene, "rRNA"))
  
primer_sets<- readRDS("primer_sets.rds") %>% 
  filter(str_detect(gene, "rRNA"))


# --- At euk level

pr2_match <- readRDS("pr2_match_18S_rRNA_mismatches_2_summary.rds")

pct_category_order <- data.frame(pct_category = c("ampli_pct","fwd_pct","rev_pct"), 
                                 pct_category_order = c(1, 3, 2))

pr2_match_long <- pr2_match %>% 
    tidyr::pivot_longer(names_to="pct_category", 
                        values_to = "pct_seq", 
                        cols = fwd_pct:ampli_pct) %>% 
    left_join(pct_category_order)

# --- At supergroup level

pr2_match_sg <- readRDS("pr2_match_18S_rRNA_mismatches_2_summary_sg.rds")


pr2_mismatches_sg <- pr2_match_sg %>% 
  tidyr::pivot_longer(names_to="mismatch_number", 
                      values_to = "mismatch_pct", 
                      cols = contains("ampli_mismatch"),
                      names_prefix = "ampli_mismatch_") %>% 
  select(-contains("ampli_mismatch")) %>% 
  mutate(mismatch_number = str_sub(mismatch_number,1 , 1)) %>% 
  mutate(mismatch_number = str_replace(mismatch_number,"5" , "5+")) %>% 
  filter(n_seq > 20) %>% 
  filter(!(supergroup %in% c("Apusozoa", "Eukaryota_X") ))

# --- At class level

# pr2_match_class <- readRDS("pr2_match_18S_rRNA_mismatches_2_summary_class.rds")
  




# Functions ---------------------------------------------------------------

# --- Function to plot matches for all primers

plot_matches <- function(specific_one = "general") {  
  df <- pr2_match_long %>% 
      filter(pct_category %in% c("ampli_pct", "fwd_pct", "rev_pct")) %>% 
      filter(specific == specific_one)

  g1 <- ggplot(df) + 
    geom_col(aes(x=primer_label, y=pct_seq, 
                 fill=fct_reorder(pct_category, pct_category_order)), width=.7, position = "dodge") +
    theme_bw() +
    labs(x = "Primer set",
         y = "Amplicon size (bp)",
         title = str_c("Primer type: ",specific_one),
         subtitle = "Percentage of sequences recovered") +
    xlab("Primer set") + ylab("% of sequences amplified") + 
    scale_fill_manual(name = "",  values = c("ampli_pct" = "black", "fwd_pct" = "grey80","rev_pct" = "grey40"), 
                      labels=c( "Amplicons", "Primer rev", "Primer fwd")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0)) +
    ylim(0,100) +
    coord_flip() + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top", legend.box = "horizontal") 
    
  df <- pr2_match %>% 
      filter(!is.nan(ampli_size_mean)) %>%  
      filter(specific == specific_one)
  
  g2 <- ggplot(df) + 
    geom_point(aes(x=primer_label, y=ampli_size_mean), colour="black") +
    geom_errorbar(aes(x=primer_label, 
                      ymax=ampli_size_mean + ampli_size_sd, 
                      ymin=ampli_size_mean - ampli_size_sd)) +
    theme_bw() +
    labs(x = "Primer set",
         y = "Amplicon size (bp)",
         title = str_c("Primer type: ",specific_one),
         subtitle = "Lines correspond to limits for Illumina 2x250 and 2x300") +
    # scale_y_continuous(breaks = (1:8)*200, limits = c(0,1500)) +
    geom_hline(yintercept = c(450,550), linetype = c(2,3) ) +
    ylim(0,850) + 
    theme(axis.text.y = element_text(angle = 0, hjust = 0)) + 
    coord_flip()  
  g1+g2
  
}

# --- Function to plot mismatches at super-group level for one primer set

plot_mismatches_sg <- function(one_primer_set){

pr2_mismatches <- pr2_match_summary_primer_set_class %>% 
  tidyr::pivot_longer(names_to="mismatch_number", 
                      values_to = "mismatch_pct", 
                      cols = contains("ampli_mismatch"),
                      names_prefix = "ampli_mismatch_") %>% 
  select(-contains("ampli_mismatch")) %>% 
  mutate(mismatch_number = str_sub(mismatch_number,1 , 1)) %>% 
  mutate(mismatch_number = str_replace(mismatch_number,"5" , "5+")) %>% 
  filter((n_seq > 20) &
         (primer_set_id == one_primer_set) & 
         !(supergroup %in% c("Apusozoa", "Eukaryota_X") ))
  
g <- ggplot(pr2_mismatches,
          aes(x=str_c(supergroup, " - n = ", n_seq),
              y = mismatch_pct,
              group = primer_set_id,
              fill = as.factor(mismatch_number)
             )
          ) +
  geom_col() +
  theme_bw()  +
  theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust=0)) +
  labs(x = "Supergroup",
       y = "% of sequences with mismatches",
       fill = "Mismatches") +
  scale_fill_viridis_d(direction = -1, alpha=0.85) + 
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.box = "horizontal") +
  coord_flip() 
}

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
        'input.dataset === "Amplification"',
          radioButtons("specificity", "Primer group:",
               c("General" = "general",
                 "Specific" = "specific")),
      ),
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset' ,
        tabPanel("Primers", dataTableOutput("table_primers")),
        tabPanel("Primer sets", dataTableOutput("table_primer_sets")),
        tabPanel("Amplification", plotOutput("plot_matches"))
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
  options = list(pageLength = 10)) 

  # Primer set table
  output$table_primer_sets <- renderDataTable({
    primer_sets[, input$show_vars_primer_sets, drop = FALSE]
  }, escape = FALSE,                                     # To keep the HTML (doi links)
  options = list(pageLength = 10))

  # Plot matches
  output$plot_matches <- renderPlot({
    plot_matches(input$specificity)
  }, width = 900, height = 600)
  
}

# Run the shiny app -------------------------------------------------------


shinyApp(ui, server)