
# Deploying to Shiny server -----------------------------------------------

# When deploying on ShinyApps this before at the console type at the console the following line
#        options(repos = BiocManager::repositories())
# See : https://aarthiramakrishnan.com/2018/01/09/deploying-shiny-apps-using-bioconductor.html
# To deploy app: rsconnect::deployApp("~/projects/shiny/app1", appName = "myapp", appTitle = "My Application")


# Libraries ---------------------------------------------------------------


library(shiny)
library(shinyFeedback) # For user feedback
library(shinycssloaders) # For the spinning wheel


library(markdown) # To display text boxes in md

library(dplyr)
library(stringr)
library(forcats)
library(rio) # For export/import

library(ggplot2)
library(patchwork)


library(Biostrings)
# library(DT)


# Constants ---------------------------------------------------------------

max_mismatch = 2
gene_selected = "18S_rRNA"  # Can be 18S_rRNA or 16_rRNA

sequence_length_min = 1500
sequence_length_min_V9 = 1650 
sequence_length_min_V4 = 1200 

taxo_levels <- c("kingdom", "supergroup", "division", "class", "genus")

# Read the data -----------------------------------------------------------

# --- Primers

  primers<- import("data/primers.rds") %>% 
    filter(str_detect(gene, "rRNA"))
    
  primer_sets<- import("data/primer_sets.rds") %>% 
    filter(str_detect(gene, "rRNA"))  %>% 
    mutate(primer_set_label_long = str_c(sprintf("%02d",  primer_set_id), "-", 
                                         gene_region, 
                                         primer_set_name, "-", 
                                         str_replace_na(specificity, "general"), 
                                         sep = " "))


# --- At euk level

  pr2_match_euk <- import("data/pr2_match_18S_rRNA_mismatches_2_summary.rds") %>% 
    left_join(select(primer_sets, primer_set_id, primer_set_label_long))
  
  pct_category_order <- data.frame(pct_category = c("ampli_pct","fwd_pct","rev_pct"), 
                                   pct_category_order = c(1, 3, 2))

# ------ This is for the first graph with matches fwd and reverse

  pr2_match_euk_long_1 <- pr2_match_euk %>% 
      tidyr::pivot_longer(names_to="pct_category", 
                          values_to = "pct_seq", 
                          cols = fwd_pct:ampli_pct) %>% 
      left_join(pct_category_order)

# ------ This is for the second graph with number of mismatches

  pr2_match_euk_long_2 <- pr2_match_euk %>% 
    tidyr::pivot_longer(names_to="mismatch_number", 
                        values_to = "mismatch_pct", 
                        cols = contains("ampli_mismatch"),
                        names_prefix = "ampli_mismatch_") %>% 
    select(-contains("ampli_mismatch")) %>% 
    mutate(mismatch_number = str_sub(mismatch_number,1 , 1)) %>% 
    mutate(mismatch_number = str_replace(mismatch_number,"5" , "5+")) 

  # ------ List of primer tested
  
  primer_sets_tested <- primer_sets %>% 
    filter(primer_set_id %in% pr2_match_euk$primer_set_id)

# --- Read pr2
  pr2 <- import("data/pr2_4.12.0.rds")

# --- Taxonomic levels

  # load("data/pr2_match_18S_rRNA_set_008_mismatches_2.rda")

  pr2_taxo <- pr2 %>% 
      select(supergroup, division, class, order, family, genus) %>% 
      arrange(supergroup, division, class, order, family, genus) %>% 
      distinct()