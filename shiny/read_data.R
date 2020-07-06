library(shiny)
library(markdown) # To display text boxes in md

library(dplyr)
library(stringr)
library(forcats)

library(ggplot2)
library(patchwork)
# library(DT)


# Constants ---------------------------------------------------------------

max_mismatch = 2
gene_selected = "18S_rRNA"  # Can be 18S_rRNA or 16_rRNA


# Read the data -----------------------------------------------------------

# --- Primers

primers<- readRDS("data/primers.rds") %>% 
  filter(str_detect(gene, "rRNA"))
  
primer_sets<- readRDS("data/primer_sets.rds") %>% 
  filter(str_detect(gene, "rRNA")) %>% 
  mutate(primer_set_label_long = str_c(gene_region, 
                                       primer_set_name, "-", 
                                       str_replace_na(specificity, "general"), 
                                       sep = " "))


# --- At euk level

pr2_match_euk <- readRDS("data/pr2_match_18S_rRNA_mismatches_2_summary.rds") %>% 
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

primer_sets_tested <- primer_sets %>% 
  filter(primer_set_id %in% pr2_match_euk$primer_set_id)

# --- At supergroup level

pr2_match_sg <- readRDS("data/pr2_match_18S_rRNA_mismatches_2_summary_sg.rds")


pr2_match_sg_long <- pr2_match_sg %>% 
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

pr2_match_class <- readRDS("data/pr2_match_18S_rRNA_mismatches_2_summary_class.rds")
  

