# ---
# This script process silva.seed_v132.fasta to incorporate into pr2
# ---


# Libraries tidyr ---------------------------------------------------------

  # library("ggplot2") 
  library("dplyr")
  library("tibble")
  library("readr")
  library("stringr")


  if(any(grepl("package:dvutils", search()))) detach("package:dvutils", unload=TRUE)
  library("dvutils")

# Read pr2 files and matching files -----------------------

  pr2 <- readRDS("../shiny/data/pr2_4.12.0.rds")
  
  match_summary <- readRDS("../shiny/data/pr2_match_18S_rRNA_mismatches_2_summary.rds")
  load("../shiny/data/pr2_match_18S_rRNA_set_001_mismatches_2.rda")


# Read the silva file -----------------------

  silva <- fasta_read ("../silva_prok/silva.seed_v132.fasta")

# Clean up the silva file -----------------------
  
  silva_pr2 <- silva %>% 
    tidyr::separate(seq_name, sep = "\t", c("silva_accession","junk", "taxonomy" )) %>% 
    tidyr::separate(silva_accession, sep = "[.]", c("pr2_accession","junk2" )) %>% 
    tidyr::separate(taxonomy, sep = ";", c("kingdom","division", "class", "order", "family", "genus" )) %>% 
    mutate(supergroup = kingdom, species = str_c(genus, "_sp."),
           sequence_length = str_length(sequence)) %>% 
    select(-contains("junk")) %>% 
    filter(kingdom != "Eukaryota",
           !is.na(genus),
           genus != "    ") %>% 
    arrange(genus)
  
# Clean up the silva file -----------------------
  
  saveRDS(silva_pr2, "../shiny/data/silva_seed_132.rds")
  
# Check that the silva file can be appended to pr2-----------------------
  
  pr2_total <- bind_rows(pr2, silva_pr2)
    
    
