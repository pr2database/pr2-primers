
# Load libraries ----------------------------------------------------------

 suppressPackageStartupMessages({

  library("dplyr")   
  library("stringr")
   
  library(pr2database)
   
  library(dvutils)
})




# Read from primer database -----------------------------------------------

  pr2_db <- db_info("pr2_google")
  pr2_db_con <- db_connect(pr2_db)
  
  primers <- tbl(pr2_db_con, "pr2_primers") %>% 
    collect()
  primer_sets <- tbl(pr2_db_con, "pr2_primer_sets") %>% 
    collect()
  
  disconnect <- db_disconnect(pr2_db_con)


# Construct primers table -------------------------------------------------

sequence_revcomp=tibble(`sequence revcomp` = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(primers$sequence))))

primers <- primers %>% 
  bind_cols(sequence_revcomp) %>% 
  relocate(`sequence revcomp`, .after = sequence) %>% 
  relocate(gene, organelle, direction, .after = primer_id) %>%
  mutate(length = str_length(sequence)) %>% 
  relocate(length, .after = `sequence revcomp`) %>%  
  arrange(gene, start_yeast) %>% 
  mutate(doi = case_when(!is.na(doi) ~ str_c('<a href="https://doi.org/', doi,'">', doi, '</a>'),
                         TRUE ~ doi)) 
    

# Construct primer set table ----------------------------------------------


  primer_sets <- primer_sets %>% 
    left_join(select(primers, 
                     primer_id, 
                     fwd_name=name,
                     fwd_seq=sequence, 
                     fwd_start_yeast= start_yeast, 
                     fwd_end_yeast= end_yeast), 
              by = c("fwd_id" = "primer_id")) %>% 
    left_join(select(primers, 
                   primer_id, 
                   rev_name=name,
                   rev_seq=sequence, 
                   rev_start_yeast= start_yeast, 
                   rev_end_yeast= end_yeast), 
            by = c("rev_id" = "primer_id")) %>% 
    mutate(amplicon_size = rev_end_yeast - fwd_start_yeast + 1) %>% 
    rename(doi = reference_doi) %>% 
    mutate(doi = case_when(!is.na(doi) ~ str_c('<a href="https://doi.org/', doi,'">', doi, '</a>'),
                           TRUE ~ doi)) %>%  
    rename_all(funs(str_replace(., "_yeast", ""))) %>% 
    select(-remark_internal, -used_in) %>%
    relocate(fwd_id, .before = fwd_name) %>% 
    relocate(rev_id, .before = rev_name) %>%
    relocate(reference, doi, remark, .after = last_col()) %>% 
    arrange(gene_region,  fwd_start, rev_start) 

# PR2 database version 4.12.0 ---------------------------------------------

gene_selected = "18S_rRNA"
sequence_length_min = 1500


load("../../versions/4.12.0/pr2.rda")  

# Only keep sequences with different sequence_hash

pr2 <- pr2 %>%  
  filter (!str_detect(sequence,"[^ATGCU]")) %>% # Remove sequences with ambiquities
  filter(! str_detect(pr2_accession, "_UC")) %>% # Remove sequences for which the introns have been removed
  filter(gene == gene_selected )  %>%  # Keep only 18S
  filter(is.na(removed_version)) %>%  # Remove sequences that have been removed
  filter(sequence_length >= sequence_length_min) %>% # Remove sequences that are too short
  select(pr2_accession, kingdom:genus, species, sequence:sequence_hash) %>% 
  distinct(sequence_hash, .keep_all = TRUE) %>% 
  select(-sequence_hash)


  

# Save the database -------------------------------------------------------

saveRDS(primers, "../shiny/data/primers.rds")
saveRDS(primer_sets, "../shiny/data/primer_sets.rds")  
saveRDS(pr2, "../shiny/data/pr2_4.12.0.rds")  
