# ---
# Aim: Update all files necessary for secondary analyses
# Author; Daniel Vaulot
# Updated : 2020-11-03
# 
# ---

# ===========================================================
#     Initialize
# ===========================================================


# Load libraries ----------------------------------------------------------

 suppressPackageStartupMessages({

  library("dplyr")   
  library("stringr")
  library("purrr") 
  library("rio") 
   
  library(dvutils)
})

# ===========================================================
#     Read primers
# ===========================================================

file_param <- "param_pr2_primers.R"
source(file_param)
cat(readChar(file_param, 1e5))

# ===========================================================
#     Primer database
# ===========================================================


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
  mutate(doi_html = case_when(!is.na(doi) ~ str_c('<a href="https://doi.org/', doi,'">', doi, '</a>'),
                         TRUE ~ doi)) %>%  
  relocate(doi_html, .after = doi)
    

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
    mutate(doi_html = case_when(!is.na(doi) ~ str_c('<a href="https://doi.org/', doi,'">', doi, '</a>'),
                           TRUE ~ doi))  %>% 
    rename_all(funs(str_replace(., "_yeast", ""))) %>% 
    select(-remark_internal, -used_in) %>%
    relocate(fwd_id, .before = fwd_name) %>% 
    relocate(rev_id, .before = rev_name) %>%
    relocate(reference, doi, remark, .after = last_col()) %>%  
    relocate(doi_html, .after = doi) %>% 
    arrange(gene_region,  fwd_start, rev_start) 

# Save the primer files -------------------------------------------------------

export(primers, "../shiny/data/primers.rds")
export(primer_sets, "../shiny/data/primer_sets.rds")  


# ===========================================================
#     PR2 database version 4.12.
# ===========================================================


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

  file_name = str_c("../output/Table_primers_",gene_selected, ".xlsx")
  export(primer_sets, file = file_name, firstActiveRow = 2)

  
  
# ===========================================================
#     Summarize matches
# ===========================================================
  
# This is to be done only when new primer sets are computed


# Build the file for all primer sets --------------------------------------

primer_sets <- primer_sets %>%
  filter(
    gene == "18S rRNA",
    !is.na(doi),
    !str_detect(gene_region, "ITS|cloning|full")
  ) %>%
  mutate(specific = ifelse(is.na(specificity), "general", "specific"))


  pr2_match_list <- list()
  
  for (i in 1:nrow(primer_sets)) {
    skip_to_next <- FALSE
    cat("i = ", i, "\n")
    
    rda_file_label <- str_c(sprintf("_set_%03d", primer_sets$primer_set_id[i]), "_mismatches_", max_mismatch)
    file_name = str_c("../shiny/data/pr2_match_", gene_selected ,rda_file_label, ".rda")
    
    tryCatch(load(file=file_name), 
             error=function(e) {
               warning(stringr::str_c("Cannot read file: ", file_name))
               skip_to_next <<- TRUE        },
		         call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    # load(file=str_c("../shiny/data/pr2_match_", gene_selected ,rda_file_label, ".rda"))
    
    if(skip_to_next) { next } 
    pr2_match_list[[i]] <- pr2_match
     
  }

  pr2_match_final <- pr2_match_list %>% 
  reduce(bind_rows)
  
  saveRDS(pr2_match_final, file=str_c("../output/pr2_match_", gene_selected ,"_mismatches_", max_mismatch, ".rds"))


# Filter the pr2_match dataframe ------------------------------------------

  

  print(str_c("Before filtration: ", nrow(pr2_match_final)))

  primer_sets_labels <- primer_sets %>% 
    mutate(primer_set_label_short = str_c(str_sub(gene_region,1,2), 
                      sprintf("%02d", primer_set_id), 
                      str_sub(str_replace_na(specificity, replacement = ""),1,3),
                      sep=" "),
           primer_set_label_long = str_c(gene_region, 
                                         primer_set_name, "-", 
                                         str_replace_na(specificity, "general"), 
                                         sep = " ")
           ) %>% 
    # Remove the last underscore if left by itself
    mutate(primer_set_label_short = str_replace(primer_set_label_short, " $", "")) %>% 
    select(primer_set_id, 
           primer_set_label_short,
           primer_set_label_long,
           # gene_region, 
           specific,
           specificity) 
  
  pr2_match_final <- pr2_match_final%>% 
      left_join(primer_sets_labels) %>% 
    # Remove sequences for which the introns have been removed
      filter(! str_detect(pr2_accession, "_UC")) %>% 
    # Remove sequence that are shorter
      filter((str_detect(gene_region, "V4") & sequence_length>= sequence_length_min_V4) |
             (str_detect(gene_region, "V9") & sequence_length>= sequence_length_min_V9 & kingdom == "Eukaryota") |
             (str_detect(gene_region, "V9") & sequence_length>= sequence_length_min & kingdom != "Eukaryota") |
             (! str_detect(gene_region, "V4|V9") & sequence_length>= sequence_length_min)) %>% 
      mutate(mismatch_number = fwd_mismatch_number + rev_mismatch_number) %>% 
    # Only keep the selected primers
      filter(primer_set_id %in% primer_sets$primer_set_id) 
  
  print(str_c("After filtration: ", nrow(pr2_match_final)))


# Summarize ---------------------------------------------------------------

## Function ---------------------------

  primer_summary <- function(pr2_match, taxo_level){
    summary <-  pr2_match  %>% 
      group_by({{taxo_level}}, 
               gene_region, 
               primer_set_label_long, 
               primer_set_label_short,
               primer_set_id, 
               specific, 
               specificity ) %>% 
      summarize (n_seq = n(),
                 fwd_number = sum(!is.na(fwd_pos)),
                 fwd_pct = fwd_number/n_seq*100, 
                 fwd_mismatch_0_pct = sum(fwd_mismatch_number==0, na.rm = TRUE)/fwd_number*100,
                 fwd_mismatch_1_pct = sum(fwd_mismatch_number==1, na.rm = TRUE)/fwd_number*100,
                 fwd_mismatch_2_pct = sum(fwd_mismatch_number==2, na.rm = TRUE)/fwd_number*100,
                 fwd_mismatch_pos = mean(fwd_mismatch_primer_position_5prime, na.rm=TRUE),
                 rev_number = sum(!is.na(rev_pos)),
                 rev_pct = rev_number/n_seq*100,
                 rev_mismatch_0_pct = sum(rev_mismatch_number==0, na.rm = TRUE)/rev_number*100,
                 rev_mismatch_1_pct = sum(rev_mismatch_number==1, na.rm = TRUE)/rev_number*100,
                 rev_mismatch_2_pct = sum(rev_mismatch_number==2, na.rm = TRUE)/rev_number*100,
                 rev_mismatch_pos = mean(rev_mismatch_primer_position_5prime, na.rm=TRUE),
                 ampli_number = sum(!is.na(ampli_size)),
                 ampli_pct = ampli_number/n_seq*100, 
                 ampli_size_mean = mean(ampli_size, na.rm=TRUE),
                 ampli_size_sd = sd(ampli_size, na.rm=TRUE),
                 ampli_size_max = max(ampli_size, na.rm=TRUE),
                 ampli_size_min = min(ampli_size, na.rm=TRUE),
                 ampli_mismatch_0_pct = sum(mismatch_number==0, na.rm = TRUE)/n_seq*100,
                 ampli_mismatch_1_pct = sum(mismatch_number==1, na.rm = TRUE)/n_seq*100,
                 ampli_mismatch_2_pct = sum(mismatch_number==2, na.rm = TRUE)/n_seq*100,
                 ampli_mismatch_3_pct = sum(mismatch_number==3, na.rm = TRUE)/n_seq*100,
                 ampli_mismatch_4_pct = sum(mismatch_number==4, na.rm = TRUE)/n_seq*100,
                 ampli_mismatch_5_pct = sum(is.na(mismatch_number), na.rm = TRUE)/n_seq*100) %>% 
      mutate(across(contains(c("mismatch", "mean")), ~ ifelse(is.nan(.x), NA, .x)))%>% 
      mutate(across(contains(c("size")), ~ ifelse(is.infinite(.x), NA, .x)))%>% 
      ungroup()
  
  # Compute the position of the primer for nearly complete sequences
  summary_pos <-  pr2_match  %>% 
      filter(sequence_length>= sequence_length_min_V9) %>% 
      group_by({{taxo_level}}, primer_set_id ) %>% 
      summarize (fwd_pos_mean = mean(fwd_pos, na.rm=TRUE),
                 rev_pos_mean = mean(rev_pos, na.rm=TRUE)) %>% 
      ungroup()
  
  summary <- summary %>% 
    left_join(summary_pos) 
    
  }

## Summarize all eukaryotes ---------------------------

  pr2_match_summary_primer_set <- primer_summary(pr2_match_final, kingdom)

  
## Summarize per supergroup ---------------------------


  pr2_match_summary_primer_set_sg <- primer_summary(pr2_match_final, supergroup)
  
  pr2_taxo <- pr2_match_final %>% 
    select(kingdom, supergroup) %>% 
    distinct()
  
  pr2_match_summary_primer_set_sg <- pr2_match_summary_primer_set_sg %>% 
    left_join(pr2_taxo) %>% 
    relocate(kingdom,  .before = supergroup)

## Summarize per class ---------------------------
 

  pr2_match_summary_primer_set_class <- primer_summary(pr2_match_final, class) 

  pr2_taxo <- pr2_match_final %>% 
    select(kingdom, supergroup, division, class) %>% 
    distinct()
  
  pr2_match_summary_primer_set_class <- pr2_match_summary_primer_set_class %>% 
    left_join(pr2_taxo) %>% 
    relocate(kingdom, supergroup, division, .before = class)
 

## Save summaries

  saveRDS(pr2_match_summary_primer_set, file=str_c("../shiny/data/pr2_match_", gene_selected ,"_mismatches_", max_mismatch, "_summary.rds"))

  saveRDS(pr2_match_summary_primer_set, file=str_c("../output/pr2_match_", gene_selected ,"_mismatches_", max_mismatch, "_summary.rds"))
  saveRDS(pr2_match_summary_primer_set_sg, file=str_c("../output/pr2_match_", gene_selected ,"_mismatches_", max_mismatch, "_summary_sg.rds"))
  saveRDS(pr2_match_summary_primer_set_class, file=str_c("../output/pr2_match_", gene_selected ,"_mismatches_", max_mismatch, "_summary_class.rds"))  