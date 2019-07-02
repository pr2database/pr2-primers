# Initialize

# Libraries for bioinfo ----------------------------------------------------

  library("Biostrings")

# Libraries tidyr ---------------------------------------------------------

  library("ggplot2") 
  library("dplyr")   
  library("readxl")
  library("tibble")
  library("readr")
  library("purrr")
  library("forcats")
  library("lubridate")
  library("stringr")

# Libraries dvutils and pr2database -------------------------------------------------------
  if(any(grepl("package:dvutils", search()))) detach("package:dvutils", unload=TRUE)
  library("dvutils")

  if(any(grepl("package:pr2database", search()))) detach("package:pr2database", unload=TRUE)
  library("pr2database")

# Options for knitting the report -------------  
  
  library(knitr)
  library(rmdformats)
  library(kableExtra) 

  knitr::opts_chunk$set(fig.width=8, 
                        fig.height=6, 
                        eval=TRUE, 
                        cache=TRUE,
                        echo=TRUE,
                        prompt=FALSE,
                        tidy=TRUE,
                        comment=NA,
                        message=FALSE,
                        warning=FALSE)
  opts_knit$set(width=90)
  options(max.print="500")  
  options(knitr.kable.NA = '')

# Read pr2
  
data(pr2)

# Remove 
# * Any sequence not in PR2
# * Any sequence with ambuiguities

pr2_active  <- pr2 %>% filter (is.na(removed_version)) %>% 
                       filter (!str_detect(sequence,"[^ATGCU]"))


pr2_db <- db_info("pr2_local")
pr2_db_con <- db_connect(pr2_db)
primers <- tbl(pr2_db_con, "pr2_primers") %>% collect()
primer_sets <- tbl(pr2_db_con, "pr2_primer_sets") %>% collect()
disconnect <- db_disconnect(pr2_db_con)

gene_regions = c("V4", "V4-specific", "V9", "Universal")

# Just keep the selected primers (V4, V9 etc..)
primer_sets <- primer_sets %>% 
    filter(gene_region %in%  gene_regions)

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
  mutate(length_yeast = rev_end_yeast - fwd_start_yeast + 1) %>% 
  select(primer_set_id:specificity, contains("fwd"), contains("rev"),length_yeast, used_in:remark)

  
# Run loops through primers  
  
pr2_match <- list()
max_mismatch <- 0

i = 1

 for (i in 1:nrow(primer_sets)) {

 print(sprintf("Primer set : %d", i))

 gene_region <- primer_sets$gene_region[[i]]  
 primer_set_id <- as.integer(primer_sets$primer_set_id[[i]])
 primer_label = str_c(str_sub(primer_sets$gene_region[[i]],1,2), 
                      str_sub(str_replace_na(primer_sets$specificity[[i]], replacement = ""),1,3), 
                      sprintf("%02d", primer_sets$primer_set_id[[i]]), sep="_")
 fwd <-  DNAString(primer_sets$fwd_seq[[i]])
 rev <-  DNAString(primer_sets$rev_seq[[i]])
 rev <-  reverseComplement(rev)

 seq <- DNAStringSet(pr2_active$sequence)
 names(seq) <- pr2_active$pr2_accession
 
 fwd_pos <- vmatchPattern(fwd, seq, max.mismatch=max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
 rev_pos <- vmatchPattern(rev, seq, max.mismatch=max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
 
# Need to unlist the MIndex position
 fwd_pos_list <- unlist(fwd_pos)
 rev_pos_list <- unlist(rev_pos)
 
# Create data frames from the list
 fwd_matches <- data.frame( fwd_pos=fwd_pos_list@start, pr2_accession =fwd_pos_list@NAMES)
 rev_matches <- data.frame( rev_pos=rev_pos_list@start +rev_pos_list@width-1 , pr2_accession =rev_pos_list@NAMES)

# Cheking whether some sequences have more than one matches fwd
 fwd_matches_pb <- fwd_matches %>% group_by(pr2_accession) %>% 
                                   summarize(n_matches=n()) %>% 
                                   filter(n_matches >1 )
 sprintf("Number of sequences with more than one match fwd : %d",nrow(fwd_matches_pb))
 
# Cheking whether some sequences have more than one matches rev 
 rev_matches_pb <- rev_matches %>% group_by(pr2_accession) %>% 
                                   summarize(n_matches=n()) %>% 
                                   filter(n_matches >1 )
 sprintf("Number of sequences with more than one match rev : %d",nrow(rev_matches_pb)) 

# If several matches, take the first match for fwd and the last match for the reverse
 fwd_matches_unique <- fwd_matches %>%  group_by(pr2_accession) %>%
                                        summarise(fwd_pos = min(fwd_pos))
 rev_matches_unique <- rev_matches %>%  group_by(pr2_accession) %>%
                                        summarise(rev_pos = max(rev_pos))
 
# Merge the fwd and reverse position, compute the length of the amplicon and check it is bigger than sum of the lengths of the two primers
 pr2_match[[i]] <- select(pr2_active, pr2_accession, kingdom:genus, species, sequence_length) %>% 
    mutate(gene_region=gene_region, primer_set_id= primer_set_id, primer_label = primer_label) %>% 
    left_join(fwd_matches_unique) %>% 
    left_join(rev_matches_unique) %>%  
    mutate(ampli_size = case_when( !is.na(fwd_pos) & 
                                   !is.na(rev_pos)  ~ rev_pos - fwd_pos + 1)) %>% 
   # Must use _NA_real_ and not NA alone because will cause an error....
    mutate(ampli_size = case_when( ampli_size < (length(fwd) + length(rev)) ~  NA_real_,  
                                   TRUE ~ ampli_size)) 

}

# Use reduce to collapse the list of dataframes into a single data frame (a very very elegant way)
 pr2_match_final <- pr2_match %>% 
   reduce(bind_rows)


## Save the big data file

  save(pr2_match_final, file="pr2_match.rda")
