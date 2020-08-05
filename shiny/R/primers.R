# Function to check primers

primer_check <- function(primer){
      primer <- str_to_upper(primer)
      ((nchar(primer) <= 30) &
       (nchar(primer) >= 15) &  
       (str_detect(primer, "[^ACGTRYSWKMBDHVN]", negate = TRUE)))   
}


# Function to summarize at any level --------------------------------------

primer_summary_detailed <- function(pr2_match, taxo_level=supergroup){
  # debug
  # load("data/pr2_match_18S_rRNA_set_004_mismatches_2.rda")
  # taxo_level=supergroup 
  
  summary <-  pr2_match  %>%  
    mutate(mismatch_number = fwd_mismatch_number + rev_mismatch_number) %>% 
    group_by(across({{taxo_level}})) %>% 
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
    mutate(across(contains(c("size")), ~ ifelse(is.infinite(.x), NA, .x)))
  
  return(summary)
  
}

# Function to summarize at any level --------------------------------------

primer_summary_simple <- function(pr2_match, taxo_level=supergroup){
  # debug
  # load("data/pr2_match_18S_rRNA_set_004_mismatches_2.rda")
  # taxo_level=supergroup 
  
  summary <-  pr2_match  %>%  
    group_by(across({{taxo_level}}))  %>% 
    summarise(fwd_pct = sum(!is.na(fwd_pos))/n()*100,
              rev_pct = sum(!is.na(rev_pos))/n()*100,
              ampli_pct = sum(!is.na(ampli_size))/n()*100,
              n_seq = n()) 

  return(summary)
  
}


# Function to compute mismatches ------------------------------------------

primer_match <- function(fwd_seq="GCCAGCAVCYGCGGTAAY", rev_seq ="CCGTCAATTHCTTYAART", 
                         fwd_max_mismatch = 0 , rev_max_mismatch = 0){
  
   print("Start primer matching")
   print(fwd_max_mismatch)
   print(rev_max_mismatch)
  

   fwd <-  DNAString(fwd_seq)
   rev <-  DNAString(rev_seq)
   rev <-  reverseComplement(rev)
  
   seq <- DNAStringSet(pr2$sequence)
   names(seq) <- pr2$pr2_accession
   
   fwd_pos <- vmatchPattern(fwd, seq, max.mismatch=fwd_max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
   rev_pos <- vmatchPattern(rev, seq, max.mismatch=rev_max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
   
  # Need to unlist the MIndex position
   fwd_pos_list <- unlist(fwd_pos)
   rev_pos_list <- unlist(rev_pos)
   
  # Use mismatch to find the position of the mismatch... 
  # https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
   
  # Create data frames from the list
   fwd_matches <- data.frame( fwd_pos=fwd_pos_list@start, pr2_accession =fwd_pos_list@NAMES)
   rev_matches <- data.frame( rev_pos=rev_pos_list@start +rev_pos_list@width-1 , pr2_accession =rev_pos_list@NAMES)
  
  # Cheking whether some sequences have more than one matches fwd
   fwd_matches_pb <- fwd_matches %>% group_by(pr2_accession) %>% 
                                     summarize(n_matches=n()) %>% 
                                     filter(n_matches >1 )
   print(sprintf("Number of sequences with more than one match fwd : %d",nrow(fwd_matches_pb)))
   
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
   pr2_match <- pr2 %>% 
      select(pr2_accession, kingdom:genus, species) %>%  
      left_join(fwd_matches_unique) %>% 
      left_join(rev_matches_unique) %>%  
      mutate(ampli_size = case_when( !is.na(fwd_pos) & 
                                     !is.na(rev_pos)  ~ rev_pos - fwd_pos + 1)) %>% 
     # Must use _NA_real_ and not NA alone because will cause an error....
      mutate(ampli_size = case_when( ampli_size < (length(fwd) + length(rev)) ~  NA_real_,  
                                     TRUE ~ ampli_size))
  
   return(pr2_match)
 }