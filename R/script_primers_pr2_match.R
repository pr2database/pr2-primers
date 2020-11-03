# ---
# Aim: Compute matches for one set of primer
# Author: D. Vaulot
# This script is much slower but find # of mismatches and location of mismatches closer to 3' end
# ---

## Programing Notes

#     * Use Biostrings
# 
# Accessor methods : In the code snippets below, x is an MIndex object.  
# 
#     * length(x): The number of patterns that matches are stored for.
#     * names(x): The names of the patterns that matches are stored for.
#     * startIndex(x): A list containing the starting positions of the matches for each pattern.
#     * endIndex(x): A list containing the ending positions of the matches for each pattern.
#     * elementNROWS(x): An integer vector containing the number of matches for each pattern.
#     * x[[i]]: Extract the matches for the i-th pattern as an IRanges object.
#     * unlist(x, recursive=TRUE, use.names=TRUE): Return all the matches in a single IRanges object. recursive and use.names are ignored.
#     * extractAllMatches(subject, mindex): Return all the matches in a single XStringViews object.
# 
# One could also use another function which does not give the position
#     * match_fwd <- vcountPattern(fwd, seq,max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")

# ===========================================================================

# Initialize

 suppressPackageStartupMessages({

   
# Libraries for bioinfo ----------------------------------------------------
  library("Biostrings")

# Libraries tidyr ---------------------------------------------------------

  # library("ggplot2") 
  library("dplyr")   
  library("readxl")
  library("tibble")
  library("readr")
  library("purrr")
  library("forcats")
  library("lubridate")
  library("stringr")
})

# Libraries for line arguments ----------------------------------------------------

  suppressPackageStartupMessages(library("optparse"))

# Parameters ---------------------------------------------

##First read in the arguments listed at the command line
  option_list = list(
    make_option(c("-m", "--mismatch"), type="integer", default=2, 
                help="Number of maximum mimsatches for each primer [default= %default]", metavar="number"),
  	make_option(c("-s", "--set"), type="integer", default=36, 
                help="Primer set [default= %default]", metavar="number"),
  	make_option(c("-g", "--gene"), type="character", default="18S_rRNA", 
                help="Gene selected, can be 18S_rRNA or 16S_rRNA [default= %default]", metavar="character"),
  	make_option(c("-t", "--test"), type="character", default=FALSE, 
                help="Test [default= %default]", metavar="logical")
  ) 
   
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  max_mismatch = opt$mismatch
  set_id  = opt$set
  gene_selected = opt$gene  # Can be 18S_rRNA or 16_rRNA
  
  # This label is added at the end of the rda file name
  rda_file_label = str_c(sprintf("_set_%03d", set_id), "_mismatches_", max_mismatch)



# Function match_primer ------------------------------------------------------------

match_primer <- function(primer, seq, max_mismatch, direction = "forward") {
  # For testing
  
  # primer = fwd <-  DNAString(primer_sets$fwd_seq[[25]])
  # seq = pr2_active$sequence[1]
  # max_mismatch = 5
  # direction = "forward"
  
  seq <- DNAString(seq) 
  matches <- matchPattern(primer, seq, max.mismatch=max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
  # print(matches)

  start <- NA_integer_
  mismatch_number <- NA_integer_
  mismatch_primer_position_5prime <- NA
  matches_number <- NA_integer_
  if (length(matches) > 0) {
    matches_number <- length(matches)
    # If there are several matches try to keep the one with lower number of mismatches or the one most upstream
    for (i in 1:matches_number) {
      start[i] <- matches@ranges@start[i]  # Take only first match
      mismatches <- unlist(mismatch(primer, matches[i], fixed=FALSE))
      mismatch_number[i] = length(mismatches)
      if(mismatch_number[i] > 0) {
      mismatch_primer_position_5prime[i] = ifelse(direction == "forward", 
                                                  max(mismatches), 
                                                  length(primer) - min(mismatches) + 1) 
      }
    i_min = which.min(mismatch_number)
    start = start[i_min]
    mismatch_number = mismatch_number[i_min]
    mismatch_primer_position_5prime = mismatch_primer_position_5prime[i_min]
    }  
  }
  
  return(data.frame(start = start, 
              matches_number = matches_number,       
              mismatch_number = mismatch_number,  
              mismatch_primer_position_5prime = mismatch_primer_position_5prime))
}

# ---------------------------------------------------------------------------------------------



  
# Read the primer file -----------------------

  primer_sets <- readRDS("input/primer_sets.rds")
  primer_sets <- primer_sets %>% 
  filter(primer_set_id == set_id)
  

# Read pr2 and silva seed -----------------------------------

  # From R object <- Use for Geisen paper
  file_name = str_c("input/pr2_4.12.rda")
  load(file_name)
  
  silva_pr2 <- readRDS("input/silva_seed_132.rds")

  # From database for other analysis
  # pr2_active  <- pr2_read()

# Remove 
# * Any sequence not in PR2
# * Any sequence with ambuiguities
#  * Choose 18S or 16S

pr2_active  <- pr2 %>% 
  # filter (is.na(removed_version)) %>% 
  filter (!str_detect(sequence,"[^ATGCU]")) %>% 
  filter(gene == gene_selected ) 

pr2_active <- bind_rows(pr2_active, silva_pr2) 

# For tests
  if (opt$test) pr2_active <- pr2_active[1:100,]
# pr2_active <- pr2_active[1,]
  
# ---------------------------------------------------------------------------------------------
  
i = 1  
  
 print(sprintf("Primer set : %d", primer_sets$primer_set_id[[i]]))

 gene_region <- primer_sets$gene_region[[i]]  
 primer_set_id <- as.integer(primer_sets$primer_set_id[[i]])
 fwd <-  DNAString(primer_sets$fwd_seq[[i]])
 rev <-  DNAString(primer_sets$rev_seq[[i]])
 rev <-  reverseComplement(rev)
 rev_length = length(rev)

 # seq <- DNAStringSet(pr2_active$sequence)
 # names(seq) <- pr2_active$pr2_accession
 
 # fwd_pos <- vmatchPattern(fwd, seq, max.mismatch=max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")
 # rev_pos <- vmatchPattern(rev, seq, max.mismatch=max_mismatch, min.mismatch=0, with.indels=FALSE, fixed=FALSE, algorithm="auto")

  fwd_pos <- map_df(pr2_active$sequence, ~ match_primer(fwd, .x, max_mismatch, "forward")) %>% 
    mutate(pos = start) %>% 
    select(-start) %>% 
    rename_with(~ str_c("fwd_",.), everything())
  rev_pos <- map_df(pr2_active$sequence, ~ match_primer(rev, .x, max_mismatch, "reverse")) %>% 
    mutate(pos = start + rev_length - 1) %>% 
    select(-start) %>%  
    rename_with(~ str_c("rev_",.), everything())
  
  pr2_match <- pr2_active %>% 
    select(pr2_accession, kingdom:genus, species, sequence_length, sequence) %>% 
    mutate(primer_set_id= primer_set_id) %>% 
    bind_cols(fwd_pos) %>% 
    bind_cols(rev_pos) %>%  
    mutate(ampli_size = case_when( !is.na(fwd_pos) & 
                                   !is.na(rev_pos)  ~ rev_pos - fwd_pos + 1)) %>% 
   # Must use _NA_real_ and not NA alone because will cause an error....
    mutate(ampli_size = case_when( ampli_size < (length(fwd) + length(rev)) ~  NA_real_,  
                                   TRUE ~ ampli_size),
   # The next line remove sequence that do not have the canonical sequence GGATCA at the end of the the V9 region
   # and are likely to be too short
           keep_sequence = case_when(gene_region == "V9" & 
                                       !is.na(fwd_pos) & 
                                       !str_detect(str_sub(sequence, start=fwd_pos), "GGATC[AT]") ~ FALSE,
                                     TRUE ~ TRUE)) %>% 
    select(-sequence) %>% 
    filter(keep_sequence == TRUE) %>% 
    select(-keep_sequence)


## Save the data file

  save(pr2_match, file=str_c("output/pr2_match_", gene_selected ,rda_file_label, ".rda"))
  

