
# Goal: Compute position of primers on yeast


# Libraries dvutils and pr2database -------------------------------------------------------

  library("dplyr")   
  library(rio)
  library("purrr")
  library("stringr")
  library("here")

  if(any(grepl("package:dvutils", search()))) detach("package:dvutils", unload=TRUE)
  library("dvutils")

  if(any(grepl("package:pr2database", search()))) detach("package:pr2database", unload=TRUE)
  library("pr2database")

# Read primers database ---------------------------------------------------

pr2_db <- db_info("pr2_google")
pr2_db_con <- db_connect(pr2_db)
primers <- tbl(pr2_db_con, "pr2_primers") %>% 
  collect() %>% 
  filter(str_detect(gene, "18S"))
disconnect <- db_disconnect(pr2_db_con)


# Reference sequence ------------------------------------------------------

ref_seq <- pr2 %>% 
  filter(pr2_accession == "FU970071.1.1799_U") 


# Matches the position of the primers on the yeast sequence ---------------

#  - use map2_dfr to get an data frame on output and not a list
#  - use ~ when defining the function

primers_pos <- map2_dfr(primers$sequence,
                        primers$direction, 
                     ~ get_primer_position(.x, ref_seq$sequence, orientation =.y, mismatches = 3))

primers <- bind_cols(primers, primers_pos) %>% 
  filter(!is.na(start)) 

primers <- primers %>% 
  filter((start_yeast != start) | is.na(start_yeast)) 

export(primers, 
       file = here(str_c("R_paper/output/primers_matches_yeast_", lubridate::today(), ".xlsx")), 
       firstActiveRow = 2)
