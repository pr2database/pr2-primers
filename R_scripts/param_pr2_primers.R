max_mismatch = 2
gene_selected = "18S_rRNA"  # Can be 18S_rRNA or 16_rRNA
rda_file_label = "_all"  # This label is added at the end of the rda file name
gene_regions = c("V4", "V9")
kingdoms = c("Eukaryota", "Bacteria", "Archaea")

# These primers should not be included (semi nested PCR or amplify very few sequences)
# primer_sets_excluded = c(73, 74, 37, 67, 91, 68) 
primer_sets_excluded = c(999) 

sequence_length_min = 1350
sequence_length_min_V9 = 1650 
sequence_length_min_V4 = 1200 