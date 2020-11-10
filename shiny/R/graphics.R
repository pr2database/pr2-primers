
# Function to plot matches for all primers --------------------------------


plot_matches <- function( kingdom_one = "Eukaryota", type = "general") {  

  # --- Plot number of sequence amplified with fwd, rev and both
  
  df <- pr2_match_summary_long_1 %>% 
      filter(pct_category %in% c("ampli_pct", "fwd_pct", "rev_pct")) %>% 
      filter(specific == type) %>% 
      filter(kingdom == kingdom_one)
  
  # print(kingdom_one)

  g1 <- ggplot(df) + 
    geom_col(aes(x=primer_set_label_long, y=pct_seq, 
                 fill=fct_reorder(pct_category, pct_category_order)), width=.7, position = "dodge") +
    theme_bw() +
    labs(x = "Primer set",
         y = "% of sequences amplified",
         title = str_c("Primer type: ",type),
         subtitle = "% of sequences amplified\n with 2 mismatches on each primer") +
    scale_fill_manual(name = "",  values = c("ampli_pct" = "black", "fwd_pct" = "blue","rev_pct" = "red"), 
                      labels=c( "Amplicons", "Primer rev", "Primer fwd")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0)) +
    ylim(0,100) +
    coord_flip() + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top", legend.box = "horizontal") 
  
  
  # --- Plot number of mismatches

  df <- pr2_match_summary_long_2 %>%  
      filter(specific == type) %>% 
      filter(kingdom == kingdom_one) 
  
  g2 <- ggplot(df,
            aes(x = primer_set_label_long,
                y = mismatch_pct,
                group = primer_set_id,
                fill = as.factor(mismatch_number)
               )
            ) +
    geom_col() +
    theme_bw()  +
    theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust=0)) +
    labs(x = "",
         y = "% of sequences with mismatches",
         fill = "Mismatches") +
    scale_fill_viridis_d(direction = -1, alpha=0.85) + 
    guides(fill = guide_legend(nrow = 1, label.position = "top")) + 
    theme(axis.text.y = element_blank()) +  #` Remove legend`
    theme(legend.position = "top", legend.box = "horizontal",
          legend.title=element_text(size=9), 
          legend.text=element_text(size=8)) +
    coord_flip() 

  # --- Plot amplicon size
    
  df <- pr2_match_summary %>% 
      filter(!is.nan(ampli_size_mean)) %>%  
      filter(specific == type) %>% 
      filter(kingdom == kingdom_one)
  
  g3 <- ggplot(df) + 
    geom_point(aes(x=primer_set_label_long, y=ampli_size_mean), colour="black") +
    geom_errorbar(aes(x=primer_set_label_long, 
                      ymax=ampli_size_mean + ampli_size_sd, 
                      ymin=ampli_size_mean - ampli_size_sd)) +
    theme_bw() +
    labs(x = "",
         y = "Amplicon size (bp)",
         title = "",
         subtitle = "Lines correspond to limits\n for Illumina 2x250 and 2x300") +
    # scale_y_continuous(breaks = (1:8)*200, limits = c(0,1500)) +
    geom_hline(yintercept = c(450,550), linetype = c(2,3) ) +
    ylim(0,850) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(axis.text.y = element_blank()) +  #` Remove legend`
    coord_flip()  

    g1+g2+g3 + 
      plot_layout(widths = c(1, 1, 1))
  
}




# Function to plot mismatches at a given taxonomic level for one primer set-----

  plot_matches_detailed_taxa <- function(one_primer_set = 4, taxo_level_quoted = "supergroup", taxo_name = "Chlorophyta"){
    
    print(cat("Taxo level : ", taxo_level_quoted))
    print(cat("Taxo name : ", taxo_name))
    
    one_primer_set_data = filter(primer_sets, primer_set_id == one_primer_set)
    fwd_name = one_primer_set_data$fwd_name
    rev_name = one_primer_set_data$rev_name        
    fwd_length = str_length(one_primer_set_data$fwd_seq)
    rev_length = str_length(one_primer_set_data$rev_seq)
    
  # --- Read file with all matches
  
    rda_file_label = str_c(sprintf("_set_%03d", as.integer(one_primer_set)), "_mismatches_", max_mismatch)
    
    # The rda file is loaded into pr2_match
    
    tryCatch(load(file=str_c("data/pr2_match_", gene_selected ,rda_file_label, ".rda")), 
             error=function(e) warning(stringr::str_c("Cannot read file: ", GB_file),
		         call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)) 

    df <- pr2_match %>% 
        filter(! str_detect(pr2_accession, "_UC")) %>%  # Remove sequences for which the introns have been removed
        filter(sequence_length >= sequence_length_min)  # Remove sequences that are too short
        
    rm(pr2_match) # To save memory 
    
  #  --- Compute summary at kingdom level
    
    summary_kingdom <- df  %>%      
      group_by(kingdom) %>% 
      summarise(pct_fwd = sum(!is.na(fwd_pos))/n()*100,
                pct_rev = sum(!is.na(rev_pos))/n()*100,
                pct_amplified = sum(!is.na(ampli_size))/n()*100,
                mean_amplicon_size = mean(ampli_size, na.rm = TRUE)) %>%
      tidyr::pivot_longer(cols = pct_fwd:mean_amplicon_size, names_to = "Parameter", values_to = "Values") %>%
      mutate(Parameter = str_replace_all(Parameter,c("_" = " ",
                                                     "pct" = "% sequences",
                                                     "fwd" = "matching forward primer",
                                                     "rev" = "matching reverse primer")))
      
    
  #  --- Filter by taxo level
    
    taxo_level=as.symbol(taxo_level_quoted)
    taxo_level_below_quoted = taxo_levels[which(taxo_levels == taxo_level_quoted) + 1]
    taxo_level_below=as.symbol(taxo_level_below_quoted)
    
    df <- df %>% 
      filter(!!taxo_level == taxo_name)  
    
    
  # --- Compute statistics according to taxonomical level
    
    df2 <- df %>% 
      primer_summary_detailed(taxo_level=!!taxo_level_below) %>% 
      tidyr::pivot_longer(names_to="mismatch_number", 
                          values_to = "mismatch_pct", 
                          cols = contains("ampli_mismatch"),
                          names_prefix = "ampli_mismatch_") %>% 
      select(-contains("ampli_mismatch")) %>% 
      mutate(mismatch_number = str_sub(mismatch_number,1 , 1)) %>% 
      mutate(mismatch_number = str_replace(mismatch_number,"5" , "5+")) 

  # --- Computer number of taxa  to adjust the height of the graph
    n_taxa = length(unique(pull(df2, !!taxo_level_below)))
    n_taxa = max(4, n_taxa)  # If one taxa, too small

  
  # --- Plot size of amplicons
  
    g2 <- ggplot(df, aes(x= !!taxo_level_below, y=ampli_size)) + 
        geom_boxplot(outlier.alpha = 0.3)+
        theme_bw() + 
        theme(axis.text.y = element_blank()) +  #` Remove legend` +
        coord_flip() +
        geom_hline(yintercept = c(450,550) , linetype= c(2,3)) +
        # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") + ylab("Amplicon size (bp)") 
    print("after g2")

    
  # --- Plot location of mismatches fwd primer 
  
    g3 <- ggplot() +
        geom_histogram(data = filter(df, !is.na(fwd_mismatch_primer_position_5prime)),
                       aes(x = fwd_mismatch_primer_position_5prime, fill = !!taxo_level_below),
                       binwidth = 1, 
                       stat = "bin",
                       alpha = 1)  +
        scale_fill_viridis_d(option = "magma") +
        theme_bw() + 
        theme(legend.position="none")   +
        scale_x_continuous(minor_breaks = 0:30, limits = c(0, fwd_length+1) ) +
        labs(y = "Density",
             x = "Position of mismatches from 5' end",
             title = str_c("fwd:", fwd_name, fwd_length, "bp", sep=" ")
             ) 
    print("after g3")
    
  # --- Plot location of mismatches rev primer 

    g4 <- ggplot() +
        geom_histogram(data = filter(df, !is.na(rev_mismatch_primer_position_5prime)),
                       aes(x = rev_mismatch_primer_position_5prime, fill = !!taxo_level_below),
                       binwidth = 1, 
                       stat = "bin",
                       alpha = 1)  +
        scale_fill_viridis_d(option = "magma") +
        theme_bw()   +
        scale_x_continuous(minor_breaks = 0:30, limits = c(0, rev_length+1)) +
        labs(y = "Density",
             x = "Position of mismatches on primer from 5' end",
             title = str_c("rev:", rev_name, rev_length, "bp", sep=" ")
             ) 
    print("after g4")

    # --- Plot % of mismatches  
    g1 <- ggplot(df2,
              aes(x=str_c(!!taxo_level_below, " - n = ", n_seq),
                  y = mismatch_pct,
                  fill = fct_rev(mismatch_number)
                 )
              ) +
      geom_col() +
      theme_bw()  +
      theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust=0)) +
      labs(x = "",
           y = "% of sequences with mismatches",
           fill = "Mismatches") +
      scale_fill_viridis_d(direction = 1, alpha=0.85) + 
      guides(fill = guide_legend(nrow = 1)) +
      theme(legend.position = "top", legend.box = "horizontal") +
      coord_flip()  
    print("after g1")
    
    gg <- (g3+g4) /(g1+g2) + plot_layout(heights = c(1, 0.3 + n_taxa/15))
    
    return(list(summary_kingdom = summary_kingdom, gg = gg, n_taxa = n_taxa))
  }
  
  
# Function to plot mismatches at a given taxonomic level for user primer match  -----

  plot_matches_simple_taxa <- function(df, taxo_level_quoted = "supergroup", taxo_name = "Chlorophyta"){
    

    taxo_level=as.symbol(taxo_level_quoted)
    taxo_level_below_quoted = taxo_levels[which(taxo_levels == taxo_level_quoted) + 1]
    taxo_level_below=as.symbol(taxo_level_below_quoted)
    
    df <- df %>% 
      filter(!!taxo_level == taxo_name)  
    
    
  # --- Compute statistics according to taxonomical level
    pct_category_order <- data.frame(pct_category = c("ampli_pct","fwd_pct","rev_pct"), 
                                     pct_category_order = c(1, 3, 2))

    
    df2 <- df %>% 
      primer_summary_simple(taxo_level=!!taxo_level_below) %>% 
      tidyr::pivot_longer(cols = ends_with("_pct"), 
                          names_to = "pct_category", 
                          values_to = "pct_seq") %>% 
      left_join(pct_category_order)  
    
    # print(df2)
      

  # --- Computer number of taxa  to adjust the height of the graph
    n_taxa = length(unique(pull(df2, !!taxo_level_below)))
    n_taxa = max(4, n_taxa) # If only 1 taxa, too small...
    
    # print(n_taxa)
  
  # --- Plot size of amplicons
  
    g2 <- ggplot(df, aes(x= !!taxo_level_below, y=ampli_size)) + 
        geom_boxplot(outlier.alpha = 0.3)+
        theme_bw() + 
        theme(axis.text.y = element_blank()) +  #` Remove legend` +
        coord_flip() +
        geom_hline(yintercept = c(450,550) , linetype= c(2,3)) +
        # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") + ylab("Amplicon size (bp)") 
    print("after g2")
    
    # --- Plot % of mismatches  

   g1 <- ggplot(df2) + 
    geom_col(aes(x=str_c(!!taxo_level_below, " - n = ", n_seq), 
                 y=pct_seq, 
                 fill=fct_reorder(pct_category, pct_category_order)), width=.7, position = "dodge") +
    theme_bw() +
    xlab("") + ylab("% of sequences amplified") + 
    scale_fill_manual(name = "",  values = c("ampli_pct" = "black", "fwd_pct" = "grey80","rev_pct" = "grey40"), 
                      labels=c( "Amplicons", "Primer rev", "Primer fwd")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0)) +
    ylim(0,100) +
    coord_flip() + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top", legend.box = "horizontal") 

    
    
    gg <- g1+g2 
    
    return(list(gg = gg, n_taxa = n_taxa))
  }  

    print("graphics.R done")