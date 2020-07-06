
# Function to summarize at any level --------------------------------------

primer_summary <- function(pr2_match, taxo_level=supergroup){
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


# Function to plot matches for all primers --------------------------------


plot_matches <- function(specific_one = "general") {  
  
  df <- pr2_match_euk_long_1 %>% 
      filter(pct_category %in% c("ampli_pct", "fwd_pct", "rev_pct")) %>% 
      filter(specific == specific_one)

  g1 <- ggplot(df) + 
    geom_col(aes(x=primer_set_label_long, y=pct_seq, 
                 fill=fct_reorder(pct_category, pct_category_order)), width=.7, position = "dodge") +
    theme_bw() +
    labs(x = "Primer set",
         y = "Amplicon size (bp)",
         title = str_c("Primer type: ",specific_one),
         subtitle = "% of sequences amplified\n with 2 mismatches on each primer") +
    xlab("Primer set") + ylab("% of sequences amplified") + 
    scale_fill_manual(name = "",  values = c("ampli_pct" = "black", "fwd_pct" = "grey80","rev_pct" = "grey40"), 
                      labels=c( "Amplicons", "Primer rev", "Primer fwd")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0)) +
    ylim(0,100) +
    coord_flip() + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top", legend.box = "horizontal") 
  
# ---

  df <- pr2_match_euk_long_2 %>%  
      filter(specific == specific_one)  
  
  g2 <- ggplot(df,
            aes(x = primer_label,
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

# ---
    
  df <- pr2_match_euk %>% 
      filter(!is.nan(ampli_size_mean)) %>%  
      filter(specific == specific_one)
  
  g3 <- ggplot(df) + 
    geom_point(aes(x=primer_label, y=ampli_size_mean), colour="black") +
    geom_errorbar(aes(x=primer_label, 
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

# Function to plot mismatches at level for one primer  --------

  plot_matches_group <- function(one_primer_set = 4, taxo_level_quoted = "supergroup", taxo_name = "Chlorophyta"){
    
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
    
    df <- pr2_match
    
    taxo_level=as.symbol(taxo_level_quoted)
    
    if (taxo_level_quoted == "class") {df <- df %>% filter(division == taxo_name) } 
    
    
    
  # --- Compute statistics according to taxonomical level
    
    df2 <- df %>% 
      primer_summary(taxo_level=!!taxo_level) %>% 
      tidyr::pivot_longer(names_to="mismatch_number", 
                          values_to = "mismatch_pct", 
                          cols = contains("ampli_mismatch"),
                          names_prefix = "ampli_mismatch_") %>% 
      select(-contains("ampli_mismatch")) %>% 
      mutate(mismatch_number = str_sub(mismatch_number,1 , 1)) %>% 
      mutate(mismatch_number = str_replace(mismatch_number,"5" , "5+")) 
  
  
  # --- Plot size of amplicons
  
    g2 <- ggplot(df, aes(x= !!taxo_level, y=ampli_size)) + 
        geom_boxplot(outlier.alpha = 0.3)+
        theme_bw() + 
        theme(axis.text.y = element_blank()) +  #` Remove legend` +
        coord_flip() +
        geom_hline(yintercept = c(450,550) , linetype= c(2,3)) +
        # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") + ylab("Amplicon size (bp)") 
    
  # --- Plot location of mismatches 
  
    g3 <- ggplot() +
        geom_histogram(data = filter(df, !is.na(fwd_mismatch_primer_position_5prime)),
                       aes(x = fwd_mismatch_primer_position_5prime, fill = !!taxo_level),
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
    

    g4 <- ggplot() +
        geom_histogram(data = filter(df, !is.na(rev_mismatch_primer_position_5prime)),
                       aes(x = rev_mismatch_primer_position_5prime, fill = !!taxo_level),
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

    # --- Plot % of mismatches  
      g1 <- ggplot(df2,
                aes(x=str_c(!!taxo_level, " - n = ", n_seq),
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
    
    
    g <- (g1+g2) /(g3+g4)
    return(g)
  }