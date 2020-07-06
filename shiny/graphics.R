
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

# Function to plot mismatches at super-group level for one primer  --------

plot_matches_sg <- function(one_primer_set = 4){

# --- Read file with Supergroup statistics

df <- pr2_match_sg_long %>%
  filter(primer_set_id == one_primer_set)

if(nrow(df) == 0)  return() 

# --- Plot % of mismatches  
g1 <- ggplot(df,
          aes(x=str_c(supergroup, " - n = ", n_seq),
              y = mismatch_pct,
              group = primer_set_id,
              fill = as.factor(mismatch_number)
             )
          ) +
  geom_col() +
  theme_bw()  +
  theme(axis.text.y = element_text(angle = 0, hjust = 0, vjust=0)) +
  labs(x = "Supergroup",
       y = "% of sequences with mismatches",
       fill = "Mismatches") +
  scale_fill_viridis_d(direction = -1, alpha=0.85) + 
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.box = "horizontal") +
  coord_flip() 

# --- Read file with all matches

rda_file_label = str_c(sprintf("_set_%03d", as.integer(one_primer_set)), "_mismatches_", max_mismatch)

# The file is loaded into pr2_match

load(file=str_c("data/pr2_match_", gene_selected ,rda_file_label, ".rda"))

df <- pr2_match %>% 
  filter(!(supergroup %in% c("Apusozoa", "Eukaryota_X")))

# --- Plot size of amplicons

g2 <- ggplot(df, aes(x= supergroup, y=ampli_size)) + 
    geom_boxplot(outlier.alpha = 0.3)+
    theme_bw()+ 
    theme(axis.text.y = element_blank()) +  #` Remove legend` +
    coord_flip()+
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("Amplicon size (bp)") 


# --- Plot location of mismatches 

g3 <- ggplot() +
    geom_histogram(data = filter(df, !is.na(fwd_mismatch_primer_position_3prime)),
                   aes(x = fwd_mismatch_primer_position_3prime, fill = supergroup),
                   binwidth = 1, 
                   stat = "bin",
                   alpha = 1)  +
    scale_fill_viridis_d(option = "inferno") +
    theme_bw() + 
    theme(legend.position="none")   +
    scale_x_continuous(minor_breaks = 0:30) +
    labs(y = "Density",
         x = "Position of mismatches from 3' end",
         title = "fwd") 


g4 <- ggplot() +
    geom_histogram(data = filter(df, !is.na(rev_mismatch_primer_position_3prime)),
                   aes(x = rev_mismatch_primer_position_3prime, fill = supergroup),
                   binwidth = 1, 
                   stat = "bin",
                   alpha = 1)  +
    scale_fill_viridis_d(option = "inferno") +
    theme_bw()   +
    scale_x_continuous(minor_breaks = 0:30) +
    labs(y = "Density",
         x = "Position of mismatches on primer from 3'end",
         title = "rev")

g <- (g1+g2) /(g3+g4)
return(g)
}