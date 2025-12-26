#!/usr/bin/env Rscript

## Visualise base-modification profiles of nuclear Trp tRNAs
## comparing all reads vs anticodon (UCA)-modified reads

library(tidyverse)

## Read input file paths from command line
args <- commandArgs(trailingOnly = TRUE)
tRNA_base_mod <- as.character(args[1])      # aggregated per-base counts (all reads)
anticodon_mod_reads <- as.character(args[2])# per-base counts for UCA reads

## Store working directory for output
dir <- getwd()

## Read and filter nuclear Trp tRNAs, excluding Arabidopsis and Balanophoraceae
trp_base_mod <- read.csv(tRNA_base_mod, header = TRUE, sep = "\t") %>% 
  filter(grepl("nuclear-Ath-Trp", GeneName)) %>%
  filter(!grepl("At", Species)) %>%
  filter(!grepl("Balano", Species)) %>% 
  select(-c(Isodecoder, Compartment))

## Label data as coming from all reads
trp_base_mod$Group <- "All reads"

## Read anticodon-modified (UCA) reads and harmonize columns
trp_base_mod <- read.csv(anticodon_mod_reads, header = TRUE, sep = "\t") %>%
  mutate(Group = "UCA reads") %>% 
  select(-c(Lib, Total, Ins)) %>% 
  filter(!grepl("At", Species)) %>%
  rbind(trp_base_mod)

## Identify tRNAs that have at least one UCA-modified read
tRNAs_wUCA_reads <- trp_base_mod %>%
  filter(Group == "UCA reads") %>%
  distinct(GeneName) %>%
  pull(GeneName)

## Keep only tRNAs present in both groups
trp_base_mod <- trp_base_mod %>%
  filter(GeneName %in% tRNAs_wUCA_reads)

## Reshape base counts to long format and classify mapping status
data2plot <- trp_base_mod %>%
  pivot_longer(cols = c(A, G, C, T, Del),
               names_to = "Base",
               values_to = "Reads") %>%
  mutate(
    MapStatus = case_when(
      Base == Reference ~ "Reference",
      Base == "Del" ~ "Deletion",
      TRUE ~ paste(Base, "Mismatch", sep = "_")
    )
  )

## Set plotting order for mapping categories
data2plot$MapStatus <- factor(
  data2plot$MapStatus,
  levels = c("Deletion", "A_Mismatch", "C_Mismatch",
             "G_Mismatch", "T_Mismatch", "Reference")
)

## Generate one PDF per species
for (i in unique(data2plot$Species)) {

  data2plot %>%
    filter(Species == i) %>%
    ggplot(aes(x = Position, y = Reads, fill = MapStatus)) +
      geom_col() +
      facet_grid(Group ~ GeneName, scales = "free_y") +
      theme_bw() +
      scale_fill_manual(values = c("black", "chartreuse4", "steelblue",
                                   "goldenrod", "firebrick", "gray70")) +
      ylab("Read Count") +
      scale_x_continuous(breaks = seq(0, 80, 5)) +
      ggtitle(i) +
      theme(
        legend.key.size = unit(0.3, "cm"),
        strip.text = element_text(size = 6, margin = margin(2,2,2,2)),
        axis.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.text = element_text(size = 5.5),
        axis.title = element_text(size = 8, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

  ## Save species-specific plot
  ggsave(paste0(dir, "/", i, "-Nuclear-Trp-base-mod.pdf"),
         device = "pdf", width = 6, height = 3)
}
