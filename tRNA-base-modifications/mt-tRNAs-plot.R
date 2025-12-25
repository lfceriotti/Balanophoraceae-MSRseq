#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
file = as.character(args[1])       

dir <- getwd()

## Load tRNA base modification data
data <- read.csv(paste(dir, "/", file, sep = ""), sep = "\t", header = T)

## Select tRNAs of interest
data_filtered <- data %>% filter(grepl("mitochondrial-Lpy", GeneName), Species == "Lp")
data_filtered <- data %>% filter(grepl("mitochondrial-Os", GeneName), Species == "Os") %>% rbind(data_filtered)
data_filtered <- data %>% filter(grepl("mitochondrial-Bla", GeneName), grepl("Balano", Species)) %>% rbind(data_filtered)

data_filtered <- data_filtered %>% filter(Cov_max > 50)


## PLOT ##
data2plot <- data_filtered %>%
  pivot_longer(cols = c(A,G,C,T,Del), names_to = "Base", values_to = "Reads") %>%
  mutate(
    MapStatus = case_when(
      Base == Reference ~ "Reference",
      Base == "Del" ~ "Deletion",
      TRUE ~ paste(Base, "Mismatch", sep = "_")
    )
  ) 

data2plot$MapStatus = factor(data2plot$MapStatus, levels = c("Deletion", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Reference"))

data2plot %>%
  ggplot(aes(x=Position, y=Reads, fill=MapStatus)) + 
  geom_col() + 
  facet_wrap(Species~GeneName, ncol=3, scales="free_y")+ 
  theme_bw() + 
  scale_fill_manual(values = c("black", "chartreuse4", "steelblue", "goldenrod", "firebrick", "gray70")) + 
  ylab ("Read Count") + 
  scale_x_continuous(breaks = seq(0,80,5)) + 
  theme(legend.key.size = unit(0.3, "cm"), 
        strip.text = element_text(size=5, margin = margin(2,2,2,2)), 
        axis.text = element_text(size=5), 
        legend.title = element_blank(), 
        legend.text=element_text(size=5.5), 
        axis.title = element_text(size=6,face="bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("Base modifications of mt tRNAs (>50 reads)")
ggsave(paste(dir, "/Mt-tRNAs-base-mod-all.pdf", sep = ""), device = "pdf", width = 9, height = 6)