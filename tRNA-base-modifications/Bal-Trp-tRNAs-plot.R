#!/usr/bin/env Rscript

## Plot base-modification profiles of nuclear Trp tRNAs in Balanophoraceae

library(tidyverse)

## Read input filename from command line
args <- commandArgs(trailingOnly = TRUE)
file <- as.character(args[1])

## Working directory
dir <- getwd()

## Load aggregated tRNA base-modification data
data <- read.csv(
  paste0(dir, "/", file),
  sep = "\t",
  header = TRUE
)

## Select nuclear Trp tRNAs from Balanophoraceae
data_filtered <- data %>%
  filter(grepl("nuclear-Ath-Trp", GeneName), grepl("Balano", Species))

data_filtered <- data %>%
  filter(grepl("nuclear-Bya", GeneName), grepl("Balano", Species)) %>%
  bind_rows(data_filtered)

## Retain positions with sufficient coverage
data_filtered <- data_filtered %>% filter(Cov_max > 50)

## Reshape data and classify mapping outcomes
data2plot <- data_filtered %>%
  pivot_longer(
    cols = c(A, G, C, T, Del),
    names_to = "Base",
    values_to = "Reads"
  ) %>%
  mutate(
    MapStatus = case_when(
      Base == Reference ~ "Reference",
      Base == "Del"     ~ "Deletion",
      TRUE              ~ paste(Base, "Mismatch", sep = "_")
    )
  )

## Define plotting order for base categories
data2plot$MapStatus <- factor(
  data2plot$MapStatus,
  levels = c(
    "Deletion",
    "A_Mismatch",
    "C_Mismatch",
    "G_Mismatch",
    "T_Mismatch",
    "Reference"
  )
)

## Plot per-position base composition
data2plot %>%
  ggplot(aes(x = Position, y = Reads, fill = MapStatus)) +
  geom_col() +
  facet_wrap(Species ~ GeneName, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(
    values = c("black", "chartreuse4", "steelblue",
               "goldenrod", "firebrick", "gray70")
  ) +
  ylab("Read Count") +
  scale_x_continuous(breaks = seq(0, 80, 5)) +
  theme(
    legend.key.size = unit(0.3, "cm"),
    strip.text = element_text(size = 6, margin = margin(2, 2, 2, 2)),
    axis.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.5),
    axis.title = element_text(size = 8, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Base modifications of nuclear Trp tRNAs")

## Save figure
ggsave(
  paste0(dir, "/Nuclear-Bal-Trp-base-mod.pdf"),
  device = "pdf",
  width = 6,
  height = 3
)
