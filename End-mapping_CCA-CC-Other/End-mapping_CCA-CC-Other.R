#!/usr/bin/env Rscript

## This script aggregates CCA, CC, and other 3′-end categories from multiple
## tRNA CCA-editing summary files and visualizes their relative abundances per library.


## Load required libraries
library(tidyverse)

## Working directory
dir <- getwd()

## List CCA summary files
filenames <- list.files(
  path = dir,
  pattern = "CCA.txt",
  full.names = TRUE
)

## Initialize dataframe
data <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(data) <- c("Ref_tRNA", "CCA_count", "CC_count", "Other_count")

## Load and merge all libraries
for (j in filenames) {
  dataj <- read.csv(j, sep = "\t", header = TRUE)
  dataj$Lib <- j %>%
    str_replace(".*/", "") %>%
    str_replace(".bbmerge.*", "")
  data <- rbind(data, dataj)
}

## Compute totals and metadata
data$Total_count <- data$CCA_count + data$CC_count + data$Other_count

data <- data %>%
  select(Lib, Ref_tRNA, Other_count, CC_count, CCA_count, Total_count)

data$Group <- data$Lib %>% str_replace(".$", "")
data$Compartment <- data$Ref_tRNA %>% str_replace("-.*", "")
data$Species <- data$Lib %>% str_replace("_.*", "") %>% str_replace(".$", "")

## Summarize CCA fractions (Total only)
cca_fraction <- data %>%
  filter(Ref_tRNA == "Total",
         Ref_tRNA != "Bacillus_trnI_control") %>%
  group_by(Lib) %>%
  summarise(
    Count = sum(Total_count),
    CCA_tail = sum(CCA_count, CC_count),
    Other = sum(Other_count)
  ) %>%
  ungroup()

cca_fraction$Species <- cca_fraction$Lib %>%
  str_replace(".$", "") %>%
  str_replace("_.*", "")

## Prepare data for stacked bar plot
fraction <- data %>%
  filter(Ref_tRNA == "Total")

fraction$Species <- fraction$Lib %>%
  str_replace(".$", "") %>%
  str_replace("_.*", "")

fraction <- fraction %>%
  pivot_longer(
    cols = c(CCA_count, CC_count, Other_count),
    names_to = "Category",
    values_to = "Fraction"
  ) %>%
  mutate(Category = factor(
    Category,
    levels = c("Other_count", "CC_count", "CCA_count")
  ))

## Plot CCA / CC / Other fractions
ggplot(
  fraction,
  aes(x = Lib, y = Fraction * 100 / Total_count, fill = Category)
) +
  geom_bar(stat = "identity") +
  labs(x = "Species", y = "Fraction (%)", fill = "Category") +
  scale_y_continuous(labels = seq(0, 100, 25)) +
  scale_fill_brewer(palette = "Greys") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Save figure
ggsave(
  file.path(dir, "cca-cc-other-fractions.pdf"),
  width = 5,
  height = 4
)
