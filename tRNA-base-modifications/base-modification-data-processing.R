#!/usr/bin/env Rscript

library(tidyverse)

dir <- getwd()

## Input files
filenames <- list.files(path = dir, pattern = "parsed.txt$", full.names = TRUE)

print(filenames)

base_mod <- tibble()

for (j in filenames) {
  base_modj <- read.table(
    j,
    sep = "\t",
    header = TRUE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )

  base_modj$Lib <- j %>%
    str_replace(".*/", "") %>%
    str_replace(".bb.*", "")

  base_mod <- bind_rows(base_mod, base_modj)
}



## Metadata parsing
base_mod <- base_mod %>%
  mutate(
    Species = str_replace(Lib, "._$", "") %>% str_replace("_.*", ""),
    Isodecoder = str_replace(name, "[A-Za-z]+-[A-Za-z]+-", "") %>% str_replace("-.*", ""),
    Compartment = str_replace(name, "-.*", ""),
    Replicate = str_replace(Lib, ".*_", ""),
    Total_wDel = Total + Del
  )

## Exclude mitochondrial stemloops
base_mod <- base_mod %>% filter(Isodecoder != "stemloop")

## Sum replicates per gene and position
base_mod_sum <- base_mod %>%
  group_by(Species, name, pos, ref, Isodecoder, Compartment) %>%
  summarise(
    A = sum(A),
    C = sum(C),
    G = sum(G),
    T = sum(T),
    Del = sum(Del),
    Total = A + C + G + T + Del,
    .groups = "drop"
  )

colnames(base_mod_sum) <- c("Species", "GeneName", "Position", "Reference", "Isodecoder", "Compartment", "A", "C", "G", "T", "Del", "Cov_max")

## Export table
write.table(
  base_mod_sum,
  file = file.path(dir, "tRNA-base-modification.txt"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
