# Anticodon Mismatch-based SAM filtering

This script filters SAM alignments to retain reads that contain exactly one mismatch in the tRNA anticodon, where the observed base matches a user-defined expected modification, and no deletions occur within the anticodon region.

It is designed for analyses of tRNA reads where anticodon-specific base changes (indicated in "mismatch-file.txt") are expected and need to be selectively retained. It was applied primarily to the analysis of nuclear-encoded Trp tRNA genes.

Annotation file consists of [Sprinzl coordinates](https://github.com/dbsloan/Arabidopsis_aminoacylation/tree/main/Sprinzl_coordinates), and includes an additional column ("Note") which is used to denote anticodon positions within the anticodon loop.

Run the script using:

`python anticodon-mismatch-SAM-filtering.py annotation-file.txt mismatch-file.txt input.sam output.sam`
