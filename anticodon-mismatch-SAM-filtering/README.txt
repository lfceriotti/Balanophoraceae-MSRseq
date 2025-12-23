# Anticodon Mismatch-based SAM filtering

This script filters SAM alignments to retain reads that contain exactly one mismatch in the tRNA anticodon, where the observed base matches a user-defined expected modification, and no deletions occur within the anticodon region.

It is designed for analyses of tRNA reads where anticodon-specific base changes (e.g., due to editing or modification) are expected and need to be selectively retained.