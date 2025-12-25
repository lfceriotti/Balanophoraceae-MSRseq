# tRNA base modifications

This analysis quantifies tRNA read coverage and identifies sequence variants associated with base modifications.

For Lophophyum and Ombrophytum nuclear-encoded tRNA-Trps, reads exhibiting specific anticodon modifications were selectively extracted following this [workflow](https://github.com/lfceriotti/Balanophoraceae-MSRseq/edit/main/anticodon-mismatch-SAM-filtering). tRNA read coverage was quantified and sequence variants identified on the outputs. 

Input files containing sequence variant data were generated following the workflow described in [Ceriotti et al. (2024)](https://github.com/dbsloan/Arabidopsis_aminoacylation).

This workflow relies on [`bam-readcount`](https://github.com/genome/bam-readcount) to extract per-base read information from BAM files. The tool was run using the `-t` option to provide a site list file (`site_list.txt`) as input. This file is a tab-delimited text file specifying reference names and coordinate ranges in the format: `ref-name\t1\tref-length`. Details on this usage are discussed in  
[this GitHub thread](https://github.com/genome/bam-readcount/issues/110).

Three R scripts generate plots for:
- mitochondrial tRNAs (mt-tRNAs)
- nuclear-encoded tRNA-Trps in Balanophora
- nuclear-encoded tRNA-Trps in Lophophytum and Ombrophytum

