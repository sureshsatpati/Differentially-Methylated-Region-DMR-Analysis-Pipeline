# Differentially-Methylated-Region-DMR-Analysis-Pipeline
This repository contains code and workflows for performing Differentially Methylated Region (DMR) analysis using DMR-seq to identify genomic regions with significant DNA methylation differences across biological conditions.

The pipeline enables robust detection of methylation changes and supports downstream annotation and biological interpretation.

# Key Features:

Input Processing: Supports whole-genome bisulfite sequencing (WGBS), RRBS, or targeted methylation data.

DMR Identification: Uses DMR-seq to detect differentially methylated regions with statistical rigor.

Region-Level Analysis: Aggregates CpG-level methylation signals into biologically meaningful regions.

Annotation of DMRs: Maps DMRs to genomic features such as promoters, enhancers, CpG islands, and gene bodies.

Visualization: Generates methylation profiles, heatmaps, and genome browser–ready tracks.

Comparative Analysis: Identifies condition- or group-specific methylation changes.

# Requirements:

R

dmrseq

bsseq

GenomicRanges

Reference genome and annotation files
