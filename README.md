# eMIRNA

eMIRNA is a comprehensive and user-friendly R-based pipeline for predicting and annotating the presence of known and novel microRNAs. This document is intended to give a technical supplementary description about how to run the eMIRNA pipeline through a detailed explanation of all the modules that form part of this program.

## Introduction

The eMIRNA pipeline makes use of a Machine Learning approach based on Support Vector Machine (SVM) algorithm to assess whether putative candidate sequences can be predicted as pre-miRNA-like structures. First, the eMIRNA model must be trained by making use of positive (microRNAs) and negative (other ncRNAs) sequence datasets. Once the SVM model has been trained, putative pre-miRNA candidates can be subjected to prediction.

## Prerequisites

The following R libraries are required for running the eMIRNA pipeline:
+ seqinr (https://CRAN.R-project.org/package=seqinr)
+ stringr (https://CRAN.R-project.org/package=stringr)
+ Biobase (https://bioconductor.org/packages/release/bioc/html/Biobase.html)
+ scales (https://CRAN.R-project.org/package=scales)
+ caret (https://CRAN.R-project.org/package=caret)
+ bimba (https://github.com/RomeroBarata/bimba)
+ LiblineaR (https://CRAN.R-project.org/package=LiblineaR)
+ biomaRt (https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

The following programs are required for running the eMIRNA pipeline:
+ RNAfold [1] (https://www.tbi.univie.ac.at/RNA/)
+ UNAFold [2] (https://github.com/rcallahan/UNAFold)
+ Triplet-SVM [3] (http://www.bioinfo.au.tsinghua.edu.cn/mirnasvm/)
+ BEDTools v2.27.0 [4] (https://bedtools.readthedocs.io/en/latest/)
+ Bowtie [5] (https://mcardle.wisc.edu/mprime/help/bowtie/manual.html)
+ Fasta_ushuffle (https://github.com/agordon/fasta_ushuffle)


All the executables should be stored at computer $PATH in order to be run properly.
