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

## Positive and Negative Datasets

Training the eMIRNA Classifier requires two FASTA files with **Positive** and **Negative** sequences.

The **Positive Sequences** must correspond to those sequences annotated as microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts can be downloaded from the Ensembl repositories available at http://www.ensembl.org/info/data/ftp/index.html.

The **Negative Sequences** must correspond to non-coding sequences other than microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts cand be downloaded from the Ensembl repositories available at http://www.ensembl.org/info/data/ftp/index.html.

In the event that no Reference Genome or no good microRNA or non-coding transcripts are available for downloading, we strongly recommend to choose sequences from the closer phylogenetically related reference species for training the model, otherwise the results can suffer from low reliability.

Both Positive and Negative datasets must be in linear FASTA format. Should you have multilinear FASTA files, they should be converted to linear FASTA files. Users can use the following perl command:

`perl -pe '/^>/ ? print "\n" : chomp' in.fa | tail -n +2 > out.fa`

where `in.fa` corresponds to multilinear FASTA file, and `out.fa` is the resulting linearized FASTA file ready to use.


## eMIRNA.Filter.by.Size

The first eMIRNA module makes use of previous Positive and Negative FASTA files, to perform an initial filtering process based on sequence length. Tipically, microRNA genes range from 50 to 150 nucleotides long. Our first aim would be to filter the selected sequences based on expected microRNA genes length. We will apply this function to both Positive and Negative FASTA files. The Positive sequences should not experience any filtering upon this process, if correctly generated. For Negative sequences, all long non-coding sequences will be removed from our Negative FASTA file, retaining only those sequences resembling microRNA genes in length, according to estabished thresholds.

This function requires four arguments:

+ PATH to Positive or Negative FASTA file.
+ String with desired output prefix name.
+ Lower length filtering threshold.
+ Upper length filtering threshold.

We recomend to set 50 nucleotides for lower length threshold, and 150 for the upper, but users can define their own limit thresholds.

Example of usage:

`eMIRNA.Filter.by.Size("PATH to Positive FASTA file", "Pos", 50, 100)`
`eMIRNA.Filter.by.Size("PATH to Negative FASTA file", "Neg", 50, 100)`

Once the eMIRNA.Filter.by.Size function has run, a new folder named eMIRNA will have been created at your computer $HOME, with a subfolder, inside eMIRNA/ folder, called FilterSize_Results/, in which a FASTA file named Pos/Neg_filter_size.fa will be generated with the results of running the function.


eMIRNA.Filter.by.Structure
The second eMIRNA module aims to estimate the secondary folding structure of selected filtered sequences both in Positive and Negative datasets, thus filtering out all sequences that do not resemble a pre-miRNA hairpin-like secondary structure. The eMIRNA.Filter.by.Structure function will make use of RNAfold [1] program to calculate the estimated secondary folding structure, which should be available in your computer $PATH to be correctly executed. Tipically, microRNA genes have a characteristic secondary structure, composed by two stems joined by complementarity and one terminal loop, forming a hairpin-like secondary structure. Some bubbles or bulges can appear within the two stems, belonging to non-paired nucleotides in the sequence.

This function requires two arguments:
PATH to eMIRNA.Filter.by.Size Positive or Negative FASTA output files.
String with desired output prefix name.

Example of usage:
> eMIRNA.Filter.by.Structure("~/eMIRNA/FilterSize_Results/Pos_filter_size.fa", "Pos")
> eMIRNA.Filter.by.Structure("~/eMIRNA/FilterSize_Results/Neg_filter_size.fa", "Neg")

Once the eMIRNA.Filter.by.Structure has run, a new folder named FilterSstructure_Results/ will be created inside eMIRNA/ folder, in which a FASTA file called Pos/Neg_filter_nloop.fa will be generated with the results of running the function.


eMIRNA.Features
The third eMIRNA module aims to calculate a series of structural, statistical and sequence-derived features from each sequence that had passed previous filterings, in order to obtain an estimated representation of their structural characteristics. Afterwards, these feature matrix will be processed by the prediction software to discriminate between microRNAs and other type of sequences.

The function requires two arguments:
PATH to eMIRNA.Filter.by.Structure Positive or Negative FASTA output files.
String with desired output prefix name.

Example of usage:
> Pos = eMIRNA.Features("~/eMIRNA/FilterStructure_Results/Pos_filter_nloop.fa”, "Pos")
> Neg = eMIRNA.Features("~/eMIRNA/FilterStructure_Results/Neg_filter_nloop.fa”, "Neg")

Once the eMIRNA.Features has run, a new folder named Features_Results/ will be created inside eMIRNA/ folder, in which a .csv file called Pos/Neg.csv will be generated with the results of running the function.

The eMIRNA.Features function includes the calculation of a total of 6 Sequence Features, comprising 55 variables, 8 Secondary Structure Features, comprising 23 variables, and 24 Structural Statistics. 
Sequence Features:
32 Triplet elements calculated with SVM-Triplets pipeline [3] (T1 - T32).
Sequence Length (Length).
Guanine+Cytosine/Length (GC).
Adenine+Uracil / Guanine+Cytosine ratio (AU.GCr).
Adenine, Uracil, Guanine, Cytosine / Length (Ar, Ur, Gr, Cr).
Dinucleotide content / Length (AAr, GGr, CCr, UUr, AGr, ACr, AUr, GAr, GCr, GUr, CAr, CGr, CUr, UAr, UGr, UCr).
Secondary Structure Features:
Terminal hairpin loop length (Hl).
5’ - 3’ Stems length (Steml5, Steml3).
Base pairs in Secondary Structure (BP).
Base pairs or matches at 5’ - 3’ Stems (BP5, BP3).
Unpaired bases or mismatches at 5’ - 3’ Stems (Mism5, Mism3).
Bulges at 5’ - 3’ Stems (B5, B3).
Bulges at 5’ - 3’ Stems of type 1,2,3,4 or 5 unpaired bases (BN1.5, BN1.3 … BN5.5, BN5.3).
A-U, G-C and G-U Pairings in sequence (AUp, GCp, GUp).
Structural Statistics:
Minimum Free Energy estimated by RNAfold [1] (MFE).
Ensemble Free Energy (EFE).
Centroid Free Energy (CFE).
Centroid Distance to Ensemble (CDE).
Maximum Expected Accuracy of Free Energy (MEAFE).
Maximum Expected Accuracy (MEA).
BP / Length (BPP).
MFE Ensemble Frequency (EFreq).
Ensemble Diversity (ED).
MFE / Length (MFEadj).
EFE / Length (EFEadj).
CDE / Length (Dadj).
Shannon Entropy / Length (SEadj).
MFE – EFE / Length (DiffMFE.EFE).
MFEadj / GC (MFEadj.GC).
MFEadj / BP (MFEadj.BP).
MFE estimated by UNAFold [2] (dG).
dG / Length (dGadj).
Structural Entropy estimated by Melt function from UNAFold [2](dS).
dS / Length (dSadj).
Structural Enthalpy estimated by Melt function (dH).
dH / Length (dHadj).
Fusion Temperature estimated by Melt function (dT).
dT / Length (dTadj).


eMIRNA.Train
The fourth eMIRNA module aims to perform the training process of a Machine Learning-based Support Vector Machine (SVM) algorithm, making use of Feature representation previously calculated, to construct a SVM model capable to distinguish between microRNAs and other non-coding sequences. The SVM algorithm is built by using a 10 k-fold cross-validation over a training set randomly selected from analyzed features.

This function requires three arguments:
Positive Features calculated by eMIRNA.Features, saved in R object.
Negative Features calculated by eMIRNA.Features, saved in R object.
Imbalance correction algorithm.


Example of usage:
> SVM = eMIRNA.Train(Pos, Neg, imbalance=”smote”)
It is important that when running the SVM training process, both Positive and Negative matrices have a balanced number of sequences to evaluate, keeping the number of positive and negative sequences to be similar, in order to avoid an overrepresentation of one of the two classes. To overcome this problem, eMIRNA.Trains implements a series of imbalance correction methods available. If required, eMIRNA.Train will first perform a Noise Reduction A Priori Synthetic correction (NRAS) of input features, as reported by Rivera W [6], followed by the preferred method to over-sampling the minority class to correct class-imbalance biases. Available methods are (adasyn, bdlsmote1, bdlsmote2, mwmote, ros, rwo, slsmote, smote):
ADASYN: Adaptive Synthetic Sampling [7]
BDLSMOTE: borderline-SMOTE1 and borderline-SMOTE2 [8]
MWMOTE: Majority Weighted Minority Over-Sampling TEchnique [9]
ROS: Random Over-Sampling
RWO: Random Walk Over-Sampling [10]
SLSMOTE: Safe-Level-SMOTE [11]
SMOTE: Synthetic Minority Over-Sampling TEchnique [12]

By default, eMIRNA.Train will not perform any class-imbalance correction, but users are imperiously advised to do so. Otherwise the training process could suffer.
Once the function has run, eMIRNA.Train will create a SVM classifier capable to differentiate between microRNAs and other structurally microRNA-like non-coding RNAs. 
