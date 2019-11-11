# eMIRNA

eMIRNA is a comprehensive and user-friendly pipeline for predicting and annotating the presence of known and novel microRNAs. This document is intended to give a technical supplementary description about how to run the eMIRNA pipeline through a detailed explanation of all the modules that form part of this program.

The eMIRNA pipeline is under active development, if you find any problem, doubt or bug while running it, please contact emilio.marmol@cragenomica.es. eMIRNA pipeline for plants will be released soon.

&nbsp;
&nbsp;

-------------------------------

## Index

- [Introduction](https://github.com/emarmolsanchez/eMIRNA/#introduction)
- [Prerrequisites](https://github.com/emarmolsanchez/eMIRNA/#prerequisites)
- [Positive, Negative and Unlabeled Data sets](https://github.com/emarmolsanchez/eMIRNA/#positive-negative-and-unlabeled-data-sets)
- [eMIRNA.Filter](https://github.com/emarmolsanchez/eMIRNA/#emirnafilter)
- [eMIRNA.Features](https://github.com/emarmolsanchez/eMIRNA/#emirnafeatures)
- [eMIRNA.Train](https://github.com/emarmolsanchez/eMIRNA/#emirnatrain)
- [eMIRNA.Hunter](https://github.com/emarmolsanchez/eMIRNA/#emirnahunter)
- [eMIRNA.Predict](https://github.com/emarmolsanchez/eMIRNA/#emirnapredict)
- [eMIRNA.Refiner](https://github.com/emarmolsanchez/eMIRNA/#emirnarefiner)
- [eMIRNA.Structural.Pvalue](https://github.com/emarmolsanchez/eMIRNA/#emirnastructuralpvalue)
- [References](https://github.com/emarmolsanchez/eMIRNA/#references)
- [Contact](https://github.com/emarmolsanchez/eMIRNA/#contact)
- [Notes](https://github.com/emarmolsanchez/eMIRNA/#notes)

&nbsp;
&nbsp;

## Introduction

The eMIRNA pipeline makes use of a Machine Learning approach based on semi-supervised transductive Graph-based algorithm [[1]] to assess whether putative candidate sequences can be predicted as pre-miRNA-like structures.

&nbsp;

![alt text](https://github.com/emarmolsanchez/eMIRNA/blob/master/bin/Figure1.jpg)

**(1)** Positive, negative and unlabeled data are filtered based on size and secondary folding structure and a set of features is extracted for each sequence. **(2)** Mature miRNA sequences from small RNA-Seq data or related reference species are mapped against the selected genome assembly and elongated to reconstruct putative pre-miRNA candidates. **(3)** Candidate precursors are filtered based on size and secondary folding structure and a set of features is extracted for each candidate sequence. Optionally, sequences showing unstable secondary structure are removed. **(4)** Candidate sequences are embedded in the semi-supervised transductive classifier and a list of putative miRNAs is predicted. **(5)** Predicted miRNAs are either assigned to already annotated miRNA loci in the selected reference assembly or classified as putative novel miRNA genes.



&nbsp;

## Prerequisites

The following R libraries are required for running the eMIRNA pipeline:
+ [stringr]
+ [seqinr]
+ [Biobase]
+ [scales]
+ [miRNAss] [[1]]
+ [PCIT] [[2]]
+ [NOISeq] [[3]]
+ [edgeR] [[4]]
+ [igraph] [[5]]

The following software programs are required for running the eMIRNA pipeline:
+ [RNAfold] [[6]]
+ [BEDTools v2.27.0] [[7]]
+ [Bowtie] [[8]]
+ [Fasta_ushuffle]


All executables should be stored at computer `$PATH` in order to be run properly (Commonly located at `/usr/bin/` or `/usr/local/bin/` in UNIX systems).

&nbsp;

## Positive, Negative and Unlabeled Data sets

Running the **eMIRNA** Classifier requires two FASTA files with **Positive** and **Negative** sequences.

The **Positive Sequences** must correspond to those sequences annotated as microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts can be downloaded from the [Ensembl repositories].

The **Negative Sequences** must correspond to non-coding sequences other than microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts can be downloaded from the [Ensembl repositories].

Additionally, other hairpin-like sequences can be extracted from the reference Genome, in order to increase the variety and diversity of sequences to be included during graph reconstruction and allow a better topological adjustment of positive and negative categories. The [HextractoR] package can be used for generating a set of hairpin-like sequences from any available genome assembly. Please be aware that quering the whole genome can be extremely time consuming and resource intensive. We reccomend establishing random blocks (1-10 Mb) within each genomic chromosome. As no prior knowledge of the identity of randomly extracted hairpins would available, they will be set as **Unlabeled sequences**.

Once positive, negative and unlabeled sequences are available, we recomend using a identity-by-sequence filtering step, in order to remove redundant sequences that may occur in any of the categories. The [CD-Hit] suite is a good resource for implementing sequence-based repetitve elements removal.

In the event that no Reference Genome or no good microRNA or non-coding transcripts are available for downloading, we strongly recommend choosing sequences from the closest phylogenetically related reference species with available genome annotation for running the classification procedure, otherwise the results can suffer from low reliability.

Positive, Negative and Unlabeled datasets must be in linear FASTA format. Should you have multilinear FASTA files, they must be converted to linear FASTA files. Users can implement the following perl command:

`perl -pe '/^>/ ? print "\n" : chomp' in.fa | tail -n +2 > out.fa`

where `in.fa` corresponds to multilinear FASTA file, and `out.fa` is the resulting linearized FASTA file ready to use.

&nbsp;

## eMIRNA.Filter

The first eMIRNA module makes use of previous Positive, Negative and Unlabeled FASTA files to perform an initial filtering process based on sequence length. Typically, microRNA genes range from 50 to 150 nucleotides long. Our first aim would be to filter the selected sequences based on expected microRNA genes length. We will apply this function to each of the FASTA files. The Positive sequences should not experience any filtering upon this process, if correctly generated. For Negative and Unlabeled sequences, all long non-coding hairpin-like sequences will be removed, retaining only those sequences resembling microRNA genes in length, according to established thresholds. 

Next, this module estimates the secondary folding structure of selected filtered sequences, thus filtering out all candidates that do not ressemble a pre-miRNA hairpin-like secondary structure. The eMIRNA.Filter function will make use of [RNAfold] software [[6]] to calculate the estimated secondary folding structure, which should be available in your computer `$PATH` to be correctly executed. Typically, microRNA genes have a characteristic secondary structure, composed by two stems joined by complementarity and one terminal loop, forming a hairpin-like secondary structure. Some bubbles or bulges can appear within the two stems, belonging to non-paired nucleotides in the sequence.

This function requires four arguments:

+ PATH to Positive or Negative FASTA file.
+ String with desired output prefix name.
+ Lower length filtering threshold.
+ Upper length filtering threshold.

We recommend setting 50 nucleotides as lower length threshold, and 150 for the upper, but users can define their own preferred thresholds.

Example of usage:

```r

eMIRNA.Filter.by.Size("PATH to Positive FASTA file", "Pos", 50, 150)

eMIRNA.Filter.by.Size("PATH to Negative FASTA file", "Neg", 50, 150)

eMIRNA.Filter.by.Size("PATH to Unlabeled FASTA file", "Unlab", 50, 150)

```

Once the eMIRNA.Filter function has run, a new folder named `eMIRNA/` will be created at your computer `$HOME`, with a subfolder called `Filter_Results/` inside, in which a FASTA file named `Pos/Neg/Unlab_filtered.fa` will be generated with the results of running the function.

&nbsp;

## eMIRNA.Features

The third eMIRNA module aims to calculate a series of structural, statistical and sequence-derived features from each sequence that had passed previous filtering, in order to obtain an estimated representation of their structural characteristics. Subsequently, these feature matrices will be processed by the prediction software to discriminate between microRNAs and other type of sequences.

A modified version of Triplet-SVM pipeline [[9]] is implemented in the eMIRNA.Features module. Triplet-SVM perl scripts 1 to 3 (available at `bin/`) should be located at computer `$PATH`, so that the function is properly executed. [RNAfold] executable must also be installed and available at `$PATH`. All Triplet-SVM perl scripts should have execution permission allowed, which can easily be set with command `chmod 777`. 

The function requires two arguments:

+ PATH to eMIRNA.Filter FASTA output files.
+ String with desired output prefix name.

Example of usage:

```r

Pos <- eMIRNA.Features("~/eMIRNA/Filter_Results/Pos_filtered.fa", "Pos")

Neg <- eMIRNA.Features("~/eMIRNA/Filter_Results/Neg_filtered.fa", "Neg")

Unlab <- eMIRNA.Features("~/eMIRNA/Filter_Results/Unlab_filtered.fa", "Unlab")

```

Once the `eMIRNA.Features` has run, a new folder named `Features_Results/` will be created inside `eMIRNA/`, in which a .txt file called `Pos/Neg/Unlab.txt` will be generated with the results of running the function.

The eMIRNA.Features function includes the calculation of a total of 6 Sequence Features, comprising 55 variables, 8 Secondary Structure Features, comprising 27 variables, and 17 Structural Statistics. 

Sequence Features:

+ 32 Triplet elements calculated with SVM-Triplets pipeline (T1 - T32).
+ Sequence Length (Length).
+ Guanine+Cytosine/Length (GC).
+ Adenine+Uracil / Guanine+Cytosine ratio (AU.GCr).
+ Adenine, Uracil, Guanine, Cytosine / Length (Ar, Ur, Gr, Cr).
+ Dinucleotide content / Length (AAr, GGr, CCr, UUr, AGr, ACr, AUr, GAr, GCr, GUr, CAr, CGr, CUr, UAr, UGr, UCr).

Secondary Structure Features:

+ Terminal hairpin loop length (Hl).
+ 5’ - 3’ Stems length (Steml5, Steml3).
+ Base pairs in Secondary Structure (BP).
+ Base pairs or matches at 5’ - 3’ Stems (BP5, BP3).
+ Unpaired bases or mismatches at 5’ - 3’ Stems (Mism5, Mism3).
+ Bulges at 5’ - 3’ Stems (Bulge5, Bulge3).
+ Bulges at 5’ - 3’ Stems of type 1, 2, 3, 4, 5, 6 or 7 unpaired bases (BN1.5, BN1.3 … BN7.5, BN7.3).
+ A-U, G-C and G-U Pairings in sequence (AUp, GCp, GUp).

Structural Statistics:

+ Minimum Free Energy estimated by RNAfold (MFE).
+ Ensemble Free Energy (EFE).
+ Centroid Free Energy (CFE).
+ Centroid Distance to Ensemble (CDE).
+ Maximum Expected Accuracy of Free Energy (MEAFE).
+ Maximum Expected Accuracy (MEA).
+ BP / Length (BPP).
+ Ensemble Diversity (ED).
+ MFE / Length (MFEadj).
+ EFE / Length (EFEadj).
+ CDE / Length (Dadj).
+ MEAFE / Length (MEAFEadj).
+ ED / Length (EDadj).
+ Shannon Entropy / Length (SEadj).
+ MFE – EFE / Length (DiffMFE.EFE).
+ MFEadj / GC (MFEadj.GC).
+ MFEadj / BP (MFEadj.BP).

&nbsp;

## eMIRNA.Hunter
The eMIRNA.Hunter module is an auxiliar BASH script developed to obtain pre-miRNA candidate sequences, For doing so, users can choose either an homology-based recovery from previously annotated microRNAs in reference species, or the use of reads derived from small RNA-Seq experiments.

The eMIRNA.Hunter script implements Bowtie (Langmead et al., 2009) for the alignment of putative mature microRNA sequences against the reference assembly selected by the user, reconstructing pre-miRNA sequence candidates and generating a FASTA and BED files for the candidates to be subsequently classified.

Users should provide a properly collapsed FASTA file with small RNA-Seq sequences from canonical FASTQ sequence files, or a list of unique annotated mature miRNA sequences in FASTA format. The FASTQ files should be quality-check filtered and sequencing adaptor trimmed before running any available collapser tool, e.g. FASTQ collapser from [FASTX-Toolkit] for collapsing FASTQ files into FASTA files with uniquely represented sequences.

Once the module has run, users can optionally perform an additional filter on generated FASTA pre-miRNA candidates according to their secondary folding structural stability by using the eMIRNA.Structural.Pvalue module.

Users are encouraged to perform a pre-filtering process of the collapsed FASTA file to retain sequences between 18-25 nucleotides in length, corresponding to the average length of mature miRNAs. The eMIRNA.Filter.by.Size module could be used for this purpose.

This module requires six arguments:

+ PATH to Reference Genome Bowtie Index from your species of interest.
+ PATH to FASTA file of selected reference organism mature microRNAs.
+ PATH to desired output.
+ Desired output prefix name.
+ Upwards number of bases for pre-miRNA reconstruction
+ Backwards number of bases for pre-miRNA reconstruction

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Hunter Usage Instructions:
eMIRNA.Hunter [options]
Input:
 -r                                     PATH to Species of interest Genome Bowtie Index
 -f                                     PATH to Reference FASTA file
 -o                                     PATH to desired output folder
 -x                                     Desired Name string for output files
 -u                                     Upwards number of bases for pre-miRNA reconstruction (30-80 bp recommended, 70-80 for CNNC motif correction)
 -b                                     Backwards number of bases for pre-miRNA reconstruction (1-10 bp recommended, 20-30 for UG motif correction)
 -h                                     Display help page
Output:
 <x>.sam                                SAM file output from Bowtie miRNA alignment
 <x>.log                                LOG file output from Bowtie miRNA alignment
 <x>_homolog_miRNAs.bed                 BED file output with candidate pre-miRNAs for prediction
 <x>_homolog_miRNAs.fa                  FASTA file output with candidate pre-miRNAs for prediction
 <x>_homolog_miRNAs_corrected.bed       BED file output with motif corrected candidate pre-miRNAs for prediction
 <x>_homolog_miRNAs_corrected.fa        FASTA file output with motif corrected candidate pre-miRNAs for prediction
 
 ```

For achieving a successful cross-species alignment, it is very important that mature microRNA sequences FASTA file from reference organism are in DNA code, with Ts for Thymine and no Us for Uracil in RNA code. Please make sure that your mature microRNA sequences are in DNA code, otherwise the alignment process will fail.

For generating Bowtie Index for your Reference Genome, please refer to Bowtie Manual.

Example of usage:

```
bash eMIRNA.Hunter -r PATH_to_Bowtie_Index -f Hsapiens_mature_miRNAs.fa -o PATH_to_output_folder -x Candidates -u 70 -b 30

```

&nbsp;

After successfully running eMIRNA.Hunter script, six files will have been created at predefined output PATH:

+ SAM file with aligned sequences.
+ Log file with alignment statistics.
+ BED file with positions of reconstructed premiRNA homologous candidates.
+ FASTA file with reconstructed premiRNA homologous candidates.
+ BED file with motif corrected positions of reconstructed pre-miRNA homologous candidates.
+ FASTA file with motif corrected reconstructed pre-miRNA homologous candidates.

Once the FASTA file with pre-miRNA candidates has been generated, users must process these sequences by following the previously described steps for eMIRNA pipeline, in order to obtain a Feature matrix representing those candidate sequences that will then be subjected to classification by the SVM trained algorithm.

Optionally, users can subject the motif corrected FASTA file for further prediction, taking into consideration that not all miRNAs would be processed following motif detection and thus some novel candidates may be missed. On the contrary, a much more accurate positioning for pre-miRNA candidates where delimiting motifs had been successfully detected wil be calculated, as reported by Auyeung *et al.* (Auyeung et al., 2013).

&nbsp;

## eMIRNA.Hunter_denovo

The eMIRNA.Hunter_denovo module is a modified versión of eMIRNA.Hunter module, developed to obtain pre-miRNA candidate sequences to perform a *de novo* discovery of novel putative miRNAs from smallRNA-seq data.

Users should provide a properly collapsed FASTA file with smallRNA-seq sequences from canonical FASTQ sequence files. The FASTQ files should be quality-check filtered and sequencing adaptor trimmed before running any available collapser tool, e.g. FASTQ collapser from FASTX-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/index.html) for collapsing FASTQ files into FASTA files with uniquely represented sequences.

Users are encouraged to perform a pre-filtering process of the collapsed FASTA file to retain sequences between 18-25 nucleotides in length, corresponding to the average length of mature miRNAs. The eMIRNA.Filter.by.Size module could be used for this purpose.

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Hunter_denovo Usage Instructions:
eMIRNA.Hunter_denovo [options]
Input:
  -r                               PATH to Species of interest Genome Bowtie Index
  -f                               PATH to Collapsed smallRNA-seq FASTA file
  -o                               PATH to desired output folder
  -x                               Desired Name string for output files
  -u                               Upwards number of bases for pre-miRNA reconstruction (30-80 bp recommended, 70-80 for CNNC motif correction)
  -b                               Backwards number of bases for pre-miRNA reconstruction (1-10 bp recommended, 20-30 for UG motif correction)
  -h                               Display help page
Output:
  <x>.sam                          SAM file output from Bowtie miRNA alignment
  <x>.log                          LOG file output from Bowtie miRNA alignment
  <x>_miRNAs.bed                   BED file output with candidate pre-miRNAs for prediction
  <x>_miRNAs.fa                    FASTA file output with candidate pre-miRNAs for prediction
  <x>_miRNAs_corrected.bed         BED file output with motif corrected candidate pre-miRNAs for prediction
  <x>_miRNAs_corrected.fa          FASTA file output with motif corrected candidate pre-miRNAs for prediction
  
  ```
Example of usage:

```

bash eMIRNA.Hunter_denovo -r PATH_to_Bowtie_Index -f Small-RNAseq_collapsed_fastq.fa -o PATH_to_output_folder -x Candidates -u 70 -b 30

```

&nbsp;
 
## eMIRNA.Predict

The fifth eMIRNA module aims to perform microRNA classification by making use of the previously trained SVM algorithm and candidates sequences generated by eMIRNA.Hunter. Users must pre-process the candidate sequences FASTA files generated by eMIRNA.Hunter module with eMIRNA.Filter.by.Size and eMIRNA.Filter.by.Structure modules, followed by the calculation of representing features matrix by eMIRNA.Features. Either the positionally corrected or not corrected FASTA files output from eMIRNA.Hunter can be selected for this task, but users must use the corresponding corrected/not corrected BED files for further steps in the eMIRNA pipeline.

This module requires three arguments:

+ SVM trained algorithm object.
+ Feature matrix representing candidate sequences to evaluate.
+ String with desired output prefix name.

Example of usage:

```r

eMIRNA.Predict(SVM, Candidates_Feature_Matirx, "Candidates")

```

Once the eMIRNA.Predict has run, a new folder named `Prediction_Results/` will be created inside `eMIRNA/` folder, in which a .txt file called `Candidates.txt` will be generated, containing a list of Sequence candidate names classified as putative miRNAs by the SVM trained algorithm.

&nbsp;

## eMIRNA.Refiner

After generating a list of putative new microRNAs by implementing the eMIRNA pipeline, users should filter these new candidates in order to select those more likely to be new putative microRNA candidates in the species of interest. The eMIRNA.Refiner module has been specifically designed to cover this task. This BASH script complements eMIRNA.Hunter and eMIRNA.Predict modules by selecting those candidate sequences that have more probability to correspond to new non-annotated microRNA sequences orthologous to the reference organism annotated microRNAs in which candidate sequences reconstruction has been performed.

The eMIRNA.Refiner module makes use of an auxiliar R script called eMIRNA_biomaRt_calc.R, which should be located at the same `$PATH` where eMIRNA.Refiner is found. The functionality of eMIRNA.Refiner is based on comparing the microRNA gene neighborhood both in the reference species and in the species under study. Many shared orthologous microRNA sequences will not only share a conserved sequence between species, but also the genome location in which they are placed, making them more reliable if sharing both sequence and location. The candidate microRNA genes neighborhood is estimated by establishing a comparison window around the sequence (typically 2-10 Mb) and contrasting orthologous genes located around the corresponding microRNA in the reference organism. By doing so, the eMIRNA.Refiner module can calculate a Neighborhood Score, measuring the ratio of similarity among annotated genes in the locus of interest between the reference organism and our species of interest. The higher this ratio would be, the most reliable the new microRNA candidate can be considered for further analyses, meaning that not only the sequence is conserved between species, but also its genetic topology, highlighting the putative metabolic importance of this miRNA gene.

This module requires nine arguments:

+ PATH to GTF annotation file from the species of interest.
+ PATH to GFF miRNA annotation file from selected reference species.
+ PATH to list of candidate sequences generated by eMIRNA.Predict module.
+ PATH to BED file with positions of reconstructed pre-miRNA orthologous candidates by eMIRNA.Hunter.
+ Window size in Mb for Neighbouring Genes between-species contrast.
+ Name of the species of interest.
+ Name of the reference species to contrast.
+ PATH to desired output.
+ Desired output prefix name.

Users should check the genome availability at biomaRt repositories by using the the following R code:

`mart = useMart('ensembl')`, followed by `listDatasets(mart)`

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Refiner Usage Instructions:
eMIRNA.Refiner [options]
Input:
 -g	                                      PATH to Species of interest Genome GTF Annotation file
 -i	                                      PATH to Reference Organism microRNA GFF Annotation file
 -p	                                      PATH to List of Predicted microRNA candidates by eMIRNA
 -b	                                      PATH to BED Ortholog microRNAs output file from eMIRNA.Hunter
 -d	                                      Mb window for between-species Neighbouring Genes Contrast (2-5 recommended)
 -s	                                      Name of species of interest (scientific name in lower case, e.g. sscrofa)
 -m	                                      Name of reference species for contrast (scientific name in lower case, e.g. hsapiens)
 -o	                                      PATH to desired output folder
 -x	                                      Desired Name string for output files
 -h	                                      Display help page
Output:
<x>_Predicted_miRNAs_annotated.bed            BED file output with eMIRNA Predicted already Annotated miRNAs
<x>_Predicted_miRNAs_NON_annotated.bed        BED file output with eMIRNA Predicted Novel miRNAs
<x>_Putative_Predicted_miRNAs.txt             TXT file output with Putative miRNA Candidates filtered by Neighbouring Score

```
Example of usage:

```
bash eMIRNA.Refiner -g Sscrofa11.1.94.gtf -i Hsapiens_miRNA_mirbase.gff -p Candidates.txt -b Candidates_homolog_miRNAs_corrected.bed -d 2 -s sscrofa -m hsapiens -o PATH_to_output_folder -x Candidates

```

After successfully running the eMIRNA.Refiner script, a .txt file will have been created at predefined output `PATH`, containing the most relevant putative novel non-annotated microRNAs, their estimated positions and corresponding Neighbouring Score. Besides, two BED files, containing novel and annotated candidates will be generated.

&nbsp;

## eMIRNA.Refiner_denovo

The eMIRNA.Refiner_denovo module is a modified versión of eMIRNA.Refiner module, developed to obtain pre-miRNA candidate sequences to perform a *de novo* discovery of novel putative miRNAs from smallRNA-seq data.

This module requires six arguments:

+ PATH to GTF annotation file from the species of interest.
+ PATH to list of candidate sequences generated by eMIRNA.Predict module.
+ PATH to BED file with positions of reconstructed pre-miRNA candidates by eMIRNA.Hunter_denovo.
+ PATH to FASTA file with sequences of reconstructed pre-miRNA candidates by eMIRNA.Hunter_denovo.
+ PATH to desired output.
+ Desired output prefix name.

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Refiner_denovo Usage Instructions:
eMIRNA.Refiner_denovo [options]
Input:
  -g                                             PATH to Species of interest Genome GTF Annotation file
  -p                                             PATH to List of Predicted microRNA candidates by eMIRNA
  -b                                             PATH to BED microRNAs output file from eMIRNA.Hunter_denovo
  -f                                             PATH to FASTA microRNAs output file from eMIRNA.Hunter_denovo
  -o                                             PATH to desired output folder
  -x                                             Desired Name string for output files
  -h                                             Display help page
Output:
  <x>_Predicted_miRNAs_annotated.bed             BED file output with eMIRNA Predicted already Annotated miRNAs
  <x>_Predicted_miRNAs_NON_annotated.bed         BED file output with eMIRNA Predicted Novel miRNAs
  <x>_Predicted_miRNAs_NON_annotated.fa          FASTA file output with eMIRNA Predicted Novel miRNAs
  
  ```

Example of usage:

```
bash eMIRNA.Refiner_denovo -g Sscrofa11.1.94.gtf -p Candidates.txt -b Candidates_miRNAs_corrected.bed -f Candidates_miRNAs_corrected.fa -o PATH_to_output_folder -x Candidates

```

After successfully running the eMIRNA.Refiner_denovo script, a BED file will have been created at predefined output `PATH`, containing the most relevant putative novel non-annotated microRNAs, and their estimated positions. Besides, a BED file containing already annotated detected candidates and FASTA file with correspondent candidate sequences will be generated.

&nbsp;

## eMIRNA.Structural.Pvalues

Finally, after having obtained a list of putative novel pre-miRNA sequences by the aforementioned eMIRNA functions, users can analyze if the structural integrity of predicted pre-miRNAs can achieve a stable conformation at a statistically significant level. The eMIRNA.Structural.Pvalues function implements a n-randomization of provided sequences while maintaining k-let counts as described by Jiang *et al*. (Jiang et al., 2008), using the fasta_ushuffle available at https://github.com/agordon/fasta_ushuffle. The fasta_ushuffle software must be downloaded, compiled and stored at computer `PATH` in order to be run properly.

This module requires three arguments:

+ PATH to FASTA file of putative novel miRNA candidates generated by eMIRNA.Refiner_denovo.
+ String with desired output prefix name.
+ Number of iterations to perform.

Example of usage:

```r

eMIRNA.Structural.Pvalues("~/eMIRNA/Candidates_Predicted_miRNAs_NON_annotated.fa", "Candidates", 100)

```

By default, eMIRNA.Structural.Pvalues will perform 100 random shuffling iterations over each provided sequence. Users can set their desired number of iterations but should be aware of computing times required for iterating and folding of secondary structures for each sequence. As computing costs can exponentially increase with higher number of iterations, we encourage users to set their desired range of iterations between 100 and 1000, depending on the number of candidate sequences to be analyzed.

Once the eMIRNA.Structural.Pvalues has run, a new .csv file called `Candidates_Structural_Pvalues.csv` will be generated at `/Prediction_Results` folder, containing calculated MFE and EFE Structural *P*-values for each candidates sequence. Users can then select those sequences with a significantly stable secondary structure folding as sequences having higher probability to have been properly profiled and predicted (i.e. *P*-values < 0.05).

&nbsp;

## References

Auyeung V.C. et al. (2013) Beyond secondary structure: primary-sequence determinants license pri-miRNA hairpins for processing. *Cell*, 152, 844–58.

Barua S. et al. (2014) MWMOTE--Majority Weighted Minority Oversampling Technique for Imbalanced Data Set Learning. *IEEE Trans. Knowl. Data Eng.*, 26, 405–25.

Bunkhumpornpat C. et al. (2013) Safe-Level-SMOTE: Safe-Level-Synthetic Minority Over-Sampling TEchnique for Handling the Class Imbalanced Problem. *ISCIT*, 570-575.

Chawla N.V. et al (2002) SMOTE: Synthetic Minority Over-sampling Technique. *J. Artif. Intell. Res.*, 16, 321–357.

Haibo H. et al. (2008) ADASYN: Adaptive synthetic sampling approach for imbalanced learning. *IEEE*, 1322–8.

Han H. et al. (2005) Borderline-SMOTE: A New Over-Sampling Method in Imbalanced Data Sets Learning. *ICIC*, 878–87.

Jiang M. et al. (2008) uShuffle: A useful tool for shuffling biological sequences while preserving the k-let counts. *BMC Bioinformatics*, 9, 192.

Langmead B. et al. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biol.*, 10, R25.

Lorenz R. et al. (2011) ViennaRNA Package 2.0. *Algorithms Mol. Biol.*, 6, 26.

Quinlan A.R. and Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26, 841–2.

Rivera W.A. (2017) Noise Reduction A Priori Synthetic Over-Sampling for class imbalanced data sets. *Inf. Sci. (Ny)*, 408, 146–61.

Xue C. et al. (2005) Classification of real and pseudo microRNA precursors using local structure-sequence features and support vector machine. *BMC Bioinformatics*, 6, 310.

Zhang H. and Li, M. (2014) RWO-Sampling: A random walk over-sampling approach to imbalanced data classification. *Inf. Fusion*, 20, 99–116.



&nbsp;

## Contact

emilio.marmol@cragenomica.es

&nbsp;

## Notes

- (01/18/2019) Matrix calculation bug was reported for some UNIX systems. New eMIRNA.features module was successfully tested and updated accordingly.

- (01/11/2019) The UNAfold software seems to be no longer available for free download. Provided this setback, features depending on UNAfold melt functions were removed from eMIRNA.Features module. SVM algorithm performance assesment reported no appreciable drawbacks due to UNAfold features removal.



[1]:https://academic.oup.com/bioinformatics/article/34/4/541/4222633
[2]:https://academic.oup.com/bioinformatics/article/26/3/411/215002
[3]:https://academic.oup.com/nar/article/43/21/e140/2468096
[4]:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/
[5]:https://igraph.org/
[6]:https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26
[7]:https://academic.oup.com/bioinformatics/article/26/6/841/244688
[8]:https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25
[9]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-310
[stringr]:https://CRAN.R-project.org/package=stringr
[seqinr]:https://CRAN.R-project.org/package=seqinr
[Biobase]:https://bioconductor.org/packages/release/bioc/html/Biobase.html
[scales]:https://CRAN.R-project.org/package=scales
[miRNAss]:https://CRAN.R-project.org/package=miRNAss
[PCIT]:https://CRAN.R-progect.org/package=PCIT
[NOISeq]:https://bioconductor.org/packages/release/bioc/html/NOISeq.html
[edgeR]:https://bioconductor.org/packages/release/bioc/html/edgeR.html
[igraph]:https://CRAN.R-project.org/package=igraph
[RNAfold]:https://www.tbi.univie.ac.at/RNA/
[BEDTools v2.27.0]:https://bedtools.readthedocs.io/en/latest/
[Bowtie]:https://mcardle.wisc.edu/mprime/help/bowtie/manual.html
[Fasta_ushuffle]:https://github.com/agordon/fasta_ushuffle
[Ensembl repositories]:http://www.ensembl.org/info/data/ftp/index.html
[HextractoR]:https://CRAN.R-project.org/package=HextractoR
[CD-Hit]:http://weizhong-lab.ucsd.edu/cdhit_suite/cgi-bin/index.cgi
[FASTX-Toolkit]:http://hannonlab.cshl.edu/fastx_toolkit/index.html
