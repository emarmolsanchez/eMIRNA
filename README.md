# eMIRNA

eMIRNA is a comprehensive and user-friendly pipeline for predicting and annotating the presence of known and novel microRNAs. This document is intended to give a technical supplementary description about how to run the eMIRNA pipeline through a detailed explanation of all the modules that form part of this tool.

If you implement this tool for miRNA prediction within your miRNA research projects, please cite:

[Mármol-Sánchez E. et al. (2019) Discovery and annotation of novel microRNAs in the porcine genome by using a semi-supervised transductive learning approach. *Genomics*, 112:2107-2118. doi: 10.1016/j.ygeno.2019.12.005.]

The eMIRNA pipeline is under active development, if you find any problem, doubt or bug while running it, please contact emilio.marmol@scilifelab.se.

&nbsp;
&nbsp;
&nbsp;


# Index

- [Introduction](https://github.com/emarmolsanchez/eMIRNA/#introduction)
- [Prerrequisites](https://github.com/emarmolsanchez/eMIRNA/#prerequisites)
- [miRNA Prediction](https://github.com/emarmolsanchez/eMIRNA/#miRNA-prediction)
    - [Positive, Negative and Unlabeled Data sets](https://github.com/emarmolsanchez/eMIRNA/#positive-negative-and-unlabeled-data-sets)
    - [eMIRNA.Filter](https://github.com/emarmolsanchez/eMIRNA/#emirnafilter)
    - [eMIRNA.Features](https://github.com/emarmolsanchez/eMIRNA/#emirnafeatures)
    - [eMIRNA.Hunter](https://github.com/emarmolsanchez/eMIRNA/#emirnahunter)
    - [eMIRNA.Structural.Pscore](https://github.com/emarmolsanchez/eMIRNA/#emirnastructuralpscore)
    - [eMIRNA.Predict](https://github.com/emarmolsanchez/eMIRNA/#emirnapredict)
    - [eMIRNA.Refiner](https://github.com/emarmolsanchez/eMIRNA/#emirnarefiner)
- [Functional Annotation](https://github.com/emarmolsanchez/eMIRNA/#functional-annotation)
    - [eMIRNA.Target](https://github.com/emarmolsanchez/eMIRNA/#emirnatarget)
    - [eMIRNA.Network](https://github.com/emarmolsanchez/eMIRNA/#emirnanetwork)
    - [eMIRNA.RIF](https://github.com/emarmolsanchez/eMIRNA/#emirnarif)
- [References](https://github.com/emarmolsanchez/eMIRNA/#references)
- [Contact](https://github.com/emarmolsanchez/eMIRNA/#contact)
- [Notes](https://github.com/emarmolsanchez/eMIRNA/#notes)

&nbsp;
&nbsp;

# Introduction

The eMIRNA pipeline makes use of a Machine Learning approach based on semi-supervised transductive Graph-based algorithm, as reported by Yones *et al.* (2018) [[1]], in order to assess whether putative candidate sequences can be predicted as novel miRNA genes. Additionally, target interactions with mRNA genes can be inferred to functionally annotate the novel miRNAs previously identified, making use of a system biology network approach.

&nbsp;

# Prerequisites

The following R libraries are required for running the eMIRNA pipeline:
+ [stringr]
+ [seqinr]
+ [Biobase]
+ [scales]
+ [PRROC]
+ [ROCR]
+ [dplyr]
+ [miRNAss] [[1]]
+ [PCIT] [[2]]
+ [NOISeq] [[3]]
+ [edgeR] [[4]]
+ [igraph]

The following software tools are required for running the eMIRNA pipeline:
+ [RNAfold] [[5]]
+ [BEDTools v2.27.0] [[6]]
+ [Bowtie] [[7]]
+ [Fasta_ushuffle]
+ [SeqKit Toolkit] [[8]] 


All executables should be stored at computer `$PATH` in order to be run properly (Commonly located at `/usr/bin/` or `/usr/local/bin/` in UNIX systems).

&nbsp;

# miRNA Prediction

The following modules are implemented for predicting novel miRNA genes in our genome of interest. All required steps are thoroughly explained and exemplified as follows:

&nbsp;

![alt text](https://github.com/emarmolsanchez/eMIRNA/blob/master/bin/eMIRNA_flowchart1.jpg)

**(1)** Positive, negative and unlabeled data are filtered based on size and secondary folding structure and a set of features is extracted for each sequence. **(2)** Mature miRNA sequences from small RNA-Seq data or related reference species are mapped against the selected genome assembly and elongated to reconstruct putative pre-miRNA candidates. **(3)** Candidate precursors are filtered based on size and secondary folding structure and a set of features is extracted for each candidate sequence. Optionally, sequences showing unstable secondary structure are removed. **(4)** Candidate sequences are embedded in the semi-supervised transductive classifier and a list of putative miRNAs is predicted. **(5)** Predicted miRNAs are either assigned to already annotated miRNA loci in the selected reference assembly or classified as putative novel miRNA genes.

&nbsp;

## Positive, Negative and Unlabeled Data sets

Running the **eMIRNA** Classifier requires two FASTA files with **Positive** and **Negative** sequences.

The **Positive Sequences** must correspond to those sequences annotated as microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts can be downloaded from the [Ensembl repositories].

The **Negative Sequences** must correspond to non-coding sequences other than microRNA genes in the available Reference Genome for the species under study. GTF annotation and FASTA files for corresponding transcripts can be downloaded from the [Ensembl repositories].

Additionally, other hairpin-like sequences can be extracted from the reference Genome, in order to increase the variety and diversity of sequences to be included during graph reconstruction, thus allowing a better topological adjustment of positive and negative categories. The [HextractoR] package can be used for generating a set of hairpin-like sequences from any available genome assembly. Please be aware that quering the whole genome can be extremely time consuming and resource intensive. We recommend establishing random blocks (1-10 Mb) within each genomic chromosome. As no prior knowledge of the identity of randomly extracted hairpins would available, they will be set as **Unlabeled sequences**.

Once positive, negative and unlabeled sequences are available, we recomend using a identity-by-sequence filtering step, in order to remove redundant sequences that may occur in any of the categories. The [CD-Hit] suite is a good resource for implementing sequence-based repetitve elements removal.

In the event that no Reference Genome or no good microRNA or non-coding transcripts are available for downloading, we strongly recommend choosing sequences from the closest phylogenetically related reference species with available genome annotation for running the classification procedure, otherwise the results can suffer from low reliability.

Positive, Negative and Unlabeled datasets must be in linear FASTA format. Should you have multilinear FASTA files, they must be converted to linear FASTA files. Users can implement the following perl command:

`perl -pe '/^>/ ? print "\n" : chomp' in.fa | tail -n +2 > out.fa`

where `in.fa` corresponds to multilinear FASTA file, and `out.fa` is the resulting linearized FASTA file ready to use.

&nbsp;

## eMIRNA.Filter

The first eMIRNA module makes use of previous Positive, Negative and Unlabeled FASTA files to perform an initial filtering process based on sequence length. Typically, microRNA genes range from 50 to 150 nucleotides long. Our first aim would be to filter the selected sequences based on expected microRNA genes length. We will apply this function to each of the FASTA files. The Positive sequences should not experience any filtering upon this process, if correctly generated. For Negative and Unlabeled sequences, all long non-coding hairpin-like sequences will be removed, retaining only those sequences resembling microRNA genes in length, according to established thresholds. 

Next, this module estimates the secondary folding structure of selected filtered sequences, thus filtering out all candidates that do not ressemble a pre-miRNA hairpin-like secondary structure. The eMIRNA.Filter function will make use of [RNAfold] software [[5]] to calculate the estimated secondary folding structure, which should be available in your computer `$PATH` to be correctly executed. Typically, microRNA genes have a characteristic secondary structure, composed by two stems joined by complementarity and one terminal loop, forming a hairpin-like secondary structure. Some bubbles or bulges can appear within the two stems, belonging to non-paired nucleotides in the sequence.

This function requires four arguments:

+ PATH to Positive or Negative FASTA file.
+ String with desired output prefix name.
+ Lower length filtering threshold.
+ Upper length filtering threshold.

We recommend setting 50 nucleotides as lower length threshold, and 150 for the upper, but users can define their own preferred thresholds.

Example of usage:

```r

eMIRNA.Filter.by.Size("PATH_to_Positive_FASTA", "Pos", 50, 150)

eMIRNA.Filter.by.Size("PATH_to_Negative_FASTA", "Neg", 50, 150)

eMIRNA.Filter.by.Size("PATH_to_Unlabeled_FASTA", "Unlab", 50, 150)

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

Once the eMIRNA.Features has run, a new folder named `Features_Results/` will be created inside `eMIRNA/`, in which a .txt file called `Pos/Neg/Unlab.txt` will be generated with the results of running the function.

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

The eMIRNA.Hunter script implements Bowtie [[7]] for the alignment of putative mature microRNA sequences against the reference assembly selected by the user, reconstructing pre-miRNA sequence candidates and generating a FASTA and BED files for the candidates to be subsequently classified.

Users should provide a properly collapsed FASTA file with small RNA-Seq sequences from canonical FASTQ sequence files, or a list of unique annotated mature miRNA sequences in FASTA format. The FASTQ files should be quality-check filtered and sequencing adaptor trimmed before running any available collapser tool, e.g. FASTQ collapser from [FASTX-Toolkit] for collapsing FASTQ files into FASTA files with uniquely represented sequences. We encourage to perform a pre-filtering process of the collapsed FASTA file to retain sequences between 18-25 nucleotides in length, corresponding to the average length of mature miRNAs.

This module requires seven arguments:

+ PATH to Reference Genome Bowtie Index from your species of interest.
+ PATH to FASTA file of selected reference organism mature microRNAs.
+ PATH to desired output.
+ Desired output prefix name.
+ Upwards number of bases for pre-miRNA reconstruction
+ Backwards number of bases for pre-miRNA reconstruction
+ Running mode (homology or *de novo*)

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Hunter Usage Instructions:
eMIRNA.Hunter [options]
Input:
 -r                                     PATH to Species of interest Genome Bowtie Index
 -f                                     PATH to Reference FASTA file
 -o                                     PATH to desired output folder
 -x                                     Desired Name string for output files
 -u                                     Upwards number of bases for pre-miRNA reconstruction (60-80 bp recommended)
 -b                                     Backwards number of bases for pre-miRNA reconstruction (1-10 bp recommended)
 -m                                     eMIRNA.Hunter Mode (homology or denovo)
 -h                                     Display help page
 
Output:
 <x>.sam                                SAM file output from Bowtie miRNA alignment
 <x>.log                                LOG file output from Bowtie miRNA alignment
 <x>_miRNAs.bed                         BED file output with candidate pre-miRNAs for prediction
 <x>_miRNAs.fa                          FASTA file output with candidate pre-miRNAs for prediction
 <x>_miRNAs_corrected.bed               BED file output with motif corrected candidate pre-miRNAs for prediction
 <x>_miRNAs_corrected.fa                FASTA file output with motif corrected candidate pre-miRNAs for prediction
 
 ```
We recommend defining a range between 60-80 bp and 15-30 bp for upwards and backwards elongation of the putative mature miRNA candidate, respectively. 

For achieving a successful cross-species alignment, it is very important that mature microRNA sequences FASTA file from reference organism are in DNA code, with Ts for Thymine and no Us for Uracil in RNA code. Please make sure that your mature microRNA sequences are in DNA code, otherwise the alignment process will fail.

For generating Bowtie Index for your Reference Genome, please refer to [Bowtie] Manual.

Example of usage:

```
bash eMIRNA.Hunter -r PATH_to_Bowtie_Index -f Small-RNAseq_collapsed_fastq.fa -o PATH_to_output_folder -x Candidates -u 60 -b 15 -m denovo
```

&nbsp;

After successfully running eMIRNA.Hunter script, six files will have been created at predefined output PATH:

+ SAM file with aligned sequences.
+ Log file with alignment statistics.
+ BED file with positions of reconstructed pre-miRNA candidates.
+ FASTA file with reconstructed pre-miRNA candidates.
+ BED file with motif corrected positions of reconstructed pre-miRNA candidates.
+ FASTA file with motif corrected reconstructed pre-miRNA candidates.

Once the module has run, users can optionally perform an additional filter on generated FASTA pre-miRNA candidates according to their secondary folding structural stability by using the eMIRNA.Structural.Pscore module. We only recommend to implement this step in the event that the number of candidate sequences to be analyzed do not surpass some hundreds or few thousands, as computing costs and time requeried for running multiple iterations for each sequence can escalate exponentially.

We recommend using the motif corrected FASTA files for subsequent steps, taking into consideration that not all miRNAs would be processed following motif detection and thus some novel candidates may be missed. On the contrary, a much less accurate positioning for pre-miRNA candidates will be estimated and fewer successfully detected novel miRNA candidates should be expected.

Next, candidate sequences derived from eMIRNA.Hunter directly or after eMIRNA.Structural.Pscore filtering, must be pre-processed as previously described. The eMIRNA.Filter module could be used for this purpose. Subsequently, users must process these sequences with the eMIRNA.Features module, in order to obtain a Feature matrix representing those candidate sequences that will then be subjected to classification.

&nbsp;

## eMIRNA.Structural.Pscore

The eMIRNA.Structural.Pscore module implements a n-randomization of provided sequences while maintaining *k*-let counts as described by Jiang *et al*. 2008 [[10]], using the [Fasta_ushuffle] package, which must be downloaded, compiled and stored at computer `PATH` in order to be run properly. Please be aware that this filtering step is optional and only useful when a small amount of candidate sequences need to be tested.

This module requires five arguments:

+ PATH to FASTA file of putative novel miRNA candidates generated by eMIRNA.Hunter
+ String with desired output prefix name.
+ Number of iterations to perform.
+ *P*-value threshold for structural stability filtering.
+ Boolean for FASTA filtering based on structural *P*-value (TRUE/FALSE).

Example of usage:

```r

eMIRNA.Structural.Pscore("~/eMIRNA/Structural_Results/Candidates_miRNAs_corrected.fa", "Candidates", iterate=100, threshold=0.1, filter=TRUE)

```

By default, eMIRNA.Structural.Pscore will perform 100 random shuffling iterations over each provided sequence. Users can set their desired number of iterations but should be aware of computing times required for iterating and folding of secondary structures for each sequence. As computing costs can exponentially increase with higher number of iterations, we encourage users to set their desired range of iterations between 100 and 1000, depending on the number of candidate sequences to be analyzed.

Once the eMIRNA.Structural.Pscore has run, a new .txt file called `Candidates_Structural_Pscore.txt` will be generated at `/Structural_Results` folder, containing estimated MFE *P*score for each candidates sequence. Besides, a new filtered FASTA file will be also generated at `/Structural_Results` folder, named `Candidates_filtered_Pscore.fa`, containing only those sequence candidates with structural *P*-scores < 0.1 (or the corresponding defined threshold).

The resulting filtered FASTA file can then be subjected to eMIRNA.Filter and eMIRNA.Features modules in order to obtain a representative feature matrix for each candidate sequence to be classified.

&nbsp;
 
## eMIRNA.Predict

The eMIRNA.Predict module aims to perform microRNA classification by making use of a semi-supervised transductive learning approach [[1]]. Candidates sequences generated by eMIRNA.Hunter and subsequently filtered and processed to generate a representative feature matrix must be used, as well as Positive, Negative and Unlabeled feature matrices previously generated and loaded into the R environment and called `Candidates`, `Pos`, `Neg` and `Unlab`, respectively.

This module requires three compulsory arguments and one additional argument if available:

+ Positive feature matrix.
+ Negative feature matrix.
+ Feature matrix representing candidate sequences to evaluate.
+ Optionally, users can included Unlabeled feature matrix if available.

Example of usage:

```r

Prediction <- eMIRNA.Predict(Pos, Neg, Unlab, target=Candidates)

```

Once the eMIRNA.Predict has run, an object of class `eMIRNA` will be saved at `Prediction`. Inside the object, users will find the following elements:

+ Performance: Vector of elements with estimated performance metrics (Sensitivity, Specifictiy, Accuracy, F1-Score, AUROC and AUPR) after running the prediction algorithm.
+ Prediction: Data.frame with candidate sequences IDs, Odd score for the probability of being a miRNA and predicted class (miRNA/Other).
+ ROC: Receiving Operating Characteristic curve data.frame with all necessary fields to plot ROC curve estimate.
+ PR: Precision-Recall curve data.frame with all necessary fields to plot PR curve estimate.

Users cand then save the Prediction data.frame at their preferred location.

```r

write.table(Prediction@Prediction, "~/eMIRNA/Prediction_results.txt", col.names=NA, quote=F, row.names=1, sep="\t")

```

&nbsp;

## eMIRNA.Refiner

After generating a list of putative microRNA candidates by implementing the eMIRNA pipeline, users should process these new candidates in order to assing them to already annotated miRNA loci or to novel miRNA gene candidates. The eMIRNA.Refiner module has been specifically designed to cover this task.

This module requires six arguments:

+ PATH to GTF annotation file from the species of interest.
+ PATH to data.frame table of candidate sequences generated by eMIRNA.Predict module.
+ PATH to BED file with positions of reconstructed pre-miRNA candidates by eMIRNA.Hunter.
+ PATH to FASTA file with sequences of reconstructed pre-miRNA candidates filtered by eMIRNA.Structural.Pscore (Optionally, raw FASTA file from eMIRNA.Hunter can be used).
+ PATH to desired output.
+ Desired output prefix name.

A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Refiner Usage Instructions:
eMIRNA.Refiner [options]
Input:
  -g                                             PATH to Species of interest Genome GTF Annotation file
  -p                                             PATH to Table of Predicted microRNA candidates by eMIRNA.Predict
  -b                                             PATH to BED microRNAs output file from eMIRNA.Hunter
  -f                                             PATH to FASTA microRNAs output file from eMIRNA.Structural.Pscore (or eMIRNA.Hunter)
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
bash eMIRNA.Refiner -g Sscrofa11.1.97.gtf -p Prediction_results.txt -b Candidates_miRNAs_corrected.bed -f Candidates_filtered_Pval.fa -o PATH_to_output_folder -x Candidates

```

After successfully running the eMIRNA.Refiner script, a BED file will have been created at predefined output `PATH`, containing the most relevant putative novel non-annotated microRNAs, and their estimated positions, as well as FASTA file with the corresponding novel candidate sequences. Besides, a BED file containing already annotated detected candidates will be also generated.

&nbsp;

# Functional Annotation

Users can further infer the functional interactions putatively occurring between the novel predicted miRNA candidates and other mRNA target genes expressed in their experimental conditions.

&nbsp;

![alt text](https://github.com/emarmolsanchez/eMIRNA/blob/master/bin/eMIRNA_flowchart2.jpg)

**(6)** mature miRNA and 3'-UTR sequences are used for predicting miRNA-to-mRNA target interactions. **(7)** Expression data from miRNA and mRNA genes belonging to the same experimental conditions are used to predict significant miRNA-to-mRNA interactions based on a Partial Correlations and Information Theory (PCIT) approach. **(8)** The Regulatory Impact Factor (RIF) of each considered miRNA is calculated based on mRNA and miRNA expression data, as well as with significant PCIT interactions.

&nbsp;

## eMIRNA.Target

The mechanism of action of miRNAs within the cell metabolism has been thoroughly described in previous reports [[11]]. Usually, seed portions of mature miRNAs (5'-2<sup>nd</sup> to -7<sup>th</sup>/8<sup>th</sup> nucleotides) are able to bind to short matching sites of 3'-UTR regions of target mRNAs by perfect or imperfect complementarity, hence triggering their degradation or difficulting their processing in the ribosomes and resulting in a downregulation of the production of proteins encoded by such targeted mRNAs.

eMIRNA.Target implements a fast seed pattern search on 3'-UTR regions from targeted mRNA transcripts, in order to detect putative mRNA targets of the novel annotated miRNAs. The [SeqKit Toolkit] [[8]] is implemented for retrieveing short sequence sites inside mRNA 3'-UTR regions of type 6mer, 7mer-A1, 7mer-m8 and 8mer.

This module requires five arguments:

+ PATH to FASTA file with 3'-UTR sequences of target mRNAs of interest.
+ PATH to FASTA file of mature miRNA sequences (annotated and/or novel).
+ Seed match pattern (6mer/7mer-A1/7mer-m8/8mer).
+ PATH to desired output.
+ Desired output prefix name.


A detailed explanation of each variable can be accessed with -h (help) option:

```
eMIRNA.Target Usage Instructions:

eMIRNA.Target [options]

Input:
  -m                                    PATH to mature miRNAs FASTA file
  -t                                    PATH to 3'-UTRs FASTA file
  -s                                    N-mer seed matching (7mer-m8 or 8mer recommended)
  -o                                    PATH to desired output folder
  -x                                    Desired Name string for output files
  -h                                    Display help page

Output:
  <x>_targeted_UTRs.txt                 Targeted 3'-UTR transcripts file"
  
```

Although any seed matching length could be defined among the four available options (6mer, 7mer-A1, 7mer-m8 and 8mer), we recommend using 7mer-m8 and/or 8mer motif search (please consider that putative interaction with 8mer match will be far more stringent than 7mer interactions), as these interaction types have been described as the most reliable when assessing the true functionality of miRNAs in regulating targeted mRNAs.

Example of usage:

```
bash eMIRNA.Target -m PATH_to_Positive_FASTA -t 3UTRs.fa -s 7mer-m8 -o PATH_to_output_folder -x Targets

```
After successfully running the eMIRNA.Target script, a `.txt` file will have been created at predefined output `PATH`, containing putative interactions between each queried miRNA and all the mRNA 3'-UTR regions provided.

Example output:

```
gene_ID                 miRNA              pattern       Start    End
ENSSSCT00000046027	ssc-miR-664-5p     CTAGCCT	 137      143
ENSSSCT00000046027	ssc-miR-383        TGCTGTG	  43       49
ENSSSCT00000046027	ssc-miR-615        GGCTCGG	  74       80
ENSSSCT00000046027	ssc-miR-236        GTGACCC	 167      173
ENSSSCT00000046027	ssc-miR-129a-3p    AAGGGCT	 410      416
ENSSSCT00000046027	ssc-miR-15b        TGCTGCT	  51       57
ENSSSCT00000046027	ssc-miR-15b        TGCTGCT	  54       60

```


&nbsp;

## eMIRNA.Network

In the event that mRNA and miRNA expression profiles from same experimental conditions are available, users can implement a system biology interaction approach to further infer the significance and reliability of miRNA-to-mRNA interactions predicted in advance by eMIRNA.Target.

To achive this purpose, we have implemented a network-oriented filtering criteria based on Partial Correlations and Information Theory ([PCIT]) approach as proposed by Reverter *et al.* (2008) [[12]]. By using first-order partial correlation coefficients estimated for each trio of genes along with an information theory approach, this tool identifies meaningful gene-to-gene nteractions. This approach aims to determine truly informative correlations between node pairs (genes in our context), once the influence of other nodes in the network has been considered.

This module recieves a total of eight arguments:

+ mRNA expression matrix (mRNA genes as rows and samples as columns).
+ miRNA expression matrix (miRNAs genes as rows and samples as columns).
+ Targets data.frame generated by eMIRNA.Target.
+ Correlation threshold (-0.5 by default).
+ Correlation inference method (pearson, spearman and kendall allowed).
+ Expression baseline threshold (1 CPM by default).
+ Sample expression ratio threshold (0.5 by default).
+ Boolean for performing normalization and expression filtering (TRUE/FALSE).

Example of usage (filtering and normalizing):

```r

Network <- eMIRNA.Network(mRNA, miRNA, Targets, cor=-0.5, type="pearson", cpm=1, percent=0.5, normalize=TRUE)

```
Example of usage (no filtering):

```r

Network <- eMIRNA.Network(mRNA, miRNA, Targets, cor=-0.5, type="pearson", normalize=FALSE)

```

Example output:

```
miRNA               mRNA                   Correlation
ssc-miR-148a-3p     ENSSSCG00000031262      -0.64
ssc-miR-148a-3p	    ENSSSCG00000014156      -0.58
ssc-miR-148a-3p	    ENSSSCG00000015334      -0.75
ssc-miR-148a-3p     ENSSSCG00000015872      -0.56
ssc-miR-148a-3p     ENSSSCG00000004332      -0.59
ssc-miR-148a-3p     ENSSSCG00000010210      -0.70

```

&nbsp;

Users should define mRNA and miRNA matrices with rownames as mRNA/miRNA names and colnames as sample names. The data.frame incorporating targets information should have mRNA names as first column.

If `normalize=FALSE` is set, users should provided already filtered, normalized and log<sub>2</sub> transformed mRNA and miRNA expression matrices. Samples names (columns) must be the same in both mRNA and miRNA expression data and the exact same number of individuals must be included so as the function can be run properly. Please be aware that [PCIT] algorithm can escalate to high time and memmory consuming requirements if a huge amount of expression data is included in mRNA and miRNA matrices. 

Tipically, the higher the number of genes, the better for [PCIT] inference. Default expression and normalization filters implemented are suited for a good representation of the miRNA-to-mRNA interaction network. For a fast running, users may focus on differentially expressed (DE) mRNA and miRNA genes, otherwise some hours should be expected for network inference with around 10 to 15 thousands of expressed genes, depending on computing resources.

&nbsp;


## eMIRNA.RIF

Once putative interactions have been inferred, users may want to estimate the relevance of each considered miRNA in regulating the expression of targeted mRNAs in their experimental conditions. For this purpose, we have implemented the eMIRNA.RIF module, which makes use of the Regulatory Impact Factor ([RIF]) algorithm described by Reverter *et al.* (2010) [[13]]. 

The [RIF] algorithm aims to identify regulator genes contributing to the observed differential expression in the analyzed contrasts. Its implementation results in two different and inter-connected RIF scores: while RIF1 score represents those transcriptional regulators that are most differentially co-expressed with the most highly abundant and highly DE genes, the RIF2 score highlights those regulators that show the most altered ability to act as predictors of the changes in the expression levels of DE genes [[13]]. Both [RIF] values capture different regulatory impact features and hence, they can be considered as two independent measurements of the putative relevance of miRNAs as gene expression regulators.

This module recieves a total of eight arguments:

+ mRNA expression matrix (mRNA genes as rows and samples as columns).
+ miRNA expression matrix (miRNAs genes as rows and samples as columns).
+ Experiment design.
+ PCIT Network data.frame from eMIRNA.Network.
+ List of DE mRNA genes.
+ Expression baseline threshold (1 CPM by default).
+ Sample expression ratio threshold (0.5 by default).
+ Boolean for performing normalization and expression filtering (TRUE/FALSE).

Expression matrices should have gene identifiers and sample identifiers as row names and column names, respectively. The exact same samples should be available for mRNA and miRNA expression data. Binary comparison is supported by [RIF] algorithm, where control group is contrasted against treatment group. Experimental design should be a data.frame format with first row having the exact same sampe identifiers than expression matrices, and second row having group identification (either control or treatment). The list of DE mRNAs should be a data.frame with one column containing all the DE mRNA genes detected by any Differential Expression analysis performed previously by users. Gene identifiers should be the same than in expression matrices.

Example of usage (filtering and normalizing):

```r

RIF <- eMIRNA.RIF(mRNA, miRNA, design, Network, DElist, cpm=1, percent=0.5, normalize=TRUE)

```
Example of usage (no filtering):

```r

RIF <- eMIRNA.RIF(mRNA, miRNA, design, Network, DElist, normalize=FALSE)

```

Example output:

```
miRNA               RIF1        RIF2
ssc-miR-32          1.7041      -0.9823
ssc-miR-136-5p	    1.3928      1.2053
ssc-miR-542-3p	    1.2969      -2.0613

 ```


Users should define mRNA and miRNA matrices with rownames as mRNA/miRNA names and colnames as sample names. The data.frame incorporating [PCIT] network information should have three columns with genes, miRNAs and correlation information, respectively. The list of DE genes should have each gene ID in one column variable.

&nbsp;
&nbsp;
&nbsp;


## References

 1. [Yones C. et al. (2018) Genome-wide pre-miRNA discovery from few labeled examples. *Bioinformatics*, 34, 541–49.]

 2. [Watson-High N.S. et al. (2010) PCIT: an R package for weighted gene co-expression networks based on partial correlation and information theory approaches. *Bioinformatics*, 26, 411-13.]

 3. [Tarazona S. et al. (2015) Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package. *Nucleic Acids Research*, 43, e140.]

 4. [Robinson M.D. et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139-40.]

 5. [Lorenz R. et al. (2011) ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26.]

 6. [Quinlan A.R. and Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26, 841–2.]

 7. [Langmead B. et al. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, 10, R25.]

 8. [Shen W. et al. (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLoS ONE*, 11, e0163962.]

 9. [Xue C. et al. (2005) Classification of real and pseudo microRNA precursors using local structure-sequence features and support vector machine. *BMC Bioinformatics*, 6, 310.]

10. [Jiang M. et al. (2008) uShuffle: A useful tool for shuffling biological sequences while preserving the k-let counts. *BMC Bioinformatics*, 9, 192.]

11. [Bartel D.P. (2018) Metazoan microRNAs. *Cell*, 173, 20-51.]

12. [Reverter A. et al. (2008) Combining partial correlation and an information theory approach to the reversed engineering of gene co-expression networks. *Bioinformatics*, 24, 2491-97.]

13. [Reverter A. et al. (2010) Regulatory impact factors: unraveling the transcriptional regulation of complex traits from expression data. *Bioinformatics*, 26, 896-904.]

&nbsp;

## Contact

emilio.marmol@cragenomica.es

&nbsp;

## Notes

- (11/11/2019) Update for functional annotation has been added.

- (06/15/2019) Prediction algorithm has been updated.

- (01/18/2019) Matrix calculation bug was reported for some UNIX systems. New eMIRNA.features module was successfully tested and updated accordingly.

- (01/11/2019) The UNAfold software seems to be no longer available for free download. Provided this setback, features depending on UNAfold melt functions were removed from eMIRNA.Features module. SVM algorithm performance assesment reported no appreciable drawbacks due to UNAfold features removal.



[1]:https://academic.oup.com/bioinformatics/article/34/4/541/4222633
[2]:https://academic.oup.com/bioinformatics/article/26/3/411/215002
[3]:https://academic.oup.com/nar/article/43/21/e140/2468096
[4]:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/
[5]:https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26
[6]:https://academic.oup.com/bioinformatics/article/26/6/841/244688
[7]:https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25
[8]:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962
[9]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-310
[10]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192
[11]:https://www.sciencedirect.com/science/article/pii/S0092867418302861?via%3Dihub
[12]:https://academic.oup.com/bioinformatics/article/24/21/2491/192682
[13]:https://academic.oup.com/bioinformatics/article/26/7/896/212064
[stringr]:https://CRAN.R-project.org/package=stringr
[seqinr]:https://CRAN.R-project.org/package=seqinr
[Biobase]:https://bioconductor.org/packages/release/bioc/html/Biobase.html
[scales]:https://CRAN.R-project.org/package=scales
[PRROC]:https://CRAN.R-project.org/package=PRROC
[ROCR]:https://CRAN.R-project.org/package=ROCR
[dplyr]:https://CRAN.R-project.org/package=dplyr
[miRNAss]:https://CRAN.R-project.org/package=miRNAss
[PCIT]:https://github.com/nathanhaigh/pcit
[RIF]:https://academic.oup.com/bioinformatics/article/26/7/896/212064
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
[SeqKit Toolkit]:https://bioinf.shenwei.me/seqkit/
[Yones C. et al. (2018) Genome-wide pre-miRNA discovery from few labeled examples. *Bioinformatics*, 34, 541–49.]:https://academic.oup.com/bioinformatics/article/34/4/541/4222633
[Watson-High N.S. et al. (2010) PCIT: an R package for weighted gene co-expression networks based on partial correlation and information theory approaches. *Bioinformatics*, 26, 411-13.]:https://academic.oup.com/bioinformatics/article/26/3/411/215002
[Tarazona S. et al. (2015) Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package. *Nucleic Acids Research*, 43, e140.]:https://academic.oup.com/nar/article/43/21/e140/2468096
[Robinson M.D. et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26, 139-40.]:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/
[Lorenz R. et al. (2011) ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26.]:https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26
[Quinlan A.R. and Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26, 841–2.]:https://academic.oup.com/bioinformatics/article/26/6/841/244688
[Langmead B. et al. (2009) Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, 10, R25.]:https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25
[Shen W. et al. (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLoS ONE*, 11, e0163962.]:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962
[Xue C. et al. (2005) Classification of real and pseudo microRNA precursors using local structure-sequence features and support vector machine. *BMC Bioinformatics*, 6, 310.]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-310
[Jiang M. et al. (2008) uShuffle: A useful tool for shuffling biological sequences while preserving the k-let counts. *BMC Bioinformatics*, 9, 192.]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192
[Bartel D.P. (2018) Metazoan microRNAs. *Cell*, 173, 20-51.]:https://www.sciencedirect.com/science/article/pii/S0092867418302861?via%3Dihub
[Reverter A. et al. (2008) Combining partial correlation and an information theory approach to the reversed engineering of gene co-expression networks. *Bioinformatics*, 24, 2491-97.]:https://academic.oup.com/bioinformatics/article/24/21/2491/192682
[Reverter A. et al. (2010) Regulatory impact factors: unraveling the transcriptional regulation of complex traits from expression data. *Bioinformatics*, 26, 896-904.]:https://academic.oup.com/bioinformatics/article/26/7/896/212064
[Mármol-Sánchez E. et al. (2019) Discovery and annotation of novel microRNAs in the porcine genome by using a semi-supervised transductive learning approach. *Genomics*, 112:2107-2118. doi: 10.1016/j.ygeno.2019.12.005.]:https://www.sciencedirect.com/science/article/pii/S0888754319304884?via%3Dihub
