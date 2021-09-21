#########################################
#                                       #
#          eMIRNA R functions           #
#                                       #
#             EMS. 2019                 #
#                                       #
#########################################


### eMIRNA.Filter ###

eMIRNA.Filter <- function(file, prefix, a, b){
  suppressMessages(require(Biobase))
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "Filter_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/Filter_Results/"
  setwd(workdir)
  message("Filtering Sequences...")
  
  #Checking Sequences and Filtering by size (50 to 150 nt)
  File0 <- unlist(lapply(file, readLines))
  n0 <- length(File0)
  ID0 <- File0[seq(1,n0,2)]
  Sequence0 <- File0[seq(2,n0,2)]
  Sequence0.split <- strsplit(Sequence0, "")
  Sequence0.length <- sapply(Sequence0.split, function(x) length(x))
  Index0.filter <- which(Sequence0.length >= a & Sequence0.length <= b)
  ID0.filter <- ID0[Index0.filter]
  Sequence0.filter <- Sequence0[Index0.filter]
  File0.filter.size <- c(rbind(ID0.filter, Sequence0.filter))
  name <- paste0(prefix, "_temp")
  write(File0.filter.size, name)
  
  #Checking Filtered Sequences and filtering by nloops
  n0 <- length(File0.filter.size)
  ID0 <- File0.filter.size[seq(1,n0,2)]
  Sequence0 <- File0.filter.size[seq(2,n0,2)]
  command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", workdir, name)
  RNAfold1 <- system(command1, intern=TRUE)
  n <- length(RNAfold1)
  SecondaryStrc1 <- RNAfold1[seq(3,n,7)]
  SecondaryStrc <- gsub( " .*$", "", SecondaryStrc1)
  Nloop1 <- strsplit(SecondaryStrc, "\\((?=(\\.+\\)))", perl = TRUE)
  Nloop <- listLen(Nloop1) - 1
  Index0.nloop <- which(Nloop == 1)
  ID0.nloop <- ID0[Index0.nloop]
  Sequence0.nloop <- Sequence0[Index0.nloop]
  File0.filter.nloop <- c(rbind(ID0.nloop, Sequence0.nloop))
  name <- paste0(prefix, "_filtered.fa")
  write(File0.filter.nloop, name)
  
  unlink("*.ps", recursive=T)
  unlink("*temp", recursive=T)
  
}


### eMIRNA.Features ###

eMIRNA.Features <- function(file, prefix){
  suppressMessages(require(seqinr))
  suppressMessages(require(stringr))
  suppressMessages(require(scales))
  suppressMessages(require(Biobase))
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "Feature_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/Feature_Results/"
  setwd(workdir)
  message("Calculating Features...")
  if(substr(file, nchar(file)-12+1, nchar(file)) == '_filtered.fa'){
    command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", file)
    RNAfold1 <- system(command1, intern=TRUE)
    
    ###########################################
    ############ Sequence Features ############
    ###########################################
    
    #Triplet estimation from SVM-Triplet (Triplets)
    command2 <- paste0("RNAfold --noPS -i ", file, " > RNAfold_pred.txt")
    system(command2)
    command2.1 <- "1_check_query_content.pl RNAfold_pred.txt RNAfold_pred_checked.txt"
    system(command2.1)
    command2.2 <- "2_get_stemloop.pl RNAfold_pred_checked.txt RNAfold_pred_stemloop.txt 5"
    system(command2.2)
    command2.3 <- "3_step_triplet_coding_for_queries.pl RNAfold_pred_stemloop.txt Triplets.txt"
    system(command2.3)
    Triplets1 <- "Triplets.txt"
    Triplets1 <- unlist(lapply(Triplets1, readLines))
    Triplets.split1 <- strsplit(Triplets1, "\\s")
    Triplets1 <- as.numeric(as.character(unlist(Triplets.split1)))
    Triplets.along <- seq_along(Triplets1)
    Triplets <- split(Triplets1, ceiling(Triplets.along/32))
    Triplets <- do.call(rbind, Triplets)
    colnames(Triplets) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                            "14","15","16","17","18","19","20","21","22","23","24","25",
                            "26","27","28","29","30","31","32")
    T1 <- as.vector(Triplets[,1])
    T2 <- as.vector(Triplets[,2])
    T3 <- as.vector(Triplets[,3])
    T4 <- as.vector(Triplets[,4])
    T5 <- as.vector(Triplets[,5])
    T6 <- as.vector(Triplets[,6])
    T7 <- as.vector(Triplets[,7])
    T8 <- as.vector(Triplets[,8])
    T9 <- as.vector(Triplets[,9])
    T10 <- as.vector(Triplets[,10])
    T11 <- as.vector(Triplets[,11])
    T12 <- as.vector(Triplets[,12])
    T13 <- as.vector(Triplets[,13])
    T14 <- as.vector(Triplets[,14])
    T15 <- as.vector(Triplets[,15])
    T16 <- as.vector(Triplets[,16])
    T17 <- as.vector(Triplets[,17])
    T18 <- as.vector(Triplets[,18])
    T19 <- as.vector(Triplets[,19])
    T20 <- as.vector(Triplets[,20])
    T21 <- as.vector(Triplets[,21])
    T22 <- as.vector(Triplets[,22])
    T23 <- as.vector(Triplets[,23])
    T24 <- as.vector(Triplets[,24])
    T25 <- as.vector(Triplets[,25])
    T26 <- as.vector(Triplets[,26])
    T27 <- as.vector(Triplets[,27])
    T28 <- as.vector(Triplets[,28])
    T29 <- as.vector(Triplets[,29])
    T30 <- as.vector(Triplets[,30])
    T31 <- as.vector(Triplets[,31])
    T32 <- as.vector(Triplets[,32])
    
    StemF <- "RNAfold_pred_stemloop.txt"
    StemF <- unlist(lapply(StemF, readLines))
    StemF.ID <- grep(">", StemF, value=TRUE)
    StemF.ID <- gsub(">", "", StemF.ID)
    StemF.ID <- strsplit(StemF.ID, "_")
    StemF.ID <- lapply(StemF.ID, "[", -2)
    StemF.ID.index <- as.numeric(unlist(lapply(StemF.ID, "[", -2)))
    
    #ID
    n <- length(RNAfold1)
    ID <- RNAfold1[seq(1,n,7)]
    ID <- gsub('>', '', ID)
    ID <- ID[StemF.ID.index]
    
    #Sequence
    Sequence <- RNAfold1[seq(2,n,7)]
    Sequence <- Sequence[StemF.ID.index]
    
    #Length of Sequences (Length)
    Length <- nchar(Sequence)
    
    #G+C ratio (GC)
    nG <- str_count(Sequence, "G")
    nC <- str_count(Sequence, "C")
    GC <- (nG + nC) / Length
    
    #G/C ratio (G.Cr)
    G.Cr <- nG / nC
    
    #A+U/G+C ratio (AU.GCr)
    nA <- str_count(Sequence, "A")
    nU <- str_count(Sequence, "U")
    AU.GCr <- (nA + nU) / (nG + nC)
    
    #Base proportions Ratios (Ar, Ur, Gr, Cr)
    Ar <- nA / Length
    Ur <- nU / Length
    Gr <- nG / Length
    Cr <- nC / Length
    
    
    #Dinucleotide Ratios (DNr)
    AAr <- str_count(Sequence, "AA")/ Length
    GGr <- str_count(Sequence, "GG")/ Length
    CCr <- str_count(Sequence, "CC")/ Length
    UUr <- str_count(Sequence, "UU")/ Length
    AGr <- str_count(Sequence, "AG")/ Length
    ACr <- str_count(Sequence, "AC")/ Length
    AUr <- str_count(Sequence, "AU")/ Length
    GAr <- str_count(Sequence, "GA")/ Length
    GCr <- str_count(Sequence, "GC")/ Length
    GUr <- str_count(Sequence, "GU")/ Length
    CAr <- str_count(Sequence, "CA")/ Length
    CGr <- str_count(Sequence, "CG")/ Length
    CUr <- str_count(Sequence, "CU")/ Length
    UAr <- str_count(Sequence, "UA")/ Length
    UGr <- str_count(Sequence, "UG")/ Length
    UCr <- str_count(Sequence, "UC")/ Length
    
    
    ###################################################
    ########### Secondary Structure Features ##########
    ###################################################
    
    #RNAfold Secondary Structure (SecondaryStrc)
    SecondaryStrc1 <- RNAfold1[seq(3,n,7)]
    SecondaryStrc1 <- SecondaryStrc1[StemF.ID.index]
    SecondaryStrc <- gsub( " .*$", "", SecondaryStrc1)
    
    #Pairing Probabilities Structure (PairProbStrc)
    PairProbStrc1 <- RNAfold1[seq(4,n,7)]
    PairProbStrc1 <- PairProbStrc1[StemF.ID.index]
    PairProbStrc <- gsub( " .*$", "", PairProbStrc1)
    
    #RNAfold centroid structure (CentroidStrc)
    CentroidStrc1 <- RNAfold1[seq(5,n,7)]
    CentroidStrc1 <- CentroidStrc1[StemF.ID.index]
    CentroidStrc <- gsub( " .*$", "", CentroidStrc1)
    
    #Maximum Expected Accuracy Structure (MEAStrc)
    MEAStrc1 <- RNAfold1[seq(6,n,7)]
    MEAStrc1 <- MEAStrc1[StemF.ID.index]
    MEAStrc <- gsub( " .*$", "", MEAStrc1)
    
    #Hairpin length (Hl)
    nl <- length(StemF)
    Hl1 <- StemF[seq(3,nl,3)]
    Hl2 <- strsplit(Hl1, "\\((?=(\\.+\\)))", perl = TRUE)
    Hl1 <- lapply(Hl2, "[", -1)
    Hl1 <- lapply(Hl1, function(x) gsub(").*", "", as.character(x)))
    Hl3 <- lapply(Hl2, nchar)
    nseq <- length(StemF.ID.index)
    Hl3.2 <- as.list(rep(1, nseq))
    Hl3 <- mapply('+', Hl3.2, Hl3, SIMPLIFY=FALSE)
    Hl1 <- lapply(Hl1, nchar)
    Hl <- as.vector(unlist(Hl1))
    
    #5' and 3' Stem length (Steml5, Steml3)
    Hl.list <- as.list(Hl)
    list0 <- as.list(rep(0, nseq))
    Steml1 <- mapply('c', list0, Hl.list, SIMPLIFY=FALSE)
    Steml <- mapply('-', Hl3, Steml1, SIMPLIFY=FALSE)
    Steml5 <- sapply(Steml, "[", 1)
    Steml3 <- sapply(Steml, "[", 2)
    Steml3 <- Steml3 - 1
    
    #Number of basepairs in Secondary Structure (BP)
    BP <- str_count(SecondaryStrc, "\\(")
    
    #Number of basepairs in 5' and 3' Stem (BP5, BP3)
    BP1 <- unlist(lapply(Hl2, "[", 1))
    BP5 <- str_count(BP1, "\\(") + 1
    BP2 <- unlist(lapply(Hl2, "[", 2))
    BP2 <- sub("^\\.*", "", BP2)
    BP3 <- str_count(BP2, "\\)")
    
    #Number of mismatches in 5' and 3' Stem (Mism5, Mism3)
    Mism5 <- str_count(BP1, "\\.")
    Mism3 <- str_count(BP2, "\\.")
    
    #Number of bulges in 5' and 3' Stem (Bulge5, Bulge3)
    Bulge1 <- strsplit(BP1, "(\\.)+.")
    Bulge5 <- listLen(Bulge1) - 1
    Bulge2 <- strsplit(BP2, "(\\.)+.")
    Bulge3 <- listLen(Bulge2) - 1
    
    #Number of bulges by type (BulgeN1, BulgeN2, BulgeN3, BulgeN4, BulgeN5)
    BulgeN1.5 <- str_count(BP1, "\\(\\.\\(")
    BulgeN1.3 <- str_count(BP2, "\\)\\.\\)")
    BulgeN2.5 <- str_count(BP1, "\\(\\.\\.\\(")
    BulgeN2.3 <- str_count(BP2, "\\)\\.\\.\\)")
    BulgeN3.5 <- str_count(BP1, "\\(\\.\\.\\.\\(")
    BulgeN3.3 <- str_count(BP2, "\\)\\.\\.\\.\\)")
    BulgeN4.5 <- str_count(BP1, "\\(\\.\\.\\.\\.\\(")
    BulgeN4.3 <- str_count(BP2, "\\)\\.\\.\\.\\.\\)")
    BulgeN5.5 <- str_count(BP1, "\\(\\.\\.\\.\\.\\.\\(")
    BulgeN5.3 <- str_count(BP2, "\\)\\.\\.\\.\\.\\.\\)")
    BulgeN6.5 <- str_count(BP1, "\\(\\.\\.\\.\\.\\.\\.\\(")
    BulgeN6.3 <- str_count(BP2, "\\)\\.\\.\\.\\.\\.\\.\\)")
    BulgeN7.5 <- str_count(BP1, "\\(\\.\\.\\.\\.\\.\\.\\.\\(")
    BulgeN7.3 <- str_count(BP2, "\\)\\.\\.\\.\\.\\.\\.\\.\\)")
    
    #Number of A-U, C-G and G-U pairs in 1 loop sequences (AUp, GCp, GUp)
    Sequence2 <- strsplit(Sequence, "")
    SecondaryStrc2 <- strsplit(SecondaryStrc, "")
    Matches5 <- sapply(SecondaryStrc2, function(x) which(x == "("))
    Matches3 <- sapply(SecondaryStrc2, function(x) which(x == ")"))
    
    if(class(Matches5) == "matrix"){
      
      Extract5 <- lapply(seq_len(ncol(Matches5)), function(i) Matches5[,i])
      Extract5 <- mapply('[', Sequence2, Extract5)
      Extract5 <- lapply(seq_len(ncol(Extract5)), function(i) Extract5[,i])
      Extract3 <- lapply(seq_len(ncol(Matches3)), function(i) Matches3[,i])
      Extract3 <- mapply('[', Sequence2, Extract3)
      Extract3 <- lapply(seq_len(ncol(Extract3)), function(i) Extract3[,i])
      Extract3 <- lapply(Extract3, function(x) rev(x))
      Pairs <- mapply('paste0', Extract5, Extract3, SIMPLIFY=FALSE)
      CGc1 <- sapply(Pairs, function(x) sum(str_count(x, "CG")))
      GCc1 <- sapply(Pairs, function(x) sum(str_count(x, "GC")))
      GCp <- CGc1 + GCc1
      AUc1 <- sapply(Pairs, function(x) sum(str_count(x, "AU")))
      UAc1 <- sapply(Pairs, function(x) sum(str_count(x, "UA")))
      AUp <- AUc1 + UAc1
      GUc1 <- sapply(Pairs, function(x) sum(str_count(x, "GU")))
      UGc1 <- sapply(Pairs, function(x) sum(str_count(x, "UG")))
      GUp <- GUc1 + UGc1
      
    } else {
      
      Extract5 <- mapply('[', Sequence2, Matches5)
      Extract3 <- mapply('[', Sequence2, Matches3)
      Extract3 <- sapply(Extract3, function(x) rev(x))
      Pairs <- mapply('paste0', Extract5, Extract3, SIMPLIFY=FALSE)
      CGc1 <- sapply(Pairs, function(x) sum(str_count(x, "CG")))
      GCc1 <- sapply(Pairs, function(x) sum(str_count(x, "GC")))
      GCp <- CGc1 + GCc1
      AUc1 <- sapply(Pairs, function(x) sum(str_count(x, "AU")))
      UAc1 <- sapply(Pairs, function(x) sum(str_count(x, "UA")))
      AUp <- AUc1 + UAc1
      GUc1 <- sapply(Pairs, function(x) sum(str_count(x, "GU")))
      UGc1 <- sapply(Pairs, function(x) sum(str_count(x, "UG")))
      GUp <- GUc1 + UGc1
      
    }
    
    Nloop1 <- strsplit(SecondaryStrc, "\\((?=(\\.+\\)))", perl = TRUE)
    Nloop <- listLen(Nloop1) - 1
    
    ##################################################
    ######### Secondary Structure Statistics #########
    ##################################################
    
    
    
    #Minimum Free Energy (MFE)
    MFE <- as.numeric(regmatches(SecondaryStrc1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          SecondaryStrc1, perl=TRUE)))
    
    #Ensemble Free Energy (EFE)
    EFE <- as.numeric(regmatches(PairProbStrc1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          PairProbStrc1, perl=TRUE)))
    
    #Centroid Free Energy (CFE)
    CFE1 <- sub(".*? (.+)", "\\1", CentroidStrc1)
    CFE1 <- gsub("\\{", "", CFE1)
    CFE1 <- gsub("\\}", "", CFE1)
    CFE1 <- gsub(" ", "", CFE1)
    CFE1 <- unlist(strsplit(CFE1, "d="))
    n2 <- length(CFE1)
    CFE <- as.numeric(CFE1[seq(1,n2,2)])
    
    #Centroid Distance to Ensemble (CDE)
    CDE <- as.numeric(CFE1[seq(2,n2,2)])
    
    #Maximum Expected Accuracy Structure Free Energy (MEAFE)
    MEAFE1 <- sub(".*? (.+)", "\\1", MEAStrc1)
    MEAFE1 <- gsub("\\{", "", MEAFE1)
    MEAFE1 <- gsub("\\}", "", MEAFE1)
    MEAFE1 <- gsub(" ", "", MEAFE1)
    MEAFE1 <- unlist(strsplit(MEAFE1, "MEA="))
    MEAFE <- as.numeric(MEAFE1[seq(1,n2,2)])
    
    #Maximum Expected Accurary (MEA)
    MEA <- as.numeric(MEAFE1[seq(2,n2,2)])
    
    #Base pair Propensity (BPP)
    BPP <- (BP / Length)
    
    #Frequency of the MFE in ensemble (EFreq)
    EFreq1 <- RNAfold1[seq(7,n,7)]
    EFreq1 <- EFreq1[StemF.ID.index]
    EFreq1 <- unlist(strsplit(EFreq1, ";"))
    EFreq2 <- EFreq1[seq(1,n2,2)]
    #EFreq <- round(as.numeric(regmatches(EFreq2,
    #                             gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
    #                                       EFreq2, perl=TRUE))), 4)
    
    
    #Ensemble Diversity (ED)
    ED1 <- EFreq1[seq(2,n2,2)]
    ED <- as.numeric(regmatches(ED1,
                                gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", 
                                         ED1, perl=TRUE)))
    
    #Adjusted MFE (MFEadj)
    MFEadj <- (MFE / Length)
    
    #Adjusted EFE (EFEadj)
    EFEadj <- (EFE / Length)
    
    #Adjusted base pair Distance (Dadj)
    Dadj <- (CDE / Length)
    
    #Adjusted Shannon Entropy (SEadj)
    SE <- -((Ar*log2(Ar))+(Ur*log2(Ur))+(Gr*log2(Gr))+(Cr*log2(Cr)))
    SEadj <- (SE /Length)
    
    #Difference between MFE and EFE Adjusted (DiffMFE.EFE)
    DiffMFE.EFE <- ((MFE - EFE) / Length)
    
    #Ratio between Adjusted MFE and GC (MFEadj.GC)
    MFEadj.GC <- (MFEadj / GC)
    
    #Ratio between Adjusted MFE and base pairs (MFEadj.BP)
    MFEadj.BP <- (MFEadj / BP)
    
    #Adjusted MEAFE
    MEAFEadj <- (MEAFE / Length)
    
    #Adjusted ED
    EDadj <- (ED / Length)
    
    
    
    table <- as.data.frame(cbind(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14,
                                 T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25,
                                 T26, T27, T28, T29, T30, T31, T32, GC, BPP, Dadj, SEadj,
                                 DiffMFE.EFE, Length, AU.GCr,
                                 Ar, Ur, Gr, Cr, AAr, GGr, CCr, UUr, AGr, ACr, AUr,
                                 GAr, GCr, GUr, CAr, CGr, CUr, UAr, UGr, UCr, Hl, Steml5,
                                 Steml3, BP, BP5, BP3, Mism5, Mism3, Bulge5, Bulge3, BulgeN1.5, 
                                 BulgeN1.3, BulgeN2.5, BulgeN2.3, BulgeN3.5, BulgeN3.3,
                                 BulgeN4.5, BulgeN4.3, BulgeN5.5, BulgeN5.3, BulgeN6.5, BulgeN6.3,
                                 BulgeN7.5, BulgeN7.3, AUp, GCp, GUp, MFE,
                                 MFEadj, EFE, EFEadj, CFE, CDE, MEAFE, MEA,
                                 ED, MFEadj.GC, MFEadj.BP, MEAFEadj, EDadj))
    
    unlink("*.ps")
    unlink("*.txt")
    unlink("*.fa*")
    
    
    rownames(table) <- ID
    table[is.na(table)] <- 0
    
    final_table <- paste0(workdir, prefix, "_features.txt")
    write.table(table, final_table, sep="\t", quote=F, col.names=NA)
    
    return(table)
    
    
  } else {
    
    stop("Input file is not filtered FASTA file.")
    
  }
  
}


### eMIRNA.Structural.Pscore ###

eMIRNA.Structural.Pscore <- function(file, prefix, iterate=100, threshold=0.1, filter=TRUE){
  suppressMessages(require(Biobase))
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "Structural_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/Structural_Results/"
  outdir <- "~/eMIRNA/Structural_Results/"
  setwd(workdir)
  message("Calculating P-scores...")
  command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", file)
  RNAfold1 <- system(command1, intern=TRUE)
  
  #ID
  n <- length(RNAfold1)
  ID <- RNAfold1[seq(1,n,7)]
  ID <- gsub('>', '', ID)
  
  #Sequence
  Sequence <- RNAfold1[seq(2,n,7)]
  
  #RNAfold Secondary Structure (SecondaryStrc)
  SecondaryStrc1 <- RNAfold1[seq(3,n,7)]
  SecondaryStrc <- gsub( " .*$", "", SecondaryStrc1)
  
  
  #Minimum Freen Energy (MFE)
  MFE <- as.numeric(regmatches(SecondaryStrc1,
                               gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                        SecondaryStrc1, perl=TRUE)))
  
  #MFE P-value (MFE.Pval)
  command3 <- paste0("fasta_ushuffle -n ", iterate, " -k 2 -s 1234 < ", 
                     file, " >", prefix, "_iterated.fa")
  seqs <- system(command3, intern=TRUE)
  
  command4 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", prefix, "_iterated.fa")
  RNAfold2 <- system(command4, intern=TRUE)
  n4 <- length(RNAfold2)
  MFEit1 <- RNAfold2[seq(3,n4,7)]
  MFEit <- as.numeric(regmatches(MFEit1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", 
                                          MFEit1, perl=TRUE)))
  N <- iterate
  MFEit.along <- seq_along(MFEit)
  MFEit.split <- split(MFEit, ceiling(MFEit.along/N))
  R <- as.vector(sapply(Map(`<=`, MFEit.split, MFE), sum))
  MFE.Pval <- (R / (N + 1))
  
  
  unlink("*.ps")
  unlink("*.fa")
  
  MFE.Pval <- round(MFE.Pval, 6)
  table <- as.data.frame(MFE.Pval)
  rownames(table) <- ID
  
  final_table <- paste0(outdir, prefix, "_Structural_Pscore.txt")
  write.table(table, final_table, quote=F, row.names=T, col.names=NA, sep="\t")
  
  if(filter == TRUE){
    
    index_thres0 <- which(table$MFE.Pval < threshhold)
    index_thres1 <- (index_thres0*2)-1
    index_thres2 <- index_thres1+1
    index_thres <- c(matrix(c(index_thres1, index_thres2), 2, byrow = T)) 
    File0 <- unlist(lapply(file, readLines))
    File_filter <- File0[index_thres]
    name <- paste0(prefix, "_filtered_Pscore.fa")
    write(File_filter, name)
    
    return(table)
    
  } else if(filter == FALSE) {
    
    return(table)
  }
  
}


### eMIRNA.Predict ###

eMIRNA.Predict <- function(pos, neg, unlab=NULL, target){
  suppressMessages(require(miRNAss))
  suppressMessages(require(PRROC))
  suppressMessages(require(ROCR))
  suppressMessages(require(dplyr))
  message("Computing Performance...")
  
  #Performance
  
  if(is.null(unlab)==FALSE){
  Table <- rbind(pos, neg, unlab)
  labels <- c(rep(1,nrow(pos)), rep(-1,nrow(neg)), rep(0,nrow(unlab)))
  
  set.seed(50)
  labels_test <- c(rep(1,nrow(pos)), rep(-1,nrow(neg)), rep(0,nrow(unlab)))
  labels_test[sample(which(labels_test > 0),nrow(pos)*0.25)] = 0
  labels_test[sample(which(labels_test < 0),nrow(neg)*0.25)] = 0
  
  p <- miRNAss(Table, labels_test, thresholdObjective = "F1")
  
  
  index_testing <- which(labels_test[1:(nrow(pos)+nrow(neg))] == 0)
  
  p_testing <- p[index_testing]
  
  npos <- length(which(index_testing <= nrow(pos)))
  nneg <- length(index_testing) - npos
  
  TP <- length(which(p_testing[1:npos] > 0)) / npos
  FN <- length(which(p_testing[1:npos] < 0)) / npos
  SE <- TP/(TP+FN)
  TN <- length(which(p_testing[(npos+1):((npos)+(nneg))] < 0)) / nneg
  FP <- length(which(p_testing[(npos+1):((npos)+(nneg))] > 0)) / nneg
  SP <- TN/(TN+FP)
  Acc <- (TP+TN)/(TP+TN+FP+FN)
  F1 <- (2*TP)/((2*TP)+FP+FN)
  
  
  final.table <- rbind(pos, neg)
  
  class <- c(rep("miRNA", nrow(pos)), rep("Other", nrow(neg)))
  
  final.table <- as.data.frame(cbind(final.table, class))
  final.table$class <- factor(final.table$class)

  testing <- final.table[index_testing,]
  
  p_testing <- rescale(p_testing, to=c(0,1))
  pos_p <- p_testing[1:(floor(nrow(pos)*0.25))]
  neg_p <- p_testing[(floor(nrow(pos)*0.25)+1):(floor(nrow(pos)*0.25)+floor(nrow(neg)*0.25))]
  prediction <- c(pos_p, neg_p)
  test.label <- as.integer(recode(testing$class,'miRNA'=1, 'Other'=0))-1
  pred <- prediction(prediction, test.label)
  roc <- performance(pred, "tpr", "fpr")
  roc2 <- roc.curve(scores.class0 = pos_p, scores.class1 = neg_p, curve = T)
  AUROC <- roc2$auc
  roc <- data.frame(roc2$curve)
  pr <- pr.curve(scores.class0 = pos_p, scores.class1 = neg_p, curve = T)
  AUPR <- pr$auc.integral
  pr <- data.frame(pr$curve)
  
  Performance <- c(SE, SP, Acc, F1, AUROC, AUPR)
  names(Performance) <- c("Sensitivity", "Specificity", "Accuracy", "F-1 score", "AUROC", "AUPR")
  
  
  # Prediction
  
  message("Predicting novel miRNAs...")
  Table <- rbind(target, pos, neg, unlab)
  labels <- c(rep(0,nrow(target)), rep(1,nrow(pos)), rep(-1,nrow(neg)), rep(0,nrow(unlab)))
  
  p <- miRNAss(Table, labels, thresholdObjective = "F1")
  
  id_target <- rownames(Table[1:nrow(target),])
  p_target <- as.data.frame(p[1:nrow(target)])
  Class <- as.factor(ifelse(p_target > 0, "miRNA", "Other"))
  Odd_score <- rescale(p_target$`p[1:nrow(target)]`, to=c(0,1))
  
  Prediction_result <- data.frame(Odd_score, Class)
  rownames(Prediction_result) <- id_target
  
  
  setClass("eMIRNA",
           slots = list(Performance = "numeric",  Prediction = "data.frame",
                        PR = "data.frame", ROC = "data.frame"))
  
  results <- new("eMIRNA", Performance = Performance, Prediction = Prediction_result,
                 PR = data.frame(Precision=pr$curve[,1], Recall=pr$curve[,2]),
                 ROC = data.frame(TPR=roc$X1, FPR=roc$X2))
  
  return(results)
  
  
  } else if(is.null(unlab)==TRUE) {
    
    Table <- rbind(pos, neg)
    labels <- c(rep(1,nrow(pos)), rep(-1,nrow(neg)))
    
    set.seed(50)
    labels_test <- c(rep(1,nrow(pos)), rep(-1,nrow(neg)))
    labels_test[sample(which(labels_test > 0),nrow(pos)*0.25)] = 0
    labels_test[sample(which(labels_test < 0),nrow(neg)*0.25)] = 0
    
    p <- miRNAss(Table, labels_test)
    
    
    index_testing <- which(labels_test[1:(nrow(pos)+nrow(neg))] == 0)
    
    p_testing <- p[index_testing]
    
    npos <- length(which(index_testing <= nrow(pos)))
    nneg <- length(index_testing) - npos
    
    TP <- length(which(p_testing[1:npos] > 0)) / npos
    FN <- length(which(p_testing[1:npos] < 0)) / npos
    SE <- TP/(TP+FN)
    TN <- length(which(p_testing[(npos+1):((npos)+(nneg))] < 0)) / nneg
    FP <- length(which(p_testing[(npos+1):((npos)+(nneg))] > 0)) / nneg
    SP <- TN/(TN+FP)
    Acc <- (TP+TN)/(TP+TN+FP+FN)
    F1 <- (2*TP)/((2*TP)+FP+FN)
    
    
    final.table <- rbind(pos, neg)
    
    class <- c(rep("miRNA", nrow(pos)), rep("Other", nrow(neg)))
    
    final.table <- as.data.frame(cbind(final.table, class))
    final.table$class <- factor(final.table$class)
    
    testing <- final.table[index_testing,]
    
    p_testing <- rescale(p_testing, to=c(0,1))
    pos_p <- p_testing[1:(floor(nrow(pos)*0.25))]
    neg_p <- p_testing[(floor(nrow(pos)*0.25)+1):(floor(nrow(pos)*0.25)+floor(nrow(neg)*0.25))]
    prediction <- c(pos_p, neg_p)
    test.label <- as.integer(recode(testing$class, 'miRNA'=1, 'Other'=0))-1
    pred <- prediction(prediction, test.label)
    roc <- performance(pred, "tpr", "fpr")
    roc2 <- roc.curve(scores.class0 = pos_p, scores.class1 = neg_p, curve = T)
    AUROC <- roc2$auc
    roc <- data.frame(roc2$curve)
    pr <- pr.curve(scores.class0 = pos_p, scores.class1 = neg_p, curve = T)
    AUPR <- pr$auc.integral
    pr <- data.frame(pr$curve)
    
    Performance <- c(SE, SP, Acc, F1, AUROC, AUPR)
    names(Performance) <- c("Sensitivity", "Specificity", "Accuracy", "F-1 score", "AUROC", "AUPR")
    
    # Prediction
    
    message("Predicting novel miRNAs...")
    Table <- rbind(target, pos, neg)
    labels <- c(rep(0,nrow(target)), rep(1,nrow(pos)), rep(-1,nrow(neg)))
    
    p <- miRNAss(Table, labels, thresholdObjective = "F1")
    
    id_target <- rownames(Table[1:nrow(target),])
    p_target <- as.data.frame(p[1:nrow(target)])
    Class <- as.factor(ifelse(p_target > 0, "miRNA", "Other"))
    Odd_score <- rescale(p_target$`p[1:nrow(target)]`, to=c(0,1))
    
    Prediction_result <- data.frame(Odd_score, Class)
    rownames(Prediction_result) <- id_target
    
    
    setClass("eMIRNA",
             slots = list(Performance = "numeric",  Prediction = "data.frame",
                          PR = "data.frame", ROC = "data.frame"))
    
    results <- new("eMIRNA", Performance = Performance,
                   Prediction = Prediction_result,
                   PR = data.frame(Precision=pr$X1, Recall=pr$X2),
                   ROC = data.frame(TPR=roc$X2, FPR=roc$X1))
    
    return(results)
  
  } 
  
}


### eMIRNA.Network ###

eMIRNA.Network <- function(mRNA, miRNA, targets, cor=-0.5, type="pearson",
                           cpm=1, percent=0.5, normalize=FALSE){
  suppressMessages(require(PCIT))
  suppressMessages(require(NOISeq))
  suppressMessages(require(edgeR))
  suppressMessages(require(igraph))
  
  #Process data
  
  if(normalize == TRUE){
    
    index_mRNA <- rowSums(cpm(mRNA)>cpm) >= ncol(mRNA)*percent
    mRNA <- mRNA[index_mRNA,]
    index_miRNA <- rowSums(cpm(miRNA)>cpm) >= ncol(miRNA)*percent
    miRNA <- miRNA[index_miRNA,]
    mRNA <- log2(tmm(mRNA)+1)
    miRNA <- log2(tmm(miRNA)+1)
    
    #Filter by targets and format inputs
    index_mRNA <- which(row.names(mRNA) %in% targets[,1])
    mRNA_filtered <- mRNA[index_mRNA,]
    index_miRNA <- which(row.names(miRNA) %in% targets[,2])
    miRNA_filtered <- miRNA[index_miRNA,]
    mRNA <- t(mRNA_filtered)
    miRNA <- t(miRNA_filtered)
    exprs <- cbind(mRNA, miRNA)
    
  } else if(normalize == FALSE){
    
    index_mRNA <- which(row.names(mRNA) %in% targets[,1])
    mRNA_filtered <- mRNA[index_mRNA,]
    index_miRNA <- which(row.names(miRNA) %in% targets[,2])
    miRNA_filtered <- miRNA[index_miRNA,]
    mRNA <- t(mRNA_filtered)
    miRNA <- t(miRNA_filtered)
    exprs <- cbind(mRNA, miRNA)
    
  }
  
  #Calculate correlations and PCIT
  message("Calculating correlations...")
  cor_exprs <- cor(exprs, method=type)
  message("Filtering correlations...")
  pcit_exprs <- pcit(cor_exprs)
  signif_pcit <- idx(pcit_exprs)
  nonsignif_pcit <- idxInvert(nrow(cor_exprs), signif_pcit)
  cor_exprs[nonsignif_pcit] <- 0
  
  #Get Edge List from PCIT
  message("Creating edge table...")
  edge_pcit <- graph_from_adjacency_matrix(cor_exprs, mode="undirected",
                                           weighted="correlation")
  edge_pcit <- as_data_frame(edge_pcit, "edges")
  
  #Get mRNA-miRNA interactions and filter by threshold
  index1 <- which(edge_pcit$from %in% rownames(miRNA_filtered))
  index2 <- which(edge_pcit$to %in% rownames(miRNA_filtered))
  index3 <- index2[! index2 %in% index1]
  index4 <- index1 %in% index2
  index4 <- which(index4==FALSE)
  index4 <- index1[index4]
  index <- c(index3, index4)
  edges <- edge_pcit[index,]
  edges <- subset(edges, correlation <= cor)
  edges.i1 <- which(edges$from %in% rownames(mRNA_filtered))
  `%notin%` <- Negate(`%in%`)
  edges.i2 <- which(edges$from %notin% rownames(mRNA_filtered))
  edges1 <- edges[edges.i1,]
  edges2 <- edges[edges.i2,]
  edges.mRNA1 <- edges1$from
  edges.miRNA1 <- edges1$to
  edges.corr1 <- edges1$correlation
  edges1 <- data.frame(edges.miRNA1, edges.mRNA1, edges.corr1)
  colnames(edges1) <- c("from", "to", "correlation")
  edges <- rbind(edges1, edges2)
  colnames(edges) <- c("miRNA", "mRNA", "Correlation")
  
  return(edges)
  
}


### eMIRNA.RIF ###

eMIRNA.RIF <- function(mRNA, miRNA, design, network, DElist,
                                     cpm=1, percent=0.5, normalize=FALSE){
  suppressMessages(require(NOISeq))
  suppressMessages(require(edgeR))
  
  #Process data
  
  if(normalize == TRUE){
    
    index_mRNA <- rowSums(cpm(mRNA)>cpm) >= ncol(mRNA)*percent
    mRNA <- mRNA[index_mRNA,]
    index_miRNA <- rowSums(cpm(miRNA)>cpm) >= ncol(miRNA)*percent
    miRNA <- miRNA[index_miRNA,]
    mRNA <- log(tmm(mRNA)+1)
    miRNA <- log(tmm(miRNA)+1)
    exprs <- rbind(mRNA, miRNA)
    
  } else if(normalize == FALSE){
    
    exprs <- rbind(mRNA, miRNA)
    
  }
  
  design <- model.matrix(~0+design[,2])
  i1 <- as.vector(which(design[,1]==1))
  i2 <- as.vector(which(design[,2]==1))
  exprs.1 <- exprs[,i1]
  exprs.2 <- exprs[,i2]
  
  #Filter data
  network1 <- data.frame(unique(network[,1]))
  network2 <- data.frame(unique(network[,2]))
  
  reg.exprs.1 <- unique(merge(exprs.1,network1,by.x="row.names",by.y="miRNA"))
  rownames(reg.exprs.1)<-reg.exprs.1[,1]
  reg.exprs.1 <- reg.exprs.1[,2:ncol(reg.exprs.1)]
  
  reg.exprs.2 <- unique(merge(exprs.2,network1,by.x="row.names",by.y="miRNA"))
  rownames(reg.exprs.2)<-reg.exprs.2[,1]
  reg.exprs.2 <- reg.exprs.2[,2:ncol(reg.exprs.2)]
  
  targ.exprs <- unique(merge(exprs,network2,by.x="row.names",by.y="mRNA"))
  rownames(targ.exprs)<-targ.exprs[,1]
  targ.exprs <- targ.exprs[,2:ncol(targ.exprs)]
  
  RIF_reg_rank1 <- data.frame(miRNA=row.names(reg.exprs.1),RIF1=FALSE)
  RIF_reg_rank2 <- data.frame(miRNA=row.names(reg.exprs.1),RIF2=FALSE)
  rownames(RIF_reg_rank1)<-RIF_reg_rank1[,1]
  rownames(RIF_reg_rank2)<-RIF_reg_rank2[,1]
  
  #DE filtering
  targ.exprs <- targ.exprs[rownames(targ.exprs) %in% DElist[,1],]
  lengthDE <- nrow(targ.exprs)
  targ.exprs.1 <- targ.exprs[,1:n.1]
  targ.exprs.2 <- targ.exprs[,(n.1+1):(n.1+n.2)]
  
  #Mean Calc by group
  e1 <- as.data.frame(apply(targ.exprs.1,1,mean))
  e2 <- as.data.frame(apply(targ.exprs.2,1,mean))
  
  
  nreg.exprs<-nrow(reg.exprs.1)
  
  #RIF1
  for( i in 1:nreg.exprs){
    cor1 <- cor(t(targ.exprs.1),t(reg.exprs.1[rownames(RIF_reg_rank1)[i],]),
                method="pearson",use="pairwise.complete.obs")
    cor2 <- cor(t(targ.exprs.2),t(reg.exprs.2[rownames(RIF_reg_rank1)[i],]),
                method="pearson",use="pairwise.complete.obs")
    PIF <- (e1^2-e2^2)*0.5
    DW2 <- (cor1-cor2)^2
    preRIF <- PIF*DW2
    preRIF <- apply(preRIF,2,sum)
    preRIFdividlength <- preRIF/lengthDE
    RIF_reg_rank1[i,2] <- preRIFdividlength[[names(preRIFdividlength)]]
    RIF1 <- RIF_reg_rank1
  }
  
  RIF1 <- RIF1[order(-abs(as.numeric(RIF1[,'RIF1']))),]
  RIF1_final <- as.data.frame(scale(RIF1$RIF1))
  RIF1_final <- data.table(RIF1$miRNA,RIF1_final$V1)
  colnames(RIF1_final) <- c("miRNA","RIF1")
  rownames(RIF1_final) <- rownames(RIF1)
  
  #RIF2
  for( i in 1:nreg.exprs){
    cor1 <- cor(t(targ.exprs.1),t(reg.exprs.1[rownames(RIF_reg_rank2)[i],]),
                method="pearson",use="pairwise.complete.obs")
    cor2 <- cor(t(targ.exprs.2),t(reg.exprs.2[rownames(RIF_reg_rank2)[i],]),
                method="pearson",use="pairwise.complete.obs")
    preRIF <- (e1*cor1)^2-(e2*cor2)^2
    preRIF <- apply(preRIF,2,sum)
    preRIFdividlength <- preRIF/lengthDE
    RIF_reg_rank2[i,2] <- preRIFdividlength[[names(preRIFdividlength)]]
    RIF2 <- RIF_reg_rank2
  }
  
  RIF2 <- RIF2[order(-abs(as.numeric(RIF2[,'RIF2']))),]
  RIF2_final <- as.data.frame(scale(RIF2$RIF2))
  RIF2_final <- data.table(RIF2$miRNA,RIF2_final$V1)
  colnames(RIF2_final) <- c("miRNA","RIF2")
  rownames(RIF2_final) <- rownames(RIF2)
  
  RIF <- merge(RIF1_final,RIF2_final,by=row.names(RIF1_final$miRNA))
  RIF <- RIF[order(-RIF1),]
  
  return(RIF)
  
}
