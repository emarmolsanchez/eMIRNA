#########################################
#                                       #
#          Scripts eMIRNA               #
#                                       #
#             EMS. 2018                 #
#                                       #
#########################################



eMIRNA.Filter.by.Size <- function(file, prefix, a, b){
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "FilterSize_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/FilterSize_Results/"
  setwd(workdir)
  
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
  name <- paste0(prefix, "_filter_size.fa")
  write(File0.filter.size, name)
  
  unlink("*.ps", recursive=T)
  
}
                             

  
eMIRNA.Filter.by.Structure <- function(file, prefix){
  suppressMessages(require(Biobase))
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "FilterStructure_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/FilterStructure_Results/"
  setwd(workdir)
  
  #Checking sequences and filtering by n loops
  message("Filtering sequences by Secondary Structure")
  File0 <- unlist(lapply(file, readLines))
  n0 <- length(File0)
  ID0 <- File0[seq(1,n0,2)]
  Sequence0 <- File0[seq(2,n0,2)]
  command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", file)
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
  name <- paste0(prefix, "_filter_nloop.fa")
  write(File0.filter.nloop, name)
  
  unlink("*.ps", recursive=T)
  
}



eMIRNA.Features <- function(file, prefix, rescale=TRUE){
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
  message("Calculating Features")
  if (substr(file, nchar(file)-16+1, nchar(file)) == '_filter_nloop.fa'){
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
  
  #Number of A-U, C-G and G-U pairs in 1 loop sequences (AUp, GCp, GUp)
  Sequence2 <- strsplit(Sequence, "")
  SecondaryStrc2 <- strsplit(SecondaryStrc, "")
  Matches5 <- sapply(SecondaryStrc2, function(x) which(x == "("))
  Matches3 <- sapply(SecondaryStrc2, function(x) which(x == ")"))
  
  if (class(Matches5) == "matrix"){
    
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
  EFreq <- round(as.numeric(regmatches(EFreq2,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          EFreq2, perl=TRUE))), 4)
  
  
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
  
  
    
 table <- as.data.frame(cbind(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14,
                                 T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25,
                                 T26, T27, T28, T29, T30, T31, T32, GC, EFreq, BPP, Dadj, SEadj,
                                 DiffMFE.EFE, Length, AU.GCr,
                                 Ar, Ur, Gr, Cr, AAr, GGr, CCr, UUr, AGr, ACr, AUr,
                                 GAr, GCr, GUr, CAr, CGr, CUr, UAr, UGr, UCr, Hl, Steml5,
                                 Steml3, BP, BP5, BP3, Mism5, Mism3, Bulge5, Bulge3, BulgeN1.5, 
                                 BulgeN1.3, BulgeN2.5, BulgeN2.3, BulgeN3.5, BulgeN3.3,
                                 BulgeN4.5, BulgeN4.3, BulgeN5.5, BulgeN5.3, AUp, GCp, GUp, MFE,
                                 MFEadj, EFE, EFEadj, CFE, CDE, MEAFE, MEA,
                                 ED, MFEadj.GC, MFEadj.BP))
    
    unlink("*.ps")
    unlink("*.txt")
    unlink("*.fa*")
    
    table_1 <- table[, 1:32]
    table_2 <- table[, 33:94]
    
    if (rescale == TRUE){
      table_2 <- as.data.frame(lapply(table_2, rescale))
      table <- cbind(table_1, table_2)
      rownames(table) <- ID
      table[is.na(table)] <- 0
    
    } else {
      table <- cbind(table_1, table_2)
      rownames(table) <- ID
      table[is.na(table)] <- 0
      
      
    }
    
    final_table <- paste0(workdir, prefix, ".csv")
    write.table(table, final_table, sep=",", quote=F, col.names=NA)
    
    return(table)
       
                   
  } else {
    
    stop("Input file is not a nloop filtered FASTA file.")
    
    
 }
    
  
}


eMIRNA.Train <- function(pos, neg, imbalance="none"){
  suppressMessages(require(caret))
  suppressMessages(require(bimba))
  suppressMessages(require(LiblineaR))
  message("Training SVM model")
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "SVM_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/SVM_Results/"
  setwd(workdir)
  
  final.table <- rbind(pos, neg)
  
  class <- c(rep("miRNA", nrow(pos)), rep("Other", nrow(neg)))
  
  final.table <- as.data.frame(cbind(final.table, class))
  final.table$class <- factor(final.table$class)
  
  if(imbalance == "smote"){
    
  set.seed(1234)
  algorithms <- c("NRAS", "SMOTE")
  final.table <- sampling_sequence(final.table, algorithms=algorithms)
  intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
  training <- final.table[intrain,]
  testing <- final.table[-intrain,]
  
  
  #SVM Linear
  grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
  
  trctrl <- trainControl(method = "cv", number = 10,
                         savePredictions = TRUE, classProbs=TRUE)
  
  svm_Linear <- train(class ~., data = training, method = "svmLinear",
                      trControl = trctrl,
                      tuneGrid = grid_lineal,
                      tuneLength = 10)
  
  
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
  
  }  else if(imbalance == "adasyn"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "ADASYN")
    final.table <- sampling_sequence(final.table, algorithms=algorithms)
    
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "bdlsmote1"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "BDLSMOTE")
    list1 <- list()
    list2 <- list(borderline=1)
    parameters <- list(list1, list2)
    final.table <- sampling_sequence(final.table, algorithms=algorithms,
                                     parameters=parameters)
    
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "bdlsmote2"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "BDLSMOTE")
    list1 <- list()
    list2 <- list(borderline=2)
    parameters <- list(list1, list2)
    final.table <- sampling_sequence(final.table, algorithms=algorithms,
                                     parameters=parameters)
    
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "mwmote"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "MWMOTE")
    final.table <- sampling_sequence(final.table, algorithms=algorithms)
  
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "ros"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "ROS")
    final.table <- sampling_sequence(final.table, algorithms=algorithms)
    
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "rwo"){
   
    set.seed(1234)
    algorithms <- c("NRAS", "RWO")
    final.table <- sampling_sequence(final.table, algorithms=algorithms)

    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
   
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "slsmote"){
    
    set.seed(1234)
    algorithms <- c("NRAS", "SLSMOTE")
    final.table <- sampling_sequence(final.table, algorithms=algorithms)

    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else if(imbalance == "none"){
    
    set.seed(1234)
    intrain <- createDataPartition(y = final.table$class, p= 0.8, list = FALSE)
    training <- final.table[intrain,]
    testing <- final.table[-intrain,]
    
    
    #SVM Linear
    grid_lineal <- expand.grid(C = seq(0.01, 1, 0.1))
    
    trctrl <- trainControl(method = "cv", number = 10,
                           savePredictions = TRUE, classProbs=TRUE)
    
    svm_Linear <- train(class ~., data = training, method = "svmLinear",
                        trControl = trctrl,
                        tuneGrid = grid_lineal,
                        tuneLength = 10)
    
    
    write.table(training, "~/eMIRNA/SVM_Results/training.csv", sep=",", quote=F, col.names=NA)
    write.table(testing, "~/eMIRNA/SVM_Results/testing.csv", sep=",", quote=F, col.names=NA)
    saveRDS(svm_Linear, "~/eMIRNA/SVM_Results/SVM.rds")
    return(svm_Linear)
    
  } else {
    
    message("SVM Training Failed. Please provide a correct algorithm for class imbalance resolution.")
  }
  
}




eMIRNA.Predict <- function(model, features, prefix){
  pred <- predict(model, newdata= features)
  
  index <- which(pred == "miRNA")
  
  pred <- as.data.frame(rownames(features[index,]))
  colnames(pred) <- "Predicted_miRNAs"
  
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "Prediction_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/Prediction_Results/"
  setwd(workdir)
  
  pred.path <- paste0(workdir, prefix, ".txt")
  write.table(pred, pred.path, quote=F, col.names=F, row.names=F)
  
  return(pred)
  
}


eMIRNA.Structural.Pvalues <- function(file, prefix, iterate=100){
  setwd("~/")
  Dir0 <- "eMIRNA"
  dir.create(file.path(Dir0), showWarnings=FALSE)
  Dir <- "Prediction_Results"
  setwd("~/eMIRNA/")
  dir.create(file.path(Dir), showWarnings = FALSE)
  workdir <- "~/eMIRNA/Feature_Results/"
  setwd(workdir)
  message("Calculating P-values")
  command1 <- paste0("RNAfold --MEA -d2 -p --noPS -i ", file)
  RNAfold1 <- system(command1, intern=TRUE)
  
  
  command2 <- paste0("RNAfold --noPS -i ", file, " > RNAfold_pred.txt")
  system(command2)
  command2.1 <- "1_check_query_content.pl RNAfold_pred.txt RNAfold_pred_checked.txt"
  system(command2.1)
  command2.2 <- "2_get_stemloop.pl RNAfold_pred_checked.txt RNAfold_pred_stemloop.txt 5"
  system(command2.2)
  
  
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
  
  
  #RNAfold Secondary Structure (SecondaryStrc)
  SecondaryStrc1 <- RNAfold1[seq(3,n,7)]
  SecondaryStrc1 <- SecondaryStrc1[StemF.ID.index]
  SecondaryStrc <- gsub( " .*$", "", SecondaryStrc1)
  
  
  #Minimum Freen Energy (MFE)
  MFE <- as.numeric(regmatches(SecondaryStrc1,
                               gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                        SecondaryStrc1, perl=TRUE)))
  
  #Pairing Probabilities Structure (PairProbStrc)
  PairProbStrc1 <- RNAfold1[seq(4,n,7)]
  PairProbStrc1 <- PairProbStrc1[StemF.ID.index]
  PairProbStrc <- gsub( " .*$", "", PairProbStrc1)
  
  #Ensemble Free Energy (EFE)
  EFE <- as.numeric(regmatches(PairProbStrc1,
                               gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                        PairProbStrc1, perl=TRUE)))
  
  
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
  
  
  #EFE P-value (EFE.Pval)
  EFEit1 <- RNAfold2[seq(4,n4,7)]
  EFEit <- as.numeric(regmatches(EFEit1,
                                 gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",
                                          EFEit1, perl=TRUE)))
  
  EFEit.along <- seq_along(EFEit)
  EFEit.split <- split(EFEit, ceiling(EFEit.along/N))
  R2 <- as.vector(sapply(Map(`<=`, EFEit.split, EFE), sum))
  EFE.Pval <- (R2 / (N + 1))
  
  unlink("*.ps")
  unlink("*.txt")
  unlink("*.fa")
  
  MFE.Pval <- round(MFE.Pval, 6)
  EFE.Pval <- round(EFE.Pval, 6)
  table <- as.data.frame(cbind(MFE.Pval, EFE.Pval))
  rownames(table) <- ID
  table <- table[order(table$MFE.Pval),]
  
  final_table <- paste0("~/eMIRNA/Prediction_Results/", prefix, "_Structural_Pvalues.csv")
  write.table(table, final_table, sep=",", quote=F, col.names=NA)
  
  return(table)
  
}

