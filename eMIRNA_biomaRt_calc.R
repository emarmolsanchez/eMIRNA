# Biomart miRNA Neighbourhood query

library(biomaRt)
args <- commandArgs(trailingOnly = TRUE)


# data
file <- paste0(args[1], args[2], "_temp7.4.bed")
d <- read.table(file)
data1 <- paste0(args[3], "_gene_ensembl")
data2 <- paste0(args[4], "_gene_ensembl")
data3 <- paste0(args[3], "_homolog_associated_gene_name")

# specify the database
ensembl = useMart("ensembl", dataset = data1)

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.
res <- cbind(
  d,
  genes = apply(d, 1, function(i){
    x <- getBM(attributes=c("external_gene_name"), 
               filters = c("chromosome_name" , "start", "end"), 
               values = list(i[1], i[2], i[3]), 
               mart = ensembl)
    
    # return genes, comma separated
    paste(x$external_gene_name, collapse = ",")
  })
)

genes <- as.vector(res$genes)
miRNAs <- as.vector(res$V4)
Results <- as.data.frame(cbind(miRNAs, genes))


# data
file2 <- paste0(args[1], args[2], "_temp14.2.bed")
d2 <- read.table(file2)

# specify the database
ensembl = useMart("ensembl", dataset = data2)

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.
res2 <- cbind(
  d2,
  genes = apply(d2, 1, function(i){
    x <- getBM(attributes=c(data3), 
               filters = c("chromosome_name" , "start", "end"), 
               values = list(i[1], i[2], i[3]), 
               mart = ensembl)
    
    # return genes, comma separated
    paste(x$sscrofa_homolog_associated_gene_name, collapse = ",")
  })
)

genes2 <- as.vector(res2$genes)
miRNAs2 <- as.vector(res2$V4)
Results2 <- as.data.frame(cbind(miRNAs2, genes2))



out <- paste0(args[1], args[2], "_temp15.txt")
out2 <- paste0(args[1], args[2], "_temp16.txt")

write.table(Results, out, quote=F, sep="\t", row.names=F)
write.table(Results2, out2, quote=F, sep="\t", row.names=F)