library(DESeq2)
library(tximport)
library(rjson)
library(readr)
library(ashr)
library(ggplot2)
source("change_PC_plotPCA.R")
samples <- read.table("samples.txt", header=TRUE)
files <- file.path("github/results",paste0(samples$Sample,"_quant.sf"))
# txOut argument avoid genelevelsummary since gene level information is not available.
txi <- tximport(files, type="salmon",txOut=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ group)

## Differential Expression Function
dds <- DESeq(ddsTxi)

## Filter Reads with less than 10 reads across 12 samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Differential Expression Function
dds <- DESeq(dds)

## Including MultiFactor
ddsMF <- dds
design(ddsMF) <- formula(~ side + group:side + group, colData(ddsMF))

## Differential Expression Function for MultiFactor
ddsMF <- DESeq(ddsMF)

resultsNames(ddsMF)
# getting the results from the DE seq function for MF design
resMF <- results(ddsMF, contrast = c("group", "Treatment","Control"))
resMF
resMFOrdered <- resMF[order(resMF$pvalue),]

#write the results to a csv
write.csv(as.data.frame(resMFOrdered), 
          file="case_control.csv")
# Data Transformation
vsd <- vst(ddsMF, blind=FALSE)

# getting the results from the DE seq function
head(assay(vsd))
# pca plot
plotPCA(vsd, intgroup="group")

# pca plot removing batch effects due to left and right side of the same mouse
vsd <- vst(ddsMF, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat,vsd$side,vsd$MouseId)
assay(vsd) <- mat
plotPCA(vsd, intgroup="group")

plotPCA.balan(vsd, intgroup="group")
