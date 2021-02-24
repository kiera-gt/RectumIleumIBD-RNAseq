library("DESeq2")
library(SummarizedExperiment)
library(pvca)
library(lme4)
genecounts <- as.matrix(read.table("~/genecounts.csv", sep=",", header=TRUE, row.names=1))
mycoldata <- read.csv("~/sampleinfo.csv", sep=",", row.names = 1, header=TRUE)
mydds = DESeqDataSetFromMatrix(countData = genecounts, colData = sampleinfo, design = ~ Ancestry + Disease + Location + Disease:Location)
mydds <- mydds[ rowMeans(counts(mydds)) > 5, ]
vsd <- vst(mydds, blind=TRUE)
mydds <- DESeq(mydds)
###PCA/PVCA###
gene.pca <- prcomp(t(assay(vsd)))
gene.pca.proportionvariances <- ((gene.pca$sdev^2) / (sum(gene.pca$sdev^2)))*100
pct_threshold <- 0.75
batch.factors <- c("Ancestry", "Location", "Disease")
geneeset <- ExpressionSet(assay(vsd),as(sampleinfo,"AnnotatedDataFrame"))
pvcaobj <- pvcaBatchAssess(geneeset, batch.factors, pct_threshold)
pvca.res <- data.frame(label=as.character(pvcaobj$label),wmpv=round(as.numeric(pvcaobj$dat),3))
immgenes <- read.table("~/immunegenes.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names = 1)
immvsd <- vsd[ rownames(vsd) %in% rownames(immgenes), ]
immgenespc <- plotPCA(immvsd, intgroup=c("Location"), returnData=TRUE)
epigenes <- read.table("~/epithelialgenes.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names = 1)
epivsd <- vsd[ rownames(vsd) %in% rownames(epigenes), ]
epigenespc <- plotPCA(epivsd, intgroup=c("Location"), returnData=TRUE)
immpca.sig <- t.test(immgenespc$PC1~immgenespc$Location)
epipca.sig <- t.test(epigenespc$PC1~epigenespc$Location) 