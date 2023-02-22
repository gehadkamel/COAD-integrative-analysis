library(BiocManager)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(NMF)

#loading the data
mirna_data <- read.csv("mirna", fill= T, header = T, sep = "")
survival_data <- read.csv("survival", fill = T, header = T, sep= "")
clinical <- read.delim("colon", fill = T, header = T)
#adjust sample names
clinical[,1] <-gsub("-", ".", clinical[,1])
rownames(clinical) <- clinical[,1]
survival_data[,1] <-gsub("-", ".", survival_data[,1])
rownames(survival_data) <- survival_data[,1]
#check for empty cells
sum(is.null(mirna_data))
sum(is.na(mirna_data))
sum(is.nan(as.matrix(mirna_data)))
#Matching the data
all(rownames(clinical) == colnames(mirna_data))

survival_data <- survival_data[which(rownames(survival_data) %in% rownames(clinical)),]
survival_data <- survival_data[match(rownames(survival_data),rownames(clinical)),]


clinical <- clinical[which(rownames(clinical) %in% colnames(mirna_data)),]
clinical <- clinical[match(rownames(clinical), colnames(mirna_data)),]

all(rownames(clinical) %in% colnames(mirna_data))
all(rownames(clinical) == colnames(mirna_data))

survival_data <- survival_data[which(rownames(survival_data) %in% colnames(mirna_data)),]
survival_data <- survival_data[match(rownames(survival_data), colnames(mirna_data)),]

all(rownames(survival_data) %in% colnames(mirna_data))
all(rownames(survival_data) == colnames(mirna_data))
#Histogram mirna
hist(as.matrix(mirna_data))
hist(log2(as.matrix(mirna_data)))
#boxplot
boxplot(log2(mirna_data[1:5,]+1))
#QQ plot for the normality
qqnorm(mirna_data[30,])
qqline(mirna_data[30,])
#save the gene names in a variable
genes_mirna <- rownames(mirna_data)
#convert the data values to integers
mirna_data <- apply(round(mirna_data),2 , as.integer)
#view the data
head(mirna_data)
#rename the rows of the data
rownames(mirna_data) <- genes_mirna
#DESEQ2 Differential exepression

#create a deseq dataset object using clinical file
cond1="LIVING" 
cond2="DECEASED"
dds_mirna <- DESeqDataSetFromMatrix(countData = mirna_data, colData = clinical,design = ~vital_status)
dds_mirna
keep <- rowSums(counts(dds_mirna)) >= 10
dds_mirna <- dds_mirna[keep,]
dds_mirna <- estimateSizeFactors(dds_mirna)
#estimateSizeFactorForMatrix
sizeFactors(dds_mirna)

normalized_mirna <- counts(dds_mirna, normalized = TRUE)
qqnorm(normalized_mirna[50,])
qqline(normalized_mirna[50,])

# working on Survival file only to make sure we're getting the same results
dds_mirna2 <- DESeqDataSetFromMatrix(countData = mirna_data, colData = survival_data,design = ~PatientID)
dds_mirna2
keep <- rowSums(counts(dds_mirna2)) >= 10
dds_mirna2 <- dds_mirna2[keep,]
dds_mirna2 <- estimateSizeFactors(dds_mirna2)
sizeFactors(dds_mirna2)

normalized_mirna2 <- counts(dds_mirna2, normalized = TRUE)
qqnorm(normalized_mirna2[50,])
qqline(normalized_mirna2[50,])

#DESeq2 workflow
cond1= "DECEASED"
cond2="LIVING"
library(DESeq2)
dds_mirna= DESeqDataSetFromMatrix( countData = mirna_data , colData = clinical, design = ~vital_status)
dds.run = DESeq(dds_mirna)
#specifying the contrast(to make a res object based on two specific conditions)
res=results(dds.run, contrast = c( "vital_status" ,cond1 ,cond2))
#remove nulls
res=as.data.frame(res[complete.cases(res), ])
#chose the statistical significant differentially expressed genes (DEGS) based
#on the p adjusted value less than 0.05 and biological significance based on Fold change >0.5
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>0.5,]
deseq.deg2=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.2,]
#export the DEGS for further analysis
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)
#draw DEGS volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs ALIVE DEGs"))
with(subset(res, padj<.05 & (log2FoldChange)>0.5), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -0.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-3.7,y=2.5,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)
dds2 <- estimateSizeFactors(dds_mirna)

#Heatmap
#normalize the data
normalized_mirna <- as.data.frame(counts(dds2, normalized=TRUE))
#extract counts values of DEGS only for each stage
exp.degs=as.matrix(normalized_mirna[rownames(normalized_mirna) %in% rownames(deseq.deg), ])
heatmap(log2(exp.degs+1), annCol =clinical$vital_status, col = rev(brewer.pal(9,"RdBu")), main="DECEASED vs ALIVE")

heatmap(log2(exp.degs+1))
