
library(BiocManager)
BiocManager::install("DESeq2",force = TRUE)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(NMF)

#loading the gene expression data 
data2 = as.matrix(read.table("exp",header = T,row.names = 1))

#loading the phenotype data in case that data isn't in excel form by read.delim
pheno2 = read.delim("colon",row.names = 1)


#explore the data, dim function give you the dimension of your data; how 
#many columns(samples) do you have and how many row(genes) do you have
dim(data2)
dim(pheno2)


#explore the data distribution using the histogram plot
hist(data2, col = "orange", main="Histogram")

#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(data2+1) ,col = "orange",main = "Histogram",breaks = 100)

# QQ plot for the normality
qqnorm(data2[1,])
qqline(data2[1])

#explore if is there any missing expression value (empty cell)
sum(is.na(data))
is.na(data) 
is.null(data) 
is.nan(data) 

#preparing data
#for rownames if it wasn't selected from the begenning
row.names(pheno)<- pheno[,1]
#for naming
rownames(pheno2)<- gsub("-",".",rownames(pheno2))

#It is absolutely critical that the columns of the "data" and the rows of 
#the "pheno" (information about samples) are in the same order.DESeq2 will
#not make guesses as to which column of the count matrix belongsto which row
#of the column data, these must be provided to DESeq2 already in consistent order

#to find the intersection between them
pheno2 <- pheno2[which(rownames(pheno2) %in% colnames(data2)), ]
#match and check if they are matching
pheno2 <- pheno2[match(colnames(data2), rownames(pheno2)), ]

#to find common between them
pheno2=pheno2[colnames(data2),]


#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes2= row.names(data2)

#convert the data values to integers
data2=apply(round(data2),2, as.integer)

#view the data
head(data2)
#rename the rows of the data
row.names(data2)=genes2
#view the data

###### DO the differential EXP analysis using DeSeq2
#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="DECEASED"
cond2="LIVING"
#for modeling
model.matrix(~0+pheno2$vital_status)

#creat a deseq dataset object from matrix
dds2= DESeqDataSetFromMatrix( countData = data2 , colData = pheno2, design = ~ vital_status)
#to make a model for storing pipeline output
dds.run2 = DESeq(dds2)
#results to be compared
res2=results(dds.run2, contrast = c("vital_status",cond1 ,cond2))
#to make a data frame and remove nulls and zeros
res2=as.data.frame(res2[complete.cases(res2), ])

boxplot(res2) 

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
deseq.deg2=res2[res2$padj < 0.05 & abs(res2$log2FoldChange)> 1.2,]
#to document in normal file
write.csv(as.matrix(deseq.deg2),file="deseq.deg2.csv", quote=F,row.names=T)

#using gprofiler for enrichment analysis
#gland development,complex,extracellelular spaces

#to visualize degs through volcanoplot
par(mfrow=c(1,1))
with(res2, plot(log2FoldChange, -log10(padj), pch=20, main=" living vs decease DEGs"))
with(subset(res2, padj<.05 & (log2FoldChange)>1.2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res2, padj<.05 & (log2FoldChange)< -1.2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-1,y=15,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

#visualization via heatmap
dds3 <- estimateSizeFactors(dds2)
normalized_counts2 <- as.data.frame(counts(dds3, normalized=TRUE))

# to find out common 
exp.degs2=as.matrix(normalized_counts2[rownames(normalized_counts2) %in% rownames(deseq.deg2), ])
heatmap(log2(exp.degs2+1), annCol =pheno2$vital_status, col = rev(brewer.pal(9,"RdBu")), main="mRNA living vs deceased")
heatmap(log2(exp.degs2+1))

# Load libraries
library(tidyverse)
library(princurve)

# Perform PCA 
pca <- prcomp(data, scale = TRUE) 

# Plot the results 
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", main="PCA of RNA seq data")

