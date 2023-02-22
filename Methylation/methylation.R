install.packages("knitr")
BiocManager::install("limma")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("missMethyl")
BiocManager::install("DMRcate")
BiocManager::install("DMRcatedata")
BiocManager::install("minfiData")
BiocManager::install("methylationArrayAnalysis")

BiocManager::install("TCGAbiolinks")
#Load library
# load packages required for analysis
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(gplots)
library(Gviz)
library(DMRcate)
library(stringr)
library(DMRcatedata)
library(TCGAbiolinks)
################################
##1. Loading Data
#Reading methylation raw data
methyl_raw = read.csv("Data/methy", fill = T, header = T, sep="")
clinical_data <- read.delim("Data/clinical/colon", fill = T, header = T)
#Adjust sample names to the sample names in methylation data
clinical_data[,1] <- gsub("-", ".", clinical_data[,1])
samples <- c(clinical_data[,1] , count = nrow(clinical_data))

rownames(clinical_data) <- clinical_data[,1]
clinical_data_f <- clinical_data[which(rownames(clinical_data)%in% colnames(methyl_raw)),] 
#match sample names with sample names in raw data
clinical_data_f <- clinical_data_f[match(colnames(methyl_raw),rownames(clinical_data_f)),]
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")

#Survival data
survival <- read.delim("Data/survival", fill = T, header = T)

#Get 450K annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

#Since data is preprocessed already so QC and normalization wouldnt be done in our pipeline
####################################

#2.Exploring Data distribution
#Since the data is in Beta values; there is no good way to normalize the ratios as all reasonable
#nomralization methods take into account the two different types of probes
#https://support.bioconductor.org/p/83204/


#Make sure no zero value or NAs in data
sum(is.na(methyl_raw))
sum(is.null(methyl_raw))

#Survival Data
TCGAanalyze_survival(
  data = clinical_data,
  clusterCol = "gender",
  main = "TCGA Set\n GBM",
  height = 10,
  width=10
)
design <- summary(factor(clinical_data_f$sample_type))

#hist(design, ylab="frequency")
#convert matrix into genomic data set so that I can apply minfi functions
#Source:https://rdrr.io/bioc/minfi/man/makeGenomicRatioSetFromMatrix.html
grset <- makeGenomicRatioSetFromMatrix(as.matrix(methyl_raw))
####################################3
##3.Filtering
#Exclude X and Y chromosomes
keep <- !(featureNames(grset) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
grset <- grset[keep,]
#Filter probes with known SNPs
grset <- dropLociWithSnps(grset)
#Filter out cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)

keep_c <- !(featureNames(grset) %in% xReactiveProbes$TargetID)
table(keep_c)
grset <- grset[keep_c,]

jpeg(file = "densityplot.jpeg")
par=c(1,2)
densityPlot(as.matrix(getBeta(grset)), sampGroups=clinical_data_f$sample_type, main = "Raw", legend = FALSE)
legend("top", legend = levels(factor(clinical_data$sample_type)), text.col=brewer.pal(8,"Dark2"))

densityPlot(as.matrix(getBeta(grset)), sampGroups=clinical_data_f$gender, main ="Beta", legend = FALSE)

legend("top", legend=levels(factor(clinical_data_f$gender)), main = "Beta", text.col=brewer.pal(8,"Dark2"))


dev.off()

jpeg(file="Histogram.jpeg")
hist(as.matrix(methyl_beta),sampGroups=clinical_data_f$sample_type, legend = FALSE, xlab="Beta")
legend("top", legend=levels(factor(clinical_data_f$sample_type)), main = "Beta", text.col = brewer.pal(8, "Dark2"))
dev.off()


#Heat map


heatmap.2(as.matrix(methyl_beta), scale = "none", col= bluered(100), trace = "none"
          ,denisty.info = "none")



#3.Data Exploration


#Using MDS plots to explore source of variability

pal <- brewer.pal(8,"Dark2")
#Perform MDS analysis suspecting sample type (Tumor Vs Normal)
plotMDS(getM(grset), top=1000, gene.selection="common",
        col=pal[factor(clinical_data_f$sample_type)], pch = 20)
legend("top", legend=levels(factor(clinical_data_f$sample_type)), text.col=pal,
       bg="white", cex=0.7)

#Plotting MDS suspecting vital status
plotMDS(getM(grset), top=1000, gene.selection="common",
        col=pal[factor(clinical_data_f$vital_status)], pch = 20)
legend("top", legend=levels(factor(clinical_data_f$vital_status)), text.col=pal,
       bg="white", cex=0.7)


#Plotting MDS suspecting gender discrimination
pal2 <- brewer.pal(12,"Set3")
plotMDS(getM(grset), top=1000, gene.selection="common",
        col=pal[factor(clinical_data_f$gender)],pch = 20)
legend("top", legend=levels(factor(clinical_data_f$gender)), text.col=pal,
       bg="white", cex=0.7)

#Plotting higher dimension 

par(mfrow=c(1,3))
plotMDS(getM(grset), top=1000, gene.selection="common", 
        col=pal[factor(clinical_data_f$sample_type)], dim=c(1,3), pch = 20)
legend("top", legend=levels(factor(clinical_data_f$sample_type)), text.col=pal, 
       cex=0.7, bg="white")



plotMDS(getM(grset), top=1000, gene.selection="common", 
        col=pal[factor(clinical_data_f$sample_type)], dim=c(2,3), pch = 20)
legend("topleft", legend=levels(factor(clinical_data_f$sample_type)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(grset), top=1000, gene.selection="common", 
        col=pal[factor(clinical_data_f$sample_type)], dim=c(3,4), pch=20)
legend("topright", legend=levels(factor(clinical_data_f$sample_type)), text.col=pal,
       cex=0.7, bg="white")


d <- dropLociWithSnps(grset)
###########################################################

bVals <- getBeta(grset)
head(bVals)
mVals <- getM(grset)

densityPlot(bVals, sampGroups = clinical_data_f$sample_type, main = "Beta values", legend = FALSE)

densityPlot(mVals, sampGroups = clinical_data_f$sample_type, main= "M-values", legend = FALSE )


#####################################################################################

#4. Probe-wise differential methylation analysis

#filter individuals with only Primary Tumor and Solid tissue normal
clinical_data_c <- clinical_data_f[clinical_data_f$sample_type == "Solid Tissue Normal" | clinical_data_f$sample_type == "Primary Tumor", ]
bVals <- bVals[, which(colnames(bVals) %in% rownames(clinical_data_c))]

sample_type <- factor(clinical_data_c$sample_type)
levels(sample_type)[levels(sample_type) == "Solid Tissue Normal"] <- "Solid_Tissue_Normal"
levels(sample_type)[levels(sample_type) == "Primary Tumor"] <- "Primary_Tumor"

clinical_data_c$sample_type[clinical_data_c$sample_type == "Solid Tissue Normal"] <- "Solid_Tissue_Normal"
clinical_data_c$sample_type[clinical_data_c$sample_type == "Primary Tumor"] <- "Primary_Tumor"

#Same for bVals 
bVals <- bVals[, which(colnames(bVals) %in% rownames(clinical_data_c))]
bVals <- bVals[, match(colnames(bVals), rownames(clinical_data_c))]
#Get M-values for statistical analysis
mVals <- getM(grset)
mVals <- mVals[,which(rownames(clinical_data_c)%in% colnames(mVals))]
#mVals <- mVals[,match(colnames(bVals),rownames(clinical_data_c))]

design <- model.matrix(~0+sample_type, data=clinical_data_c)
colnames(design) <- c("Primary_Tumor", "Solid_Tissue_Normal")
#Fit linear model
fit <- lmFit(mVals, design)

contMatrix <- makeContrasts(Primary_Tumor-Solid_Tissue_Normal,
                            levels=design)

#Fit the contrasts

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

#Plot pie chart for the  summary of the fitted model
piepercent <- round(100*summary(decideTests(fit2)) / sum(summary(decideTests(fit2))), 1)
summary_m <- summary(decideTests(fit2))
png(file="Methylation profile percentage.jpg")

pie(summary_m, labels= piepercent, main="Methylation levels",col=rainbow(length(summary_m)))
legend("topright", c("Down regulated", "Not significant","Up regulated"),cex=0.8, 
       fill = rainbow(length(summary_m)))

dev.off()

#get the table of result of the contrast( Ordered by B-statistic which is log(odds ratio that 
#the probe is differentially methylated))

ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
#Using topTable 
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
#write csv file
write.table("DMPs", file="DMPs.csv", sep=",", row.names = FALSE)

#Plot top differentially methylated Probes
par(mfrow=c(2,2))

sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=clinical_data_c$sample_type, ylab = "Beta values")
})

#######################

#Identify differentially methylated regions based on dmrcate of limma package

#Annotate regions based on genomic position and gene name
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "Primary_Tumor - Solid_Tissue_Normal", arraytype = "450K")
#Identify differentially methylated regions

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
#Visualize the regions using dmr plot

groups <- pal[1:length(unique(clinical_data_c$sample_type))]
names(groups) <- levels(factor(clinical_data_c$sample_type))
cols <- groups[as.character(factor(clinical_data_c$sample_type))]

tiff("plot.tiff", width=10, height=12,res = 200, units="in")
par(mfrow=c(1,1))

DMR.plot(ranges= results.ranges, dmr = 1, CpGs = bVals, phen.col=cols
         , what = "Beta", arraytype="450", genome="hg19")
dev.off()


#Enrichment analysis

#Get the most significant CpGs
sigDMRs  <- results.ranges[results.ranges$HMFDR < 0.05]
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
length(sigCpGs)

#Explore top 10 probes
sigCpGs[1:10]
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
length(all)
par(mfrow=c(1,1))


#Kegg pathways invloved with CpGs
gst <- gometh(sig.cpg=sigCpGs, collection = 'KEGG' ,plot.bias=TRUE)

topGSA(gst, number = 20)
write.csv(topGSA(gst, number=20), "Top 5 pathways.csv")
gst.region.kegg2 <- goregion(sigDMRs , all.cpg=rownames(methyl_raw), 
                             collection="GO", array.type="450K")
topGSA(gst.region.kegg2, number = 20)
