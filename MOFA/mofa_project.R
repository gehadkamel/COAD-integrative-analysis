
#Upload libraries for MOFA 
library(MOFA2)
library(MOFAdata)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(NMF)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(Seuray)
#set local working directory 
setwd("D:/Bioinformatics/Nile University/Integrative/Project")
#Upload
#########1 upload data 
#methyldata
methyl_data <- read.csv("Data/methy", fill = T, header = T, sep = "")
#mRNA 
mRNA_data <- read.csv("Data/exp", fill = T, header = T, sep = "")
#mirna
mirna_data <- read.csv("Data/mirna", fill = T, header = T, sep ="")
#Load clinical data and survival data
clinical_mofa <- read.delim("Data/clinical/colon", fill = T , header = T)
#adjust sample names
clinical_mofa[,1] <- gsub("-", ".", clinical_mofa[,1])
rownames(clinical_mofa) <- clinical_mofa[ ,1]
#load survival data
survival_mofa <- read.csv("Data/survival", fill =T , header = T, sep = "")
survival_mofa[,1] <- gsub("-", ".", survival_mofa[,1])
rownames(survival_mofa) <- survival_mofa[,1]
#get intersecting samples with clinical data
survival_mofa <- survival_mofa[which(rownames(survival_mofa) %in% rownames(clinical_mofa)),]
survival_mofa <- survival_mofa[match(rownames(survival_mofa),rownames(clinical_mofa)), ]

all(rownames(survival_mofa) == rownames(clinical_mofa))


#extract intersecting sample names are the same for all omics data and clinical data

samp_names <- intersect(colnames(mRNA_data), colnames(mirna_data))
samp_names <- intersect(samp_names, colnames(methyl_data))
clinical_mofa_f <- clinical_mofa[rownames(clinical_mofa)%in% samp_names,] 
#intersecting samples between 3 OMICS
survival_mofa_f <- survival_mofa[which(rownames(survival_mofa) %in% rownames(clinical_mofa_f)),]
survival_mofa_f <- survival_mofa_f[match(rownames(survival_mofa_f), rownames(clinical_mofa_f)), ] 

#########2 Explore and clean each OMIC data
#mRNA data
gene_mRNA <- rownames(mRNA_data)
#Convert data values to integers
mRNA_data=apply(round(mRNA_data),2,as.integer)
rownames(mRNA_data) <- gene_mRNA
#Check for empty cells
sum(is.null(mRNA_data))
sum(is.na(mRNA_data))
sum(is.nan(as.matrix(mRNA_data)))
#all(rownames(clinical_data_f) == colnames(methyl_data))
hist(as.matrix(mRNA_data))
hist(log2(as.matrix(mRNA_data+1)))
# QQ plot for the normality
qqnorm(mRNA_data[60,])
qqline(mRNA_data[60,])

#Normalize GE data
#match clinical data to GFE data
clinical_mofa <- clinical_mofa[which(rownames(clinical_mofa) %in% colnames(mRNA_data)), ]
#match and check if they are matching
clinical_mofa <- clinical_mofa[match(colnames(mRNA_data), rownames(clinical_mofa)), ]
all(rownames(clinical_mofa) == colnames(mRNA_data))
#Pre-filter
#remove low count genes
#keep <- rowsum(counts(as.matrix(mRNA_data)) >= 10)
#make DEseq data object for normalization
#Source http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering


#nomralized_mRNA <- NormalizeData(mRNA_data, normalization.method = "LogNormalize", assay = "RNA")
#normalized_mRNA <- ScaleData(mRNA_data, do.center = TRUE, do.scale = FALSE)

#normalized_mRNA <- FindVariableFeatures(normalized_mRNA, 
                                
  #                             nfeatures = 5000,
 # )


#Remove 

dds_mofa <- DESeqDataSetFromMatrix(countData = mRNA_data,
                              colData = clinical_mofa,
                             design = ~ vital_status)
dds_mofa
#filter low count genes
keep <- rowSums(counts(dds_mofa)) >= 10
dds_mofa <- dds_mofa[keep,]

dds_mofa <- estimateSizeFactors(dds_mofa)
#m <- estimateSizeFactorsForMatrix(mRNA_data)
sizeFactors(dds_mofa)

#m <- counts(m , normalized = TRUE)
normalized_mRNA <- counts(dds_mofa, normalized=TRUE)
hist(log2(normalized_mRNA+1))
hist(normalized_mRNA)
qqnorm(normalized_mRNA[60,])
qqline(normalized_mRNA[60,])

##########################
#miRNA  data
View(mirna_data)
genes_mirna <- rownames(mirna_data)
mirna_data <- apply(round(mirna_data),2 , as.integer)
rownames(mirna_data) <- genes_mirna
#make sure samples names are in the same order
all(rownames(clinical_mofa_f) == colnames(mirna_data))
qqnorm(mirna_data[60,])
qqline(mirna_data[60,])

dds_mirna <- DESeqDataSetFromMatrix(countData = mirna_data,
                                   colData = clinical_mofa_f,
                                   design = ~ vital_status)
dds_mirna
#filter low count genes
keep <- rowSums(counts(dds_mirna)) >= 10
dds_mirna <- dds_mirna[keep,]

dds_mirna <- estimateSizeFactors(dds_mirna)
#m <- estimateSizeFactorsForMatrix(mRNA_data)
sizeFactors(dds_mirna)

normalized_mirna <- counts(dds_mirna, normalized = TRUE)
qqnorm(normalized_mirna[50,])
qqline(normalized_mirna[50,])

#######Working on methyl data 

#Get intersecting samples 
methyl_data <- methyl_data[, which(colnames(methyl_data) %in% rownames(clinical_mofa_f))]
methyl_data <- methyl_data[, match(colnames(methyl_data), rownames(clinical_mofa_f))]
all(colnames(methyl_data) == rownames(clinical_mofa_f))
#Filter methyl data
sum(is.na(methyl_data))

gr_methyl <- makeGenomicRatioSetFromMatrix(as.matrix(methyl_data))

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(gr_methyl) %in% ann450k$Name[ann450k$chr %in% 
                                                  c("chrX","chrY")])
gr_methyl <- gr_methyl[keep,]
#Filter probes with known SNPs
gr_methyl <- dropLociWithSnps(gr_methyl)
#Filter out cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)

keep_c <- !(featureNames(gr_methyl) %in% xReactiveProbes$TargetID)
table(keep_c)
gr_methyl <- gr_methyl[keep_c,]
methyl_beta <- getBeta(gr_methyl)
methyl_m <- getM(gr_methyl)                                           

#############################################################################################
############################################################################################
#Start integrative analysis 

#Explore OMICs data
normalized_mRNA <- normalized_mRNA[,which(colnames(normalized_mRNA) %in% rownames(clinical_mofa_f))]
normalized_mRNA <- normalized_mRNA[, match(colnames(normalized_mRNA), rownames(clinical_mofa_f))]
hist(normalized_mRNA)
hist(normalized_mirna)
hist(methyl_beta)
#Ensuring all data have the same samples in the same order
all(colnames(normalized_mRNA) == colnames(normalized_mirna))
all(colnames(normalized_mRNA) == colnames(methyl_m))
all(colnames(methyl_beta) == rownames(clinical_mofa_f))
################################################################
mofa_data <- make_example_data(
  n_views = 2, 
  n_samples = 200, 
  n_features = 1000, 
  n_factors = 10
)[[1]]
lapply(mofa_data,dim)
#Meta data
which(colnames(clinical_mofa_f) == "vital_status")

clinical_meta <- clinical_mofa_f[,c(1,117)]
colnames(clinical_mofa_f)[1] <- "sample"
colnames(clinical_meta) <- c("sample", "vital_status") 
#load data into mofa_data
rownames(normalized_mRNA) <- gsub("\\|.", "", rownames(normalized_mRNA))
mofa_data[[1]] <- normalized_mRNA
mofa_data[[2]] <- mirna_data
mofa_data[[3]] <- methyl_m
names(mofa_data) <- c("mRNA", "miRNA", "Methylation")
#Create mofa object
MOFA_object <- create_mofa(mofa_data, extract_metadata = F)
#samples_metadata(MOFA_object) <- clinical_meta
#Adjust column data so that it can accept sample meta dfata
#colnames(survival_mofa_f) <- c("sample", "Survival", "Death")  
#samples_metadata(MOFA_object) <- survival_mofa_f
#omic_data <- makeExampleDESeqDataSet(gene = normalized_MRNA, 
 #                 rna = normalized_mirna,
  #                meth = methyl_ma, 
   #               surv = survival_mofa)

plot_data_overview(MOFA_object)
#data options
data_opts <- get_default_data_options(MOFA_object)
data_opts
#Model opts
model_opts <- get_default_model_options(MOFA_object)
model_opts$num_factors <- 10
#Define training options
train_opts <- get_default_training_options(MOFA_object)
train_opts$convergence_mode <- "slow" #after exploration of my model
#prepare MOFA Object
MOFA_object <- prepare_mofa(
  object = MOFA_object,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Run mofa
outfile = file.path(getwd(),"model.hd5")
MOFA_object <- run_mofa(MOFA_object, outfile, use_basilisk = T)

#MOFA_object_2 <- run_mofa(MOFA_object, outfile, use_basilisk = T)
#sample_met
samples_metadata(MOFA_object) <- clinical_mofa_f
#Calculate the variance explained per factor
r2 <- get_variance_explained(MOFA_object)
r2$r2_total

#plot variance
plot_variance_explained(MOFA_object)
#Show factor

plot_factor(MOFA_object,
            factor = 1:3,
            shape_by= "vital_status")


plot_factor(MOFA_object, 
                 factors = c(1,2,3),
                 color_by = "vital_status",
                 dot_size = 3,        
                 dodge = T,           
                 legend = F,          
                 add_violin = T,      
                 violin_alpha = 0.25  
)

#Plot heatmap
plot_data_heatmap(
  MOFA_object, 
  view = "Methylation", 
  factor = 1, 
  features = 50, 
  show_rownames = FALSE,
  denoise = T,
  annotation_samples = "pathologic_stage"
)


#Visualization of multiple factors
plot_factor(MOFA_object,
             factors = 1:3,
             color_by = "vital_status")


#Plot top features associated with each factor
plot_top_weights(MOFA_object,
                 view = "Methylation", 
                 factor = 1,
                 nfeatures = 10)
#top factors for  mRNA
plot_top_weights(MOFA_object,
                 view = "mRNA", 
                 factor = 3,
                 nfeatures = 10)
#Top factors for miRNA
plot_top_weights(MOFA_object,
                 view = "miRNA", 
                 factor = 3,
                 nfeatures = 10)
#plot scatter plot to see direct relationship betweeen weights and factor

plot_data_scatter(MOFA_object,
                  view = "miRNA",         # view of interest
                  factor = 3,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "vital_status"
)

#Enrichment analysis