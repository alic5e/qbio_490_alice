setwd("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/analysis_data")

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
  BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE)) 
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
  if (! require(survival)){
    install.packages("survival")
  }
  library(survival)
  if (!require(survminer)){
    install.packages("survminer")
  }
  library(survminer)
  if (!require(ggplot2)){
    install.packages("ggplot2")
  }
  library(ggplot2)

##import clinical data
clin_query <-clin_query <- GDCquery(project = "TCGA-BRCA", 
                                    data.category = "Clinical", 
                                    file.type = "xml")
clinic <- GDCprepare_clinic(clin_query, 
                            clinical.info = "patient")

##import rna data
rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

#GDCdownload(rna_query)

rna_se <- GDCprepare(rna_query)

#filter all the NA values of race
race_mask<-ifelse(clinic$race_list=="",F,T)
clinic<-clinic[race_mask,]
race_mask_2<-ifelse(!is.na(clinic$race_list),T,F)
clinic<-clinic[race_mask_2,]
#create and save boxplot of age of diagonosis for different races
jpeg("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/boxplot.jpg")
boxplot(formula=clinic$age_at_initial_pathologic_diagnosis~ clinic$race_list,
        data=clinic,
        xlab="Race",
        ylab="Age at Diagnosis",
        main="Age at Diagnosis vs Race")
dev.off()
#create a KM-survival plot for different races
clinic$survival_time<-ifelse(is.na(clinic$days_to_death),
                                         clinic$survival_time <-clinic$days_to_last_followup,
                                         clinic$survival_time<-clinic$days_to_death)
inf_mask <- ifelse(clinic$survival_time =="-Inf", F, T)
clinic<- clinic[inf_mask,]

clinic$death_event <-ifelse(clinic$vital_status=="Alive",clinic$death_event<-F,clinic$death_event<-T)
survival_na_mask<-ifelse(is.na(clinic$survival_time),F,T)
clinic<-clinic[survival_na_mask,]

survival_object <- Surv(time = clinic$survival_time,event = clinic$death_event)
fit_object <- survfit( survival_object~ clinic$race_list, data =
                         clinic)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = 
                         "right")
KM_plot_race <- survplot$plot + theme_bw() + theme(axis.title = 
                                                    element_text(size=20), axis.text = element_text(size=16),
                                                  legend.title = element_text(size=14), legend.text = 
                                                    element_text(size=12))
KM_plot_race
jpeg("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/KM_plot_race.jpg")
KM_plot_race<-survplot$plot + theme_bw() + theme(axis.title= element_text (size=20),
                                                axis.text = element_text (size=16),
                                                legend.title = element_text(size=14),
                                                legend.text = element_text(size=12))
KM_plot_race
dev.off()

##import MAF
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

#filter the NA values
na_mask<-ifelse(maf_object@clinical.data$race_list=='',F,T)
maf_object@clinical.data<-maf_object@clinical.data[na_mask,]
na_mask_2<-ifelse(!is.na(maf_object@clinical.data$race_list),T,F)
maf_object@clinical.data<-maf_object@clinical.data[na_mask_2,]
#create oncoplot to show top 10 mutated genes for all patients
oncoplot(maf = maf_object,
         top = 10,
         clinicalFeatures = "race_list",
         borderCol = NA)

#subset small MAF according to races
asian_mask<-ifelse(maf_object@clinical.data$race_list=='ASIAN',T,F)
asian_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[asian_mask]
asian_maf <- subsetMaf(maf = maf_object,
                     tsb = asian_patient_barcodes)

white_mask<-ifelse(maf_object@clinical.data$race_list=='WHITE',T,F)
white_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[white_mask]
white_maf <- subsetMaf(maf = maf_object,
                       tsb = white_patient_barcodes)

black_mask<-ifelse(maf_object@clinical.data$race_list=='BLACK OR AFRICAN AMERICAN',T,F)
black_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[black_mask]
black_maf <- subsetMaf(maf = maf_object,
                       tsb = black_patient_barcodes)

native_mask<-ifelse(maf_object@clinical.data$race_list=='AMERICAN INDIAN OR ALASKA NATIVE',T,F)
native_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[native_mask]
native_maf <- subsetMaf(maf = maf_object,
                       tsb = native_patient_barcodes)

#create oncoplot to show top 10 mutated genes for White patients
oncoplot(maf = white_maf,
         top = 10,
         borderCol = NA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/oncoplot_white.png")

#create oncoplot to show top 10 mutated genes for Asian patients
oncoplot(maf = asian_maf,
         top = 10,
         borderCol = NA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/oncoplot_asian.png")
#create oncoplot to show top 10 mutated genes for Black patients
oncoplot(maf = black_maf,
         top = 10,
         borderCol = NA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/oncoplot_black.png")
#create oncoplot to show top 10 mutated genes for Native patients
oncoplot(maf = native_maf,
         top = 10,
         borderCol = NA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/midterm_project_wang/outputs/oncoplot_native.png")
##process rna
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)

row.names(rna_genes)<-rna_genes$gene_id
colnames(rna_counts)<-rna_clinical$barcode
rownames(rna_counts)<-rna_genes$gene_id

tissue_mask<-ifelse(rna_clinical$definition=='Solid Tissue Normal',F,T)
rna_clinical<-rna_clinical[tissue_mask,]
rna_counts<-rna_counts[,tissue_mask]

rna_clinical$race <- factor(rna_clinical$race)
head(rna_clinical$race)
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)
rna_clinical$definition <- factor(rna_clinical$definition)
head(rna_clinical$ajcc_pathologic_stage)
head(rna_clinical$definition)


ajcc_mask<-ifelse(!is.na(rna_clinical$ajcc_pathologic_stage),T,F)
rna_clinical<-rna_clinical[ajcc_mask,]
rna_counts<-rna_counts[,ajcc_mask]

sum(is.na(rna_clinical$age_category))
sum(is.na(rna_clinical$ajcc_pathologic_stage))
sum(is.na(rna_clinical$definition))
#create category to distinguish white and non-white
rna_clinical$race_white<-ifelse(rna_clinical$race=='white','white','non-white')

row_sums <-rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums<10,F,T)
rna_counts <- rna_counts[low_counts_mask,]
rna_genes <-rna_genes[low_counts_mask,]
BiocManager::install("DESeq2")
library(DESeq2)
#apply differential expression analysis
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ajcc_pathologic_stage + definition + race_white)

?DESeq
dds_obj <- DESeq(dds) 
?resultsNames
resultsNames(dds_obj)  

?results
results <- results(dds_obj, format = "DataFrame", contrast = c("race_white", "non-white", "white"))

genename_mask<-ifelse(rna_genes$gene_id %in% results@rownames,T,F)
rna_genes<-rna_genes[genename_mask,]
results <- data.frame(rna_genes$gene_name,rna_genes$gene_id,results$log2FoldChange,results$pvalue,results$padj,-log10(results$padj))
colnames(results) <- c("gene_name","gene_id","log2foldchange","pvalue","padj","-log10padj")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
#draw volcano plot to find out gene expression differences between groups
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2foldchange',
                y = '-log10padj')
               

