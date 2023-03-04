setwd("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/analysis_data")
if(!require(TCGAbiolinks)){
  install.packages("TCGAbiolinks")
}
library(TCGAbiolinks) 
if(!require(maftools)){
  install.packages("maftools")
}
library(maftools)
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

clinical <- read.csv("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/analysis_data/brca_clinical_data.csv")

maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) 

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

clinical$gender<-factor(clinical$gender)


oncoplot(maf = maf_object,
         top = 15,
         clinicalFeatures = "gender",
         borderCol = NA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/week7_maf/oncoplot_sex.png")

#A large portion of male patients have mutations in PIK3CA gene.A gene that makes one of the proteins in an enzyme called PI3K, which is involved in many important functions in a cell. Mutations (changes) in the PIK3CA gene may cause the PI3K enzyme to become overactive, which may cause cancer cells to grow.
#reason: PIK3CA correlates more with male breast cancer than female breast cancer.

genePIK3CA_maf <- subsetMaf(maf = maf_object,
                       genes = 'PIK3CA')
PIK3CA_clinical<-genePIK3CA_maf@clinical.data
PIK3CAgender_mask<-ifelse(PIK3CA_clinical$gender=='FEMALE',T,F)
female_PIK3CA_clinical<-PIK3CA_clinical[PIK3CAgender_mask,]
male_PIK3CA<-nrow(PIK3CA_clinical)-nrow(female_PIK3CA_clinical)
female_PIK3CA<-nrow(female_PIK3CA_clinical)

gender_mask<-ifelse(clinical$gender=='FEMALE',T,F)
female_clinical<-clinical[gender_mask,]
female_without_PIK3CA<-nrow(female_clinical)-nrow(female_PIK3CA_clinical)
male_without_PIK3CA<-nrow(clinical)-nrow(female_clinical)-male_PIK3CA

contig_sex_PIK3CA <- matrix(c(female_PIK3CA, 
                   female_without_PIK3CA,
                   male_PIK3CA,
                   male_without_PIK3CA), 
                 nrow=2)

contig_sex_PIK3CA

mosaicplot(contig_sex_PIK3CA)
ggsave("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/week7_maf/contig_sex_PIK3CA.png")

fisher_test <- fisher.test(contig_sex_PIK3CA)
fisher_test
#p-value=1 no significance for male in the PIK3CA gene mutation
#odds ratio 1.29556 - being a woman patient means 1.5x higher rate in PIK3CA gene mutation

female_mask<-ifelse(maf_object@clinical.data$gender=='FEMALE',T,F)
female_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[female_mask]

female_maf <- subsetMaf(maf = maf_object,
                       tsb = female_patient_barcodes)

male_mask<-ifelse(maf_object@clinical.data$gender=='MALE',T,F)
male_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[male_mask]

male_maf <- subsetMaf(maf = maf_object,
                        tsb = male_patient_barcodes)


lollipopPlot2(m1 = female_maf, 
              m2 = male_maf, 
              m1_name ='Female',
              m2_name ='Male',
              gene = "PIK3CA")
#male have less mutation locations than female, the mutation locations of male are the most frequent locations of female.

maf_object@clinical.data$Overall_Survival_Status<-ifelse(maf_object@clinical.data$vital_status=='Alive',T,F)
mafSurvival(maf = maf_object,
            genes = "PIK3CA", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)
#Patients without mutations in PIK3CA gene have a higher survival rate after 2000 days.It's probably because PIK3 enzyme involves in patients' long term suvival.



