Alice


library(BiocManager) 
library(TCGAbiolinks) 
library(maftools)

setwd("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/analysis_data")
clinic <- read.csv("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/analysis_data/brca_clinical_data.csv")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")
##1. age at initial diagnosis
##2. discrete
##3. therapy_types:different therpy types received by the patients
##4. categorical 
##5. Older patients may choose more conservative therpy.
## The younger patients at the time of diagnosis, the higher the survival rate.
## The more conservative methods may result in lower survival rates.

clinical_full <- merge(clinic,clinical_drug,by="bcr_patient_barcode")

boxplot(formula=clinical_full$age_at_initial_pathologic_diagnosis~ clinical_full$therapy_types,
        data=clinical_full,
        xlab="Therapy Types",
        ylab="Age at Diagnosis",
        main="Age at Diagnosis vs Therapy Types")

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

age_na_mask <- ifelse(!is.na(clinic$age_at_initial_pathologic_diagnosis),TRUE,FALSE)
age_cleaned_clinic<-clinic[age_na_mask,]

young_mask<- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis<=35,TRUE,FALSE)
middle_mask<- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis>35&age_cleaned_clinic$age_at_initial_pathologic_diagnosis<50,TRUE,FALSE)
old_mask<- ifelse(age_cleaned_clinic$age_at_initial_pathologic_diagnosis>=50,TRUE,FALSE)
age_cleaned_clinic$age_status<- ifelse(young_mask,"Young",ifelse(middle_mask,"Middle","Old"))

age_cleaned_clinic$survival_time<-ifelse(is.na(age_cleaned_clinic$days_to_death),
                                           age_cleaned_clinic$survival_time <-age_cleaned_clinic$days_to_last_followup,
                                           age_cleaned_clinic$survival_time<-age_cleaned_clinic$days_to_death)

inf_mask <- ifelse(age_cleaned_clinic$survival_time =="-Inf", F, T)
age_cleaned_clinic<- age_cleaned_clinic[inf_mask,]

age_cleaned_clinic$death_event <-ifelse(age_cleaned_clinic$vital_status=="Alive",age_cleaned_clinic$death_event<-F,age_cleaned_clinic$death_event<-T)
survival_na_mask<-ifelse(is.na(age_cleaned_clinic$survival_time),F,T)
age_cleaned_clinic<-age_cleaned_clinic[survival_na_mask,]


survival_object <- Surv(time = age_cleaned_clinic$survival_time,event = age_cleaned_clinic$death_event)
fit_object <- survfit( survival_object~ age_cleaned_clinic$age_status, data =
                        age_cleaned_clinic)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = 
                         "right")
KM_plot_age <- survplot$plot + theme_bw() + theme(axis.title = 
                                                  element_text(size=20), axis.text = element_text(size=16),
                                                legend.title = element_text(size=14), legend.text = 
                                                  element_text(size=12))
KM_plot_age
jpeg("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/KM_plot_age.jpg")
KM_plot_age<-survplot$plot + theme_bw() + theme(axis.title= element_text (size=20),
                                                axis.text = element_text (size=16),
                                                legend.title = element_text(size=14),
                                                legend.text = element_text(size=12))
KM_plot_age
dev.off()

therapy_na_mask <- ifelse(clinical_full$therapy_types==""|is.na(clinical_full$therapy_types),F,T)
therapy_cleaned_clinic<-clinical_full[therpy_na_mask,]

therapy_cleaned_clinic$survival_time<-ifelse(is.na(therapy_cleaned_clinic$days_to_death),
                                            therapy_cleaned_clinic$survival_time <-therapy_cleaned_clinic$days_to_last_followup,
                                            therapy_cleaned_clinic$survival_time<-therapy_cleaned_clinic$days_to_death)

inf_mask <- ifelse(therapy_cleaned_clinic$survival_time =="-Inf", F, T)
therapy_cleaned_clinic<- therapy_cleaned_clinic[inf_mask,]

therapy_cleaned_clinic$death_event <-ifelse(therapy_cleaned_clinic$vital_status=="Alive",therapy_cleaned_clinic$death_event<-F,therpy_cleaned_clinic$death_event<-T)
survival_na_mask_2<-ifelse(is.na(therapy_cleaned_clinic$survival_time),F,T)
therapy_cleaned_clinic<-therapy_cleaned_clinic[survival_na_mask_2,]
therapy_cleaned_clinic$therapy_types<-as.character(therapy_cleaned_clinic$therapy_types)

survival_object <- Surv(time = therapy_cleaned_clinic$survival_time,event = therapy_cleaned_clinic$death_event)
fit_object <- survfit( survival_object~ therapy_cleaned_clinic$therapy_types, data =
                         therapy_cleaned_clinic)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme = 
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = 
                         "right")
KM_plot_therapy <- survplot$plot + theme_bw() + theme(axis.title = 
                                                    element_text(size=20), axis.text = element_text(size=16),
                                                  legend.title = element_text(size=14), legend.text = 
                                                    element_text(size=12))
KM_plot_therapy
jpeg("/Users/aliceeeee/Desktop/Spring23/QBIO490/qbio_490_alice/KM_plot_therapy.jpg")
KM_plot_therapy<-survplot$plot + theme_bw() + theme(axis.title= element_text (size=20),
                                                axis.text = element_text (size=16),
                                                legend.title = element_text(size=14),
                                                legend.text = element_text(size=12))
KM_plot_therapy
dev.off()






