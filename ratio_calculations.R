library(ggplot2)
library(dplyr)
library(cowplot)

ratio_dat<-data.frame(DATA)

ratio_dat$HT_01_04_mean_TPM<- rowMeans(ratio_dat[,c("HT1A_01_04_TPM","HT2A_01_04_TPM","HT3A_01_04_TPM")])
ratio_dat$LT_01_04_mean_TPM<- rowMeans(ratio_dat[,c("LT1A_01_04_TPM","LT2A_01_04_TPM","LT3A_01_04_TPM")])
ratio_dat$HT_18_05_mean_TPM<- rowMeans(ratio_dat[,c("HT1A_18_05_TPM","HT2A_18_05_TPM","HT3A_18_05_TPM")])
ratio_dat$LT_18_05_mean_TPM<- rowMeans(ratio_dat[,c("LT1A_18_05_TPM","LT2A_18_05_TPM","LT3A_18_05_TPM")])

ratio_dat$HT_01_04_mean_NumReads<- rowMeans(ratio_dat[,c("HT1A_01_04_NumReads","HT2A_01_04_NumReads","HT3A_01_04_NumReads")])
ratio_dat$LT_01_04_mean_NumReads<- rowMeans(ratio_dat[,c("LT1A_01_04_NumReads","LT2A_01_04_NumReads","LT3A_01_04_NumReads")])
ratio_dat$HT_18_05_mean_NumReads<- rowMeans(ratio_dat[,c("HT1A_18_05_NumReads","HT2A_18_05_NumReads","HT3A_18_05_NumReads")])
ratio_dat$LT_18_05_mean_NumReads<- rowMeans(ratio_dat[,c("LT1A_18_05_NumReads","LT2A_18_05_NumReads","LT3A_18_05_NumReads")])

ratio_dat$ratio_01_04_HT_LT<-ratio_dat$HT_01_04_mean_TPM/ratio_dat$LT_01_04_mean_TPM
ratio_dat$ratio_01_04_LT_HT<-ratio_dat$LT_01_04_mean_TPM/ratio_dat$HT_01_04_mean_TPM
ratio_dat$ratio_18_05_HT_LT<-ratio_dat$HT_18_05_mean_TPM/ratio_dat$LT_18_05_mean_TPM
ratio_dat$ratio_18_05_LT_HT<-ratio_dat$LT_18_05_mean_TPM/ratio_dat$HT_18_05_mean_TPM

ratio_dat <- subset(ratio_dat,HT_01_04_mean_NumReads > 100 | 
                               LT_01_04_mean_NumReads > 100 |
                              HT_18_05_mean_NumReads > 100 |
                                LT_18_05_mean_NumReads > 100)

ratio_dat <- subset(ratio_dat,ratio_01_04_HT_LT > 2 | 
                      ratio_01_04_LT_HT > 2 |
                      ratio_18_05_HT_LT > 2 |
                      ratio_18_05_LT_HT > 2)
ratio_dat_asco <- subset(ratio_dat,phylum=="Ascomycota")