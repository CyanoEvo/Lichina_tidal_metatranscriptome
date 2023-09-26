#use the quant files from Salmon imported as HT1A_01_04_dat etc
gene_name<-HT1A_01_04_dat$Name
HT1A_01_04_asc_raw<-HT1A_01_04_dat$HT1A_01_04_asc_NumReads
HT2A_01_04_asc_raw<-HT2A_01_04_dat$HT2A_01_04_asc_NumReads
HT3A_01_04_asc_raw<-HT3A_01_04_dat$HT3A_01_04_asc_NumReads
LT1A_01_04_asc_raw<-LT1A_01_04_dat$LT1A_01_04_asc_NumReads
LT2A_01_04_asc_raw<-LT2A_01_04_dat$LT2A_01_04_asc_NumReads
LT3A_01_04_asc_raw<-LT3A_01_04_dat$LT3A_01_04_asc_NumReads
asc_deseq<-data.frame(gene_name,HT1A_01_04_asc_raw,HT2A_01_04_asc_raw,HT3A_01_04_asc_raw,LT1A_01_04_asc_raw,LT2A_01_04_asc_raw,LT3A_01_04_asc_raw)
asc_deseq2 <- asc_deseq[,-1]
rownames(asc_deseq2) <- asc_deseq[,1]
asc_deseq<-asc_deseq2

nos_gene_name<-HT1A_01_04_nos_dat$Name
HT1A_01_04_nos_raw<-HT1A_01_04_nos_dat$HT1A_01_04_nos_NumReads
HT2A_01_04_nos_raw<-HT2A_01_04_nos_dat$HT2A_01_04_nos_NumReads
HT3A_01_04_nos_raw<-HT3A_01_04_nos_dat$HT3A_01_04_nos_NumReads
LT1A_01_04_nos_raw<-LT1A_01_04_nos_dat$LT1A_01_04_nos_NumReads
LT2A_01_04_nos_raw<-LT2A_01_04_nos_dat$LT2A_01_04_nos_NumReads
LT3A_01_04_nos_raw<-LT3A_01_04_nos_dat$LT3A_01_04_nos_NumReads
nos_deseq<-data.frame(nos_gene_name,HT1A_01_04_nos_raw,HT2A_01_04_nos_raw,HT3A_01_04_nos_raw,LT1A_01_04_nos_raw,LT2A_01_04_nos_raw,LT3A_01_04_nos_raw)
nos_deseq2 <- nos_deseq[,-1]
rownames(nos_deseq2) <- nos_deseq[,1]
nos_deseq<-nos_deseq2

ple_gene_name<-HT1A_01_04_ple_dat$Name
HT1A_01_04_ple_raw<-HT1A_01_04_ple_dat$HT1A_01_04_ple_NumReads
HT2A_01_04_ple_raw<-HT2A_01_04_ple_dat$HT2A_01_04_ple_NumReads
HT3A_01_04_ple_raw<-HT3A_01_04_ple_dat$HT3A_01_04_ple_NumReads
LT1A_01_04_ple_raw<-LT1A_01_04_ple_dat$LT1A_01_04_ple_NumReads
LT2A_01_04_ple_raw<-LT2A_01_04_ple_dat$LT2A_01_04_ple_NumReads
LT3A_01_04_ple_raw<-LT3A_01_04_ple_dat$LT3A_01_04_ple_NumReads
ple_deseq<-data.frame(ple_gene_name,HT1A_01_04_ple_raw,HT2A_01_04_ple_raw,HT3A_01_04_ple_raw,LT1A_01_04_ple_raw,LT2A_01_04_ple_raw,LT3A_01_04_ple_raw)
ple_deseq2 <- ple_deseq[,-1]
rownames(ple_deseq2) <- ple_deseq[,1]
ple_deseq<-ple_deseq2


#DESEQ
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)

deseq_cols <- data.frame (name  = c("HT1A_01_04_asc_NumReads", "HT2A_01_04_asc_NumReads", "HT2A_01_04_asc_NumReads","LT1A_01_04_asc_NumReads", "LT2A_01_04_asc_NumReads", "LT2A_01_04_asc_NumReads"),
                          treatment = c("HT", "HT", "HT","LT","LT","LT"))

dds <- DESeqDataSetFromMatrix(countData = round(asc_deseq),
                              colData = deseq_cols,
                              design = ~ treatment)
dds <- DESeq(dds)

res <- results(dds)
res<-data.frame(res)
res$name<-row.names(res)

asc_deseq_plot<-res %>% 
  mutate(Significant = padj < 0.05) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col=Significant, label=name)) + 
  geom_point() +
  scale_colour_manual(values = c("grey", "#d60066"))+
  geom_vline(xintercept = 0, linetype="dashed", color = "red", size=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", color = "black", size=0.5) +
  annotate("text", x = -0.7883858-1.5, y = -log10(0.01304296), label = "GH3.1") +
  geom_point(aes(x= -0.7883858, -log10(0.01304296)), colour="black",shape=4) +
  annotate("text", x = -0.967371-1.5, y = -log10(0.000260902), label = "GH3.2") +
  geom_point(aes(x= -0.967371, -log10(0.000260902)), colour="black",shape=4) +
  annotate("text", x = -1.623542-1.5, y = -log10(0.00001176753), label = "GH3.3") +
  geom_point(aes(x= -1.623542, -log10(0.00001176753)), colour="black",shape=4) +
  annotate("text", x = -1.452382-1.5, y = -log10(0.001377535), label = "GH16.2") +
  geom_point(aes(x= -1.452382, -log10(0.001377535)), colour="black",shape=4) +
  annotate("text", x = -1.47644-1.5, y = -log10(0.002627178), label = "GH16.3") +
  geom_point(aes(x= -1.47644, -log10(0.002627178)), colour="black",shape=4) +
  annotate("text", x = 0.7729482-1.5, y = -log10(0.002518139), label = "GH17.2") +
  geom_point(aes(x= 0.7729482, -log10(0.002518139)), colour="black",shape=4) +
  annotate("text", x = -1.508467-1.5, y = -log10(3.279238E-14), label = "GH65") +
  geom_point(aes(x= -1.508467, -log10(3.279238E-14)), colour="black",shape=4) +
  annotate("text", x = 0.7624655-1.5, y = -log10(0.007989664), label = "GH115") +
  geom_point(aes(x= 0.7624655, -log10(0.007989664)), colour="black",shape=4) +
  annotate("text", x = -1.336588-1.5, y = -log10(0.00000009669841), label = "ST.1") +
  geom_point(aes(x= -1.336588, -log10(0.00000009669841)), colour="black",shape=4) +
  annotate("text", x = 0.6479942-1.5, y = -log10(0.01046288), label = "ST.4") +
  geom_point(aes(x= 0.6479942, -log10(0.01046288)), colour="black",shape=4) +
  annotate("text", x = 0.6512495-1.5, y = -log10(0.002992756), label = "ST.6") +
  geom_point(aes(x= 0.6512495, -log10(0.002992756)), colour="black",shape=4) +
  xlim(-12, 12)+
  theme_linedraw() + 
  theme(aspect.ratio=1, legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#rivularia
deseq_nos_cols <- data.frame (name  = c("HT1A_01_04_nos_NumReads", "HT2A_01_04_nos_NumReads", "HT2A_01_04_nos_NumReads","LT1A_01_04_nos_NumReads", "LT2A_01_04_nos_NumReads", "LT2A_01_04_nos_NumReads"),
                          treatment = c("HT", "HT", "HT","LT","LT","LT"))

dds_nos <- DESeqDataSetFromMatrix(countData = round(nos_deseq),
                              colData = deseq_nos_cols,
                              design = ~ treatment)
dds_nos <- DESeq(dds_nos)

res_nos <- results(dds_nos)
res_nos<-data.frame(res_nos)
res_nos$name<-row.names(res_nos)

nos_deseq_plot<-res_nos %>% 
  mutate(Significant = padj < 0.05) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col=Significant, label=name)) + 
  geom_point() +
  scale_colour_manual(values = c("grey", "#008827"))+
  geom_vline(xintercept = 0, linetype="dashed", color = "red", size=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", color = "black", size=0.5) +
  annotate("text", x = -1.636132-1.5, y = -log10(0.0002017103), label = "mtdh.1") +
  geom_point(aes(x= -1.623542, -log10(0.0002017103)), colour="black",shape=4) +
  annotate("text", x = 2.244183-1.5, y = -log10(0.04421577), label = "kpsT.4") +
  geom_point(aes(x= 2.244183, -log10(0.04421577)), colour="black",shape=4) +
  annotate("text", x = 1.715334-1.5, y = -log10(0.04816992), label = "wzx.3") +
  geom_point(aes(x= 1.715334, -log10(0.04816992)), colour="black",shape=4) +
#  annotate("text", x = -1.336588-1.5, y = -log10(0.00000009669841), label = "ST.1") +
#  geom_point(aes(x= -1.336588, -log10(0.00000009669841)), colour="black",shape=4) +
  #xlim(-12, 12)+
  theme_linedraw() + 
  theme(aspect.ratio=1, legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#pleurocapsa
deseq_ple_cols <- data.frame (name  = c("HT1A_01_04_ple_NumReads", "HT2A_01_04_ple_NumReads", "HT2A_01_04_ple_NumReads","LT1A_01_04_ple_NumReads", "LT2A_01_04_ple_NumReads", "LT2A_01_04_ple_NumReads"),
                              treatment = c("HT", "HT", "HT","LT","LT","LT"))

dds_ple <- DESeqDataSetFromMatrix(countData = round(ple_deseq),
                                  colData = deseq_ple_cols,
                                  design = ~ treatment)
dds_ple <- DESeq(dds_ple)

res_ple <- results(dds_ple)
res_ple<-data.frame(res_ple)
res_ple$name<-row.names(res_ple)

ple_deseq_plot<-res_ple %>% 
  mutate(Significant = padj < 0.05) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col=Significant, label=name)) + 
  geom_point() +
  scale_colour_manual(values = c("grey", "#0a50a1"))+
  geom_vline(xintercept = 0, linetype="dashed", color = "red", size=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", color = "black", size=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", color = "black", size=0.5) +
  annotate("text", x = 1.331556-1.5, y = -log10(0.00002437791), label = "psbA.2") +
  geom_point(aes(x= 1.331556, -log10(0.00002437791)), colour="black",shape=4) +
  #annotate("text", x = -0.7858143-1.5, y = -log10(0.02399926), label = "Ggps.3") +
  #geom_point(aes(x= -0.7858143, -log10(0.02399926)), colour="black",shape=4) +
    annotate("text", x = 2.823285-1.5, y = -log10(0.001074782), label = "spsA.2") +
    geom_point(aes(x= 2.823285, -log10(0.001074782)), colour="black",shape=4) +
  annotate("text", x = 3.849179-1.5, y = -log10(0.00006196524), label = "spsA.3") +
  geom_point(aes(x= 3.849179, -log10(0.00006196524)), colour="black",shape=4) +
  annotate("text", x = 1.141301-1.5, y = -log10(0.005771914), label = "spsA.6") +
  geom_point(aes(x= 1.141301, -log10(0.005771914)), colour="black",shape=4) +
  annotate("text", x = -0.7858143-1.5, y = -log10(0.02399926), label = "ggpsA.3") +
  geom_point(aes(x= -0.7858143, -log10(0.02399926)), colour="black",shape=4) +
  annotate("text", x = 2.148677-1.5, y = -log10(4.399692E-10), label = "all1823.09") +
  geom_point(aes(x= 2.148677, -log10(4.399692E-10)), colour="black",shape=4) +
  annotate("text", x = 1.036511-1.5, y = -log10(0.02904725), label = "all1823.14") +
  geom_point(aes(x= 1.036511, -log10(0.02904725)), colour="black",shape=4) +
  annotate("text", x = -1.068749-1.5, y = -log10(0.002889949), label = "kpsD.1") +
  geom_point(aes(x= -1.068749, -log10(0.002889949)), colour="black",shape=4) +
  xlim(-9, 9)+
  theme_linedraw() + 
  theme(aspect.ratio=1, legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


library(cowplot)
plot_grid(nos_deseq_plot,ple_deseq_plot,asc_deseq_plot, nrow=1)