if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")

#https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/06_DGE_visualizing_results.html

n_occur <- data.frame(table(ple_diff_ex$Name))
n_occur[n_occur$Freq > 1,]

#PLEUROCAPSALES
count_matrix<-as.matrix(ple_diff_ex_READS)
count_matrix <- count_matrix[, -7]
count_matrix <- floor(count_matrix)
head(ple_diff_mat)
coldata <- data.frame(
  sample = c( "HT1A_01_04_ple_NumReads", "HT2A_01_04_ple_NumReads", "HT3A_01_04_ple_NumReads", "LT1A_01_04_ple_NumReads", "LT2A_01_04_ple_NumReads", "LT3A_01_04_ple_NumReads" ),
  condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)
all(rownames(coldata) %in% colnames(count_matrix))
all(rownames(coldata) == colnames(count_matrix))
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, 
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,]
# set control condition as reference
dds$condition <- relevel(dds$condition, ref = "LT")
dds <- DESeq(dds)

# see all comparisons (here there is only one)
resultsNames(dds)

res <- results(dds)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
res
res[order(res$padj),]  
summary(results(dds, alpha=0.05))

padj.cutoff <- 0.005

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE <- res_tb %>%
  filter(padj < padj.cutoff)

mov10_meta <- coldata %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

norm_OEsig <- normalized_counts[,c(1:4,5:7)] %>% 
  filter(gene %in% sigOE$gene) 

heat_colors <- inferno(256)

norm_OEsig2 <- norm_OEsig[,-1]
rownames(norm_OEsig2) <- norm_OEsig[,1]

top<-pheatmap(norm_OEsig2, 
         color = heat_colors, 
         cluster_cols = F, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = coldata, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


res_tb <- res_tb %>% 
  mutate(threshold_OE = padj < 0.005 & abs(log2FoldChange) >= 0.58)



ggplot(res_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  geom_text_repel(data=subset(res_tb, padj <0.001),
                  aes(x = log2FoldChange, y = -log10(padj), label = gene),size=3)+
 # ggtitle("Mov10 overexpression") +
 # xlab("log2 fold change") + 
#  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



#######RAINY
  count_matrix_rain<-as.matrix(ple_diff_ex_READS)
  count_matrix_rain <- count_matrix_rain[, -1]
  count_matrix_rain <- floor(count_matrix_rain)
  coldata_rain <- data.frame(
    sample = c( "HT1A_18_05_ple_NumReads", "HT2A_18_05_ple_NumReads", "HT3A_18_05_ple_NumReads", "LT1A_18_05_ple_NumReads", "LT2A_18_05_ple_NumReads", "LT3A_18_05_ple_NumReads" ),
    condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
    row.names = "sample" )
  coldata_rain$condition <- as.factor(coldata_rain$condition)
  all(rownames(coldata_rain) %in% colnames(count_matrix_rain))
  all(rownames(coldata_rain) == colnames(count_matrix_rain))
  dds_rain <- DESeqDataSetFromMatrix(countData = count_matrix_rain, colData = coldata_rain, 
                                design = ~ condition)
  dds_rain <- dds_rain[rowSums(counts(dds_rain)) >= 10,]
  # set control condition as reference
  dds_rain$condition <- relevel(dds_rain$condition, ref = "LT")
  dds_rain <- DESeq(dds_rain)
  
  # see all comparisons (here there is only one)
  resultsNames(dds_rain)
  
  res_rain <- results(dds_rain)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
  res_rain[order(res_rain$padj),]  
  summary(results(dds_rain, alpha=0.05))
  
  padj.cutoff <- 0.005
  
  res_tb_rain <- res_rain %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  sigOE_rain <- res_tb_rain %>%
    filter(padj < padj.cutoff)
  
  mov10_meta_rain <- coldata_rain %>% 
    rownames_to_column(var="samplename") %>% 
    as_tibble()
  
  normalized_counts_rain <- counts(dds_rain, normalized=T) %>% 
    data.frame() %>%
    rownames_to_column(var="gene") 
  
  norm_OEsig_rain <- normalized_counts_rain[,c(1:4,5:7)] %>% 
    filter(gene %in% sigOE_rain$gene) 
  
  heat_colors <- inferno(256)
  
  norm_OEsig_rain2 <- norm_OEsig_rain[,-1]
  rownames(norm_OEsig_rain2) <- norm_OEsig_rain[,1]
  
  bottom<-pheatmap(norm_OEsig_rain2, 
           color = heat_colors,
           cluster_cols = F, 
           cluster_rows = T, 
           show_rownames = T,
           annotation = coldata_rain, 
           border_color = NA, 
           fontsize = 10, 
           scale = "row", 
           fontsize_row = 10, 
           height = 20)
  
plot_grid(as.grob(top),as.grob(bottom),ncol=1)

#NOSTOCALES
nos_count_matrix<-as.matrix(nos_diff_ex_READS)
nos_count_matrix <- nos_count_matrix[, -7]
nos_count_matrix <- floor(nos_count_matrix)
nos_coldata <- data.frame(
  sample = c( "HT1A_01_04_nos_NumReads", "HT2A_01_04_nos_NumReads", "HT3A_01_04_nos_NumReads", "LT1A_01_04_nos_NumReads", "LT2A_01_04_nos_NumReads", "LT3A_01_04_nos_NumReads" ),
  condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
  row.names = "sample" )
nos_coldata$condition <- as.factor(nos_coldata$condition)
all(rownames(nos_coldata) %in% colnames(nos_count_matrix))
all(rownames(nos_coldata) == colnames(nos_count_matrix))
nos_dds <- DESeqDataSetFromMatrix(countData = nos_count_matrix, colData = nos_coldata, 
                              design = ~ condition)
nos_dds <- nos_dds[rowSums(counts(nos_dds)) >= 10,]
# set control condition as reference
nos_dds$condition <- relevel(nos_dds$condition, ref = "LT")
nos_dds <- DESeq(nos_dds)

# see all comparisons (here there is only one)
resultsNames(nos_dds)

nos_res <- results(nos_dds)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
nos_res
nos_res[order(nos_res$padj),]  
summary(results(nos_dds, alpha=0.05))

padj.cutoff <- 0.005

nos_res_tb <- nos_res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

nos_sigOE <- nos_res_tb %>%
  filter(padj < padj.cutoff)

nos_mov10_meta <- nos_coldata %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

nos_normalized_counts <- counts(nos_dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

nos_norm_OEsig <- nos_normalized_counts[,c(1:4,5:7)] %>% 
  filter(gene %in% nos_sigOE$gene) 

heat_colors <- inferno(256)

nos_norm_OEsig2 <- nos_norm_OEsig[,-1]
rownames(nos_norm_OEsig2) <- nos_norm_OEsig[,1]

nos_top<-pheatmap(nos_norm_OEsig2, 
         color = heat_colors, 
         cluster_cols = F, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = nos_coldata, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)


nos_res_tb1 <- nos_res_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

nos_very_sig<-subset(nos_res_tb1, padj <0.05)

ggplot(nos_res_tb1) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  geom_text_repel(data=subset(nos_res_tb1, padj <0.001),
                  aes(x = log2FoldChange, y = -log10(padj), label = gene),size=3)+
 # ggtitle("Mov10 overexpression") +
 # xlab("log2 fold change") + 
#  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



#######RAINY
  count_matrix_rain<-as.matrix(ple_diff_ex_READS)
  count_matrix_rain <- count_matrix_rain[, -1]
  count_matrix_rain <- floor(count_matrix_rain)
  coldata_rain <- data.frame(
    sample = c( "HT1A_18_05_ple_NumReads", "HT2A_18_05_ple_NumReads", "HT3A_18_05_ple_NumReads", "LT1A_18_05_ple_NumReads", "LT2A_18_05_ple_NumReads", "LT3A_18_05_ple_NumReads" ),
    condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
    row.names = "sample" )
  coldata_rain$condition <- as.factor(coldata_rain$condition)
  all(rownames(coldata_rain) %in% colnames(count_matrix_rain))
  all(rownames(coldata_rain) == colnames(count_matrix_rain))
  dds_rain <- DESeqDataSetFromMatrix(countData = count_matrix_rain, colData = coldata_rain, 
                                design = ~ condition)
  dds_rain <- dds_rain[rowSums(counts(dds_rain)) >= 10,]
  # set control condition as reference
  dds_rain$condition <- relevel(dds_rain$condition, ref = "LT")
  dds_rain <- DESeq(dds_rain)
  
  # see all comparisons (here there is only one)
  resultsNames(dds_rain)
  
  res_rain <- results(dds_rain)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
  res_rain[order(res_rain$padj),]  
  summary(results(dds_rain, alpha=0.05))
  
  padj.cutoff <- 0.005
  
  res_tb_rain <- res_rain %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  sigOE_rain <- res_tb_rain %>%
    filter(padj < padj.cutoff)
  
  mov10_meta_rain <- coldata_rain %>% 
    rownames_to_column(var="samplename") %>% 
    as_tibble()
  
  normalized_counts_rain <- counts(dds_rain, normalized=T) %>% 
    data.frame() %>%
    rownames_to_column(var="gene") 
  
  norm_OEsig_rain <- normalized_counts_rain[,c(1:4,5:7)] %>% 
    filter(gene %in% sigOE_rain$gene) 
  
  heat_colors <- inferno(256)
  
  norm_OEsig_rain2 <- norm_OEsig_rain[,-1]
  rownames(norm_OEsig_rain2) <- norm_OEsig_rain[,1]
  
  bottom<-pheatmap(norm_OEsig_rain2, 
           color = heat_colors,
           cluster_cols = F, 
           cluster_rows = T, 
           show_rownames = T,
           annotation = coldata_rain, 
           border_color = NA, 
           fontsize = 10, 
           scale = "row", 
           fontsize_row = 10, 
           height = 20)
  
plot_grid(as.grob(top),as.grob(nos_top),ncol=1)

#Ascomycota
asc_count_matrix<-as.matrix(asc_diff_ex_READS)
asc_count_matrix <- asc_count_matrix[, -7]
asc_count_matrix <- floor(asc_count_matrix)
asc_coldata <- data.frame(
  sample = c( "HT1A_01_04_asc_NumReads", "HT2A_01_04_asc_NumReads", "HT3A_01_04_asc_NumReads", "LT1A_01_04_asc_NumReads", "LT2A_01_04_asc_NumReads", "LT3A_01_04_asc_NumReads" ),
  condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
  row.names = "sample" )
asc_coldata$condition <- as.factor(asc_coldata$condition)
all(rownames(asc_coldata) %in% colnames(asc_count_matrix))
all(rownames(asc_coldata) == colnames(asc_count_matrix))
asc_dds <- DESeqDataSetFromMatrix(countData = asc_count_matrix, colData = asc_coldata, 
                                  design = ~ condition)
asc_dds <- asc_dds[rowSums(counts(asc_dds)) >= 10,]
# set control condition as reference
asc_dds$condition <- relevel(asc_dds$condition, ref = "LT")
asc_dds <- DESeq(asc_dds)

# see all comparisons (here there is only one)
resultsNames(asc_dds)

asc_res <- results(asc_dds)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
asc_res
asc_res[order(asc_res$padj),]  
summary(results(asc_dds, alpha=0.05))

padj.cutoff <- 0.005

asc_res_tb <- asc_res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

asc_sigOE <- asc_res_tb %>%
  filter(padj < padj.cutoff)

asc_mov10_meta <- asc_coldata %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

asc_normalized_counts <- counts(asc_dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

asc_norm_OEsig <- asc_normalized_counts[,c(1:4,5:7)] %>% 
  filter(gene %in% asc_sigOE$gene) 

heat_colors <- inferno(256)

asc_norm_OEsig2 <- asc_norm_OEsig[,-1]
rownames(asc_norm_OEsig2) <- asc_norm_OEsig[,1]

asc_top<-pheatmap(asc_norm_OEsig2, 
                  color = heat_colors, 
                  cluster_cols = F, 
                  cluster_rows = T, 
                  show_rownames = F,
                  annotation = asc_coldata, 
                  border_color = NA, 
                  fontsize = 10, 
                  scale = "row", 
                  fontsize_row = 10, 
                  height = 20)


nos_res_tb1 <- nos_res_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

nos_very_sig<-subset(nos_res_tb1, padj <0.05)

ggplot(nos_res_tb1) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  geom_text_repel(data=subset(nos_res_tb1, padj <0.001),
                  aes(x = log2FoldChange, y = -log10(padj), label = gene),size=3)+
  # ggtitle("Mov10 overexpression") +
  # xlab("log2 fold change") + 
  #  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



#######RAINY
count_matrix_rain<-as.matrix(ple_diff_ex_READS)
count_matrix_rain <- count_matrix_rain[, -1]
count_matrix_rain <- floor(count_matrix_rain)
coldata_rain <- data.frame(
  sample = c( "HT1A_18_05_ple_NumReads", "HT2A_18_05_ple_NumReads", "HT3A_18_05_ple_NumReads", "LT1A_18_05_ple_NumReads", "LT2A_18_05_ple_NumReads", "LT3A_18_05_ple_NumReads" ),
  condition = c( "HT", "HT",  "HT", "LT", "LT", "LT" ), 
  row.names = "sample" )
coldata_rain$condition <- as.factor(coldata_rain$condition)
all(rownames(coldata_rain) %in% colnames(count_matrix_rain))
all(rownames(coldata_rain) == colnames(count_matrix_rain))
dds_rain <- DESeqDataSetFromMatrix(countData = count_matrix_rain, colData = coldata_rain, 
                                   design = ~ condition)
dds_rain <- dds_rain[rowSums(counts(dds_rain)) >= 10,]
# set control condition as reference
dds_rain$condition <- relevel(dds_rain$condition, ref = "LT")
dds_rain <- DESeq(dds_rain)

# see all comparisons (here there is only one)
resultsNames(dds_rain)

res_rain <- results(dds_rain)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
res_rain[order(res_rain$padj),]  
summary(results(dds_rain, alpha=0.05))

padj.cutoff <- 0.005

res_tb_rain <- res_rain %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_rain <- res_tb_rain %>%
  filter(padj < padj.cutoff)

mov10_meta_rain <- coldata_rain %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts_rain <- counts(dds_rain, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

norm_OEsig_rain <- normalized_counts_rain[,c(1:4,5:7)] %>% 
  filter(gene %in% sigOE_rain$gene) 

heat_colors <- inferno(256)

norm_OEsig_rain2 <- norm_OEsig_rain[,-1]
rownames(norm_OEsig_rain2) <- norm_OEsig_rain[,1]

bottom<-pheatmap(norm_OEsig_rain2, 
                 color = heat_colors,
                 cluster_cols = F, 
                 cluster_rows = T, 
                 show_rownames = T,
                 annotation = coldata_rain, 
                 border_color = NA, 
                 fontsize = 10, 
                 scale = "row", 
                 fontsize_row = 10, 
                 height = 20)

plot_grid(as.grob(top),as.grob(nos_top),ncol=1)