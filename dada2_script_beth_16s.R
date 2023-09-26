#FUNCTION FOR COMPUTING COMPLEX HULL
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
library(dada2)
setwd("~/Projects/MRes/16S/Beth_16s")

path <- "~/Projects/MRes/16S/Beth_16s/fastq"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(28,32), truncLen=c(150,150),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/DBs/Silva/silva_nr_v128_train_set.fa.gz", multithread=TRUE, tryRC=TRUE,taxLevels = c("Kingdom", "Phylum","Class", "Order", "Family", "Genus", "Species"))
taxa <- addSpecies(taxa, "~/DBs/Silva/silva_species_assignment_v128.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "16s_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "16s_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "16s_ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


library("phyloseq")
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)

#Ro's theme ------------------------------------------------------------
theme_ro <- function(){
  theme_bw()+
    theme(
      #legend.position = "none",
          #          panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10))
}
#make a phyloseq object looking at everything
ASV_16 <- otu_table(`16s_ASVs_counts`, taxa_are_rows = TRUE)
TAX_16 <- tax_table(as.matrix(`16s_ASVs_taxonomy`))
sample_info_tab_phy_16 <- sample_data(`16S_meta_new`)
physeq_16 <- phyloseq(ASV_16, TAX_16,sample_info_tab_phy_16)
physeq_16 <- subset_samples(physeq_16, 
                                SAMPLE_ID == "01_04_HT1A" |
                                SAMPLE_ID == "01_04_HT2A" |
                                SAMPLE_ID == "01_04_HT3A" |
                                SAMPLE_ID == "01_04_LT1A" |
                                SAMPLE_ID == "01_04_LT2A" |
                                SAMPLE_ID == "01_04_LT3A" |
                                SAMPLE_ID == "18_05_HT1A" |
                                SAMPLE_ID == "18_05_HT2A" |
                                SAMPLE_ID == "18_05_HT3A" |
                                SAMPLE_ID == "18_05_LT1A" |
                                SAMPLE_ID == "18_05_LT2A" |
                                SAMPLE_ID == "18_05_LT3A")

physeq_16_1 <- subset_taxa(physeq_16,Class!="Chloroplast")
physeq_16_1 <- subset_taxa(physeq_16,Family!="Mitochondria")
#physeq_16_1<-tax_glom(physeq_16_1,taxrank = "Phylum")
physeq_16_1 <-  transform_sample_counts(physeq_16_1, function(x) x / sum(x))

#######Check abundance of non-cyanos###
physeq_16_1_nocyano <- subset_taxa(physeq_16_1,Phylum!="Cyanobacteria")
physeq_16_1_nocyano_df<-psmelt(physeq_16_1_nocyano)
#######################################
#physeq_1 = filter_taxa(physeq_1, function(x) mean(x) > 0.05, TRUE)
physeq_16_1df<-psmelt(physeq_16_1)
physeq_16_1df$Genus<-as.character(physeq_16_1df$Genus)
physeq_16_1df$Genus = replace(x = physeq_16_1df$Genus, list = !physeq_16_1df$Genus %in% c('Rivularia', 'Pleurocapsa'), values =  'Other')
physeq_16_1df = aggregate(Abundance ~ SAMPLE_ID + OTU + Genus +DATE, data = physeq_16_1df, FUN = sum, na.rm = TRUE)
physeq_16_1df<- arrange(physeq_16_1df) %>% mutate(Genus = factor(Genus, levels=c("Rivularia","Pleurocapsa","Other")))

# New facet label names for variable
pro_date.labs <- c("Prokaryotes: High tide - Low tide (dry)", "Prokaryotes: High tide - Low tide (rain)")
names(pro_date.labs) <- c("01_04", "18_05")

pro_abundance<-ggplot(data = physeq_16_1df, aes(x = SAMPLE_ID,y = Abundance*100)) + geom_bar(colour="black",aes(fill=Genus), stat="identity", position="stack")+ 
  theme_bw() + 
  scale_x_discrete(labels= c("HT1A","HT2A","HT3A","LT1A","LT2A","LT3A","HT1A","HT2A","HT3A","LT1A","LT2A","LT3A")) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        aspect.ratio = 0.5) +
  ylab(label = "Total prokaryote rRNA (%)")+
  facet_wrap(facets = "DATE", scales = "free",labeller = labeller(DATE = pro_date.labs))+
  coord_cartesian(xlim=c(1,6), ylim=c(0,100), clip="off") +
  annotate("segment", x = 0.5, xend = 3.5, y = -15, yend = -15,linetype=1,color = "blue")+
  annotate("segment", x = 3.5, xend = 6.5, y = -20, yend = -20,linetype=2,color = "blue")+
  scale_fill_manual(values = c("Rivularia" = "#3690c0","Pleurocapsa"="#016c59","Other"="grey"))

#NEXT subset out the cyanos
cyanobacteria_only <- subset_taxa(physeq_16, 
                             Phylum=="Cyanobacteria")
cyanobacteria_only <- subset_taxa(cyanobacteria_only, 
                                  Class!="Chloroplast")
#write out the taxonomy table. This now needs to be manually updated os make sure the first stage is correct
write.csv(tax_table(cyanobacteria_only), file="cyano_tax.csv")
#read back in the new updated 
asv_tax_new<-cyano_tax_ASV
ASV <- otu_table(`16s_ASVs_counts`, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(asv_tax_new))
sample_info_tab_phy <- sample_data(`16S_meta_new`)
physeq <- phyloseq(ASV, TAX,sample_info_tab_phy)
cyanobacteria_only <- subset_taxa(physeq_16, 
                                  Phylum=="Cyanobacteria")

cyanobacteria_only <-  transform_sample_counts(cyanobacteria_only, function(x) x / sum(x))
cyanobacteria_only = filter_taxa(cyanobacteria_only, function(x) mean(x) > 0.05, TRUE)
cyanobacteria_only<-tax_glom(cyanobacteria_only,taxrank = "Genus")
cyanobacteria_only <- subset_samples(cyanobacteria_only, 
                                       SAMPLE_ID == "01_04_HT1A" |
                                       SAMPLE_ID == "01_04_HT2A" |
                                       SAMPLE_ID == "01_04_HT3A" |
                                       SAMPLE_ID == "01_04_LT1A" |
                                       SAMPLE_ID == "01_04_LT2A" |
                                       SAMPLE_ID == "01_04_LT3A" |
                                       SAMPLE_ID == "18_05_HT1A" |
                                       SAMPLE_ID == "18_05_HT2A" |
                                       SAMPLE_ID == "18_05_HT3A" |
                                       SAMPLE_ID == "18_05_LT1A" |
                                       SAMPLE_ID == "18_05_LT2A" |
                                       SAMPLE_ID == "18_05_LT3A")

ggplot(data = psmelt(cyanobacteria_only), mapping = aes_string(x = "SAMPLE_ID",y = "Abundance")) + geom_bar(colour="black",aes(fill=Genus), stat="identity", position="fill")+ 
  theme_bw() + 
  scale_x_discrete(labels= c("HT1A","HT2A","HT3A","LT1A","LT2A","LT3A","HT1A","HT2A","HT3A","LT1A","LT2A","LT3A")) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
        scale_fill_manual(values=c("#3690c0","#016c59")) +
  facet_wrap(facets = "DATE", scales = "free")

#stats for pleurocapsa
cyanobacteria_only <- subset_taxa(physeq_16,Phylum=="Cyanobacteria")
cyanobacteria_only <-  transform_sample_counts(cyanobacteria_only, function(x) x / sum(x))
#cyanobacteria_only <- filter_taxa(cyanobacteria_only, function(x) mean(x) > 0.05, TRUE)
Pleurocapsa_only <- subset_taxa(cyanobacteria_only, Genus=="Pleurocapsa")
Pleurocapsa_only <- tax_glom(Pleurocapsa_only, taxrank="Genus")
my_comparisons <- list( c("HIGH", "LOW"))

Pleurocapsa_only_01_04 <- subset_samples(Pleurocapsa_only, SAMPLE_ID == "01_04_HT1A" | SAMPLE_ID == "01_04_HT2A"| SAMPLE_ID == "01_04_HT3A" | SAMPLE_ID == "01_04_LT1A" | SAMPLE_ID == "01_04_LT2A"| SAMPLE_ID == "01_04_LT3A") 
Pleurocapsa_only_01_04_df<-psmelt(Pleurocapsa_only_01_04 )
Pleurocapsa_only_01_04_df$TIDE_1 <- as.factor(Pleurocapsa_only_01_04_df$TIDE_1)
Pleurocapsa_only_01_04_df$TIDE_1 <- factor(Pleurocapsa_only_01_04_df$TIDE_1, levels=c("HIGH", "LOW"))
Pleurocapsa_only_01_04_df$Abundance <- Pleurocapsa_only_01_04_df$Abundance*100
Pleurocapsa_01_04p<-ggbarplot(Pleurocapsa_only_01_04_df, x = "TIDE_1", y = "Abundance", add = "mean_sd", fill="#016c59")+
   stat_compare_means(method = "t.test",paired = TRUE)+ 
  ggtitle("Pleurocapsa 01/04/21") + 
  scale_y_continuous(labels=function(x)x*100) +
  theme(aspect.ratio = 2) +
  labs(y="% 16S rRNA transcripts",x="")+ 
  ylim(c(0,100))

Pleurocapsa_only_18_05 <- subset_samples(Pleurocapsa_only, SAMPLE_ID == "18_05_HT1A" | SAMPLE_ID == "18_05_HT2A"| SAMPLE_ID == "18_05_HT3A" | SAMPLE_ID == "18_05_LT1A" | SAMPLE_ID == "18_05_LT2A"| SAMPLE_ID == "18_05_LT3A") 
Pleurocapsa_only_18_05_df<-psmelt(Pleurocapsa_only_18_05 )
Pleurocapsa_only_18_05_df$TIDE_1 <- as.factor(Pleurocapsa_only_18_05_df$TIDE_1)
Pleurocapsa_only_18_05_df$TIDE_1 <- factor(Pleurocapsa_only_18_05_df$TIDE_1, levels=c("HIGH", "LOW"))
Pleurocapsa_only_18_05_df$Abundance <- Pleurocapsa_only_18_05_df$Abundance*100
Pleurocapsa_18_05p<-ggbarplot(Pleurocapsa_only_18_05_df, x = "TIDE_1", y = "Abundance", add = "mean_sd", fill="#016c59")+
  stat_compare_means(method = "t.test",paired = TRUE)+ 
  ggtitle("Pleurocapsa 18/05/21 (rain event)") + 
  scale_y_continuous(labels=function(x)x*100) +
  theme(aspect.ratio = 2) +
  labs(y="% 16S rRNA transcripts",x="")+ 
  ylim(c(0,100))

#stats for Rivularia
cyanobacteria_only <- physeq_16
cyanobacteria_only <-  transform_sample_counts(cyanobacteria_only, function(x) x / sum(x))
#cyanobacteria_only <- filter_taxa(cyanobacteria_only, function(x) mean(x) > 0.05, TRUE)
Rivularia_only <- subset_taxa(cyanobacteria_only, Genus=="Rivularia")
Rivularia_only <- tax_glom(Rivularia_only, taxrank="Genus")
my_comparisons <- list( c("HIGH", "LOW"))

Rivularia_only_01_04 <- subset_samples(Rivularia_only, SAMPLE_ID == "01_04_HT1A" | SAMPLE_ID == "01_04_HT2A"| SAMPLE_ID == "01_04_HT3A" | SAMPLE_ID == "01_04_LT1A" | SAMPLE_ID == "01_04_LT2A"| SAMPLE_ID == "01_04_LT3A") 
Rivularia_only_01_04_df<-psmelt(Rivularia_only_01_04 )
Rivularia_only_01_04_df$TIDE_1 <- as.factor(Rivularia_only_01_04_df$TIDE_1)
Rivularia_only_01_04_df$TIDE_1 <- factor(Rivularia_only_01_04_df$TIDE_1, levels=c("HIGH", "LOW"))
Rivularia_only_01_04_df$Abundance <- Rivularia_only_01_04_df$Abundance*100
Rivularia_01_04p<-ggbarplot(Rivularia_only_01_04_df, x = "TIDE_1", y = "Abundance", add = "mean_sd", fill="#3690c0")+
  stat_compare_means(method = "t.test",paired = TRUE)+ 
  ggtitle("Rivularia 01/04/21") + 
  scale_y_continuous(labels=function(x)x*100) +
  theme(aspect.ratio = 2) +
  labs(y="% 16S rRNA transcripts",x="")+ 
  ylim(c(0,100))

Rivularia_only_18_05 <- subset_samples(Rivularia_only, SAMPLE_ID == "18_05_HT1A" | SAMPLE_ID == "18_05_HT2A"| SAMPLE_ID == "18_05_HT3A" | SAMPLE_ID == "18_05_LT1A" | SAMPLE_ID == "18_05_LT2A"| SAMPLE_ID == "18_05_LT3A") 
Rivularia_only_18_05_df<-psmelt(Rivularia_only_18_05 )
Rivularia_only_18_05_df$TIDE_1 <- as.factor(Rivularia_only_18_05_df$TIDE_1)
Rivularia_only_18_05_df$TIDE_1 <- factor(Rivularia_only_18_05_df$TIDE_1, levels=c("HIGH", "LOW"))
Rivularia_only_18_05_df$Abundance <- Rivularia_only_18_05_df$Abundance*100
Rivularia_18_05p<-ggbarplot(Rivularia_only_18_05_df, x = "TIDE_1", y = "Abundance", add = "mean_sd", fill="#3690c0")+
  stat_compare_means(method = "t.test",paired = TRUE) +
  ggtitle("Rivularia 18/05/21 (rain event)") + 
  scale_y_continuous(labels=function(x)x*100) +
  theme(aspect.ratio = 2) +
  labs(y="% 16S rRNA transcripts",x="")+ 
  ylim(c(0,100))

plot_grid(Pleurocapsa_01_04p,Pleurocapsa_18_05p,
          Rivularia_01_04p,Rivularia_18_05p,ncol=2)


ple_riv_only <- subset_taxa(physeq,Genus=="Pleurocapsa"|
                                         Genus=="Rivularia")
ple_riv_only_18_05 <- subset_samples(ple_riv_only, SAMPLE_ID == "18_05_HT1A" | SAMPLE_ID == "18_05_HT2A"| SAMPLE_ID == "18_05_HT3A" | SAMPLE_ID == "18_05_LT1A" | SAMPLE_ID == "18_05_LT2A"| SAMPLE_ID == "18_05_LT3A") 
ple_riv_only_18_05_df<-psmelt(ple_riv_only_18_05 )
ple_riv_only_18_05_df$TIDE_1 <- as.factor(ple_riv_only_18_05_df$TIDE_1)
ple_riv_only_18_05_df$TIDE_1 <- factor(ple_riv_only_18_05_df$TIDE_1, levels=c("HIGH", "LOW"))
ggbarplot(ple_riv_only_18_05_df, x = "TIDE_1", y = "Abundance", add = "mean_se", fill="skyblue")+
  facet_wrap("Genus")

cyanobacteria_only.ord <- ordinate(cyanobacteria_only, "NMDS", "bray")

plot_ordination(cyanobacteria_only, cyanobacteria_only.ord, color="TIDE")+ 
  geom_point(size = 4, alpha=1)+
  theme_bw()+ 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 1)


# Calculate bray curtis distance matrix
cyanobacteria_bray <- phyloseq::distance(cyanobacteria_only, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(cyanobacteria_only))

# Adonis test
ad<-adonis(cyanobacteria_bray ~ TIDE_1, data = sampledf)
TukeyHSD(sampledf,TIDE_1)



library(vegan)
permanova <- adonis(t(otu) ~ group,
                    data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])


#calculate pairwise distances between thalli 
Thallus_1_HT <- subset_samples(cyanobacteria_only, 
                                SAMPLE_ID=="01_04_HT1A" | 
                                SAMPLE_ID=="27_05_HT1A")
Thallus_1_HT_dist<-distance(Thallus_1_HT,"bray")
Thallus_1_LT <- subset_samples(cyanobacteria_only, 
                                SAMPLE_ID=="01_04_LT1A" |
                                SAMPLE_ID=="27_05_LT1A")
Thallus_1_LT_dist<-distance(Thallus_1_LT,"bray")

Thallus_2_HT <- subset_samples(cyanobacteria_only, 
                               SAMPLE_ID=="01_04_HT2A" | 
                                 SAMPLE_ID=="27_05_HT2A")
Thallus_2_HT_dist<-distance(Thallus_2_HT,"bray")
Thallus_2_LT <- subset_samples(cyanobacteria_only, 
                               SAMPLE_ID=="01_04_LT2A" |
                                 SAMPLE_ID=="27_05_LT2A")
Thallus_3_LT_dist<-distance(Thallus_2_LT,"bray")

Thallus_3_HT <- subset_samples(cyanobacteria_only, 
                               SAMPLE_ID=="01_04_HT3A" | 
                                 SAMPLE_ID=="27_05_HT3A")
Thallus_3_HT_dist<-distance(Thallus_3_HT,"bray")
Thallus_3_LT <- subset_samples(cyanobacteria_only, 
                               SAMPLE_ID=="01_04_LT3A" |
                                 SAMPLE_ID=="27_05_LT3A")
Thallus_3_LT_dist<-distance(Thallus_3_LT,"bray")


Pleurocapsa<- subset_taxa(physeq_16, 
               Genus=="Pleurocapsa")
plot_bar(cyanobacteria_only, x="SAMPLE_ID", fill="Species") +geom_bar(stat="identity",position="stack")




cyanobacteria_only_genus<-tax_glom(cyanobacteria_only,taxrank = "Genus")
write.csv(tax_table(cyanobacteria_only_genus),file="cyano_tax.csv")
write.csv(otu_table(cyanobacteria_only_genus),file="cyano_otu.csv")

dim(tax_table(physeq))
saveRDS(physeq, "Beth_16S_phyloseq.rds")
rowSums(otu_table(physeq))
mean(rowSums(otu_table(physeq)))
min(rowSums(otu_table(physeq))) 
max(rowSums(otu_table(physeq)))  
set.seed(711) 
physeq.rare <- rarefy_even_depth(physeq, sample.size = 167324, trimOTUs = TRUE)  

saveRDS(ps.rare,"Beth_16S_phyloseq_rare.rds")


#For supplementary material we would like to show everything in the 18s dataset
physeq <- phyloseq(ASV, TAX,sample_info_tab_phy)
physeq1 <- subset_samples(physeq, DATE == '01_04'|DATE=='18_05')
physeq1_18_ALL <- physeq1
physeq1_18_ALL <- transform_sample_counts(physeq1_18_ALL, function(x) x / sum(x))
cyanobacteria_only_ALLdf <- psmelt(cyanobacteria_only)

# convert Phylum to a character vector from a factor because R
cyanobacteria_only_ALLdf$Genus <- as.character(cyanobacteria_only_ALLdf$Genus)
cyanobacteria_only_ALLdf$Genus <- cyanobacteria_only_ALLdf$Genus %>% replace_na('Unknown')

cyanobacteria_only_ALLdf <- aggregate(Abundance ~ SAMPLE_ID + Genus + DATE, data = cyanobacteria_only_ALLdf, FUN = sum, na.rm = TRUE)

# group dataframe by Phylum, calculate median rel. abundance
medians_cyano <- ddply(cyanobacteria_only_ALLdf, ~Genus, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
Other_cyano <- medians[medians$median <= 0.01,]$Genus

# change their name to "Remainder"
cyanobacteria_only_ALLdf[cyanobacteria_only_ALLdf$Class %in% Other,]$Genus<- 'Other'

# New facet label names for variable
euk_date.labs <- c("Eukaryotes: High tide - Low tide (dry)", "Eukaryotes: High tide - Low tide (rain)")
names(euk_date.labs) <- c("01_04", "18_05")

# Create the plot
euk_abundance<-ggplot(data = cyanobacteria_only_ALLdf, aes(x = SAMPLE_ID,y = Abundance*100)) + geom_bar(colour="black",aes(fill=Genus), stat="identity", position="stack")+ 
  theme_bw() + 
  scale_x_discrete(labels= c("HT1A","HT2A","HT3A","LT1A","LT2A","LT3A","HT1A","HT2A","HT3A","LT1A","LT2A","LT3A")) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        aspect.ratio = 0.5) +
  ylab(label = "Total eukaryote rRNA (%)")+
  facet_wrap(facets = "DATE", scales = "free",labeller = labeller(DATE = euk_date.labs))+
  coord_cartesian(xlim=c(1,6), ylim=c(0,100), clip="off") +
  annotate("segment", x = 0.5, xend = 3.5, y = -15, yend = -15,linetype=1,color = "blue")+
  annotate("segment", x = 3.5, xend = 6.5, y = -20, yend = -20,linetype=2,color = "blue")+
  scale_fill_manual(values = c("Arthropoda"="#377eb8","Ascomycota"="#e41a1c","Bacillariophyta"="#4daf4a","Discosea-Flabellinia"="#984ea3","Heterolobosea"="#ff7f00","Oligohymenophorea"="#ffff33","Other"="grey","Phyllopharyngea"="#a65628","Unknown"="white"))

euk_abundance