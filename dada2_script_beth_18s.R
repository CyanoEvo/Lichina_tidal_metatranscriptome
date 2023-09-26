library(dada2)
setwd("~/Projects/MRes/18S/Beth_18s")

path <- "~/Projects/MRes/18S/Beth_18s/Beth"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(28,32), truncLen=c(130,130),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
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

taxa <- assignTaxonomy(seqtab.nochim, "~/DBs/PR2/pr2_version_4.11.1_dada2.fasta.gz", multithread=TRUE, tryRC=TRUE,taxLevels = c("Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species"))
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
write(asv_fasta, "Euk_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Euk_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "Euk_ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


library("phyloseq")
library(ggplot2)
library(tidyr)
theme_set(theme_bw())
#make a phyloseq object
ASV <- otu_table(Euk_ASVs_counts, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(Euk_ASVs_taxonomy))
sample_info_tab_phy <- sample_data(`18S_meta`)
physeq <- phyloseq(ASV, TAX,sample_info_tab_phy)
physeq1 <- subset_samples(physeq, DATE == '01_04'|DATE=='18_05')
#Ulvales_only <- subset_taxa(physeq, 
#                                  Family=="Ulvales-relatives_X")
physeq1.ord <- ordinate(physeq1, "NMDS", "bray")
plot_ordination(physeq1, physeq1.ord, color="DATE", shape="TIDE")+ geom_point(size = 4)
plot_bar(physeq1, x="SAMPLE_ID", fill="Class") +geom_bar(stat="identity")
Ulvales_only_genus<-tax_glom(Ulvales_only,taxrank = "Genus")
write.csv(tax_table(Ulvales_only),file="Ulvales_tax_all.csv")
write.csv(otu_table(Ulvales_only),file="Ulvales_otu_all.csv")


#physeq1_18df<-psmelt(physeq1)
#physeq1_18_1<-tax_glom(physeq1,taxrank = "Class")
physeq1_18_1<-physeq1
physeq1_18_1 <-  transform_sample_counts(physeq1_18_1, function(x) x / sum(x))
#physeq_1 = filter_taxa(physeq_1, function(x) mean(x) > 0.05, TRUE)
physeq1_18_1df<-psmelt(physeq1_18_1)

# convert Phylum to a character vector from a factor because R
physeq1_18_1df$OTU <- as.character(physeq1_18_1df$OTU)

# group dataframe by Phylum, calculate median rel. abundance
#medians <- ddply(physeq1_18_1df, ~Class, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
#Other <- medians[medians$median <= 0.01,]$Class

# change their name to "Remainder"
#physeq1_18_1df[physeq1_18_1df$Class %in% Other,]$Class<- 'Other'

physeq1_18_1df$OTU[physeq1_18_1df$OTU != 'ASV_1'] <- 'Other'
physeq1_18_1df = aggregate(Abundance ~ SAMPLE_ID + OTU + DATE, data = physeq1_18_1df, FUN = sum, na.rm = TRUE)

# New facet label names for variable
euk_date.labs <- c("Eukaryotes: High tide - Low tide (dry)", "Eukaryotes: High tide - Low tide (rain)")
names(euk_date.labs) <- c("01_04", "18_05")

# Create the plot
euk_abundance<-ggplot(data = physeq1_18_1df, aes(x = SAMPLE_ID,y = Abundance*100)) + geom_bar(colour="black",aes(fill=OTU), stat="identity", position="stack")+ 
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
  scale_fill_manual(values = c("ASV_1"="#ef3b2c","Other"="grey"))


#stats for Ascomycota
Ascomycota_only <- subset_taxa(physeq1_18_1, Class=="Ascomycota")
Ascomycota_only <- tax_glom(Ascomycota_only, taxrank="Class")
my_comparisons <- list( c("HIGH", "LOW"))

Ascomycota_only_01_04 <- subset_samples(Ascomycota_only, SAMPLE_ID == "01_04_HT1A" | SAMPLE_ID == "01_04_HT2A"| SAMPLE_ID == "01_04_HT3A" | SAMPLE_ID == "01_04_LT1A" | SAMPLE_ID == "01_04_LT2A"| SAMPLE_ID == "01_04_LT3A") 
Ascomycota_only_01_04_df<-psmelt(Ascomycota_only_01_04 )
Ascomycota_only_01_04_df$TIDE <- as.factor(Ascomycota_only_01_04_df$TIDE)
Ascomycota_only_01_04_df$TIDE <- factor(Ascomycota_only_01_04_df$TIDE, levels=c("HIGH", "LOW"))
Ascomycota_only_01_04_df$Abundance <- Ascomycota_only_01_04_df$Abundance*100
Ascomycota_01_04p<-ggbarplot(Ascomycota_only_01_04_df, x = "TIDE", y = "Abundance", add = "mean_sd", fill="#ef3b2c")+
  stat_compare_means(method = "t.test",paired = TRUE)+ 
  ggtitle("Ascomycota 01/04/21") + 
  scale_y_continuous(labels=function(x)x*100) +
  theme(aspect.ratio = 2) +
  labs(y="% 18S rRNA transcripts",x="")+ 
  ylim(c(0,100))

Ascomycota_only_18_05 <- subset_samples(Ascomycota_only, SAMPLE_ID == "18_05_HT1A" | SAMPLE_ID == "18_05_HT2A"| SAMPLE_ID == "18_05_HT3A" | SAMPLE_ID == "18_05_LT1A" | SAMPLE_ID == "18_05_LT2A"| SAMPLE_ID == "18_05_LT3A") 
Ascomycota_only_18_05_df<-psmelt(Ascomycota_only_18_05 )
Ascomycota_only_18_05_df$TIDE <- as.factor(Ascomycota_only_18_05_df$TIDE)
Ascomycota_only_18_05_df$TIDE <- factor(Ascomycota_only_18_05_df$TIDE, levels=c("HIGH", "LOW"))
Ascomycota_only_18_05_df$Abundance <- Ascomycota_only_18_05_df$Abundance*100
Ascomycota_18_05p<-ggbarplot(Ascomycota_only_18_05_df, x = "TIDE", y = "Abundance", add = "mean_sd", fill="#ef3b2c")+
  stat_compare_means(method = "t.test",paired = TRUE) +
  ggtitle("Ascomycota 18/05/21 (rain event)") + 
  scale_y_continuous(labels=function(x)x*100)+
  theme(aspect.ratio = 2) + 
  labs(y="% 18S rRNA transcripts",x="")+ 
  ylim(c(0,100))

asco_stats<-plot_grid(Ascomycota_01_04p,Ascomycota_18_05p,ncol=2)


#For supplementary material we would like to show everything in the 18s dataset
physeq <- phyloseq(ASV, TAX,sample_info_tab_phy)
physeq1 <- subset_samples(physeq, DATE == '01_04'|DATE=='18_05')
physeq1_18_ALL <- physeq1
physeq1_18_ALL <- transform_sample_counts(physeq1_18_ALL, function(x) x / sum(x))
physeq1_18_ALLdf <- psmelt(physeq1_18_ALL)

# convert Phylum to a character vector from a factor because R
physeq1_18_ALLdf$Class <- as.character(physeq1_18_ALLdf$Class)
physeq1_18_ALLdf$Class <- physeq1_18_ALLdf$Class %>% replace_na('Unknown')

physeq1_18_ALLdf <- aggregate(Abundance ~ SAMPLE_ID + Class + DATE, data = physeq1_18_ALLdf, FUN = sum, na.rm = TRUE)

# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(physeq1_18_ALLdf, ~Class, function(x) c(median=median(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
Other <- medians[medians$median <= 0.01,]$Class

# change their name to "Remainder"
physeq1_18_ALLdf[physeq1_18_ALLdf$Class %in% Other,]$Class<- 'Other'

# New facet label names for variable
euk_date.labs <- c("Eukaryotes: High tide - Low tide (dry)", "Eukaryotes: High tide - Low tide (rain)")
names(euk_date.labs) <- c("01_04", "18_05")

# Create the plot
euk_abundance<-ggplot(data = physeq1_18_ALLdf, aes(x = SAMPLE_ID,y = Abundance*100)) + geom_bar(colour="black",aes(fill=Class), stat="identity", position="stack")+ 
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