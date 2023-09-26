#GET GC for UNIGENES in bash with
# infoseq -auto -only -name -pgc UNIGENES.fasta > UNIGENES_GC

DATA$GC<-UNIGENES_GC$X.GC
DATA_EDIT<-DATA
DATA_EDIT<-DATA_EDIT %>%
  mutate(GC_LABEL = case_when(
    phylum == "Ascomycota"~ "Ascomycota",
    order == "Pleurocapsales"~ "Pleurocapsales",
    order == "Nostocales"~ "Nostocales"
  ))

GC_PLOT_DATA<-data.frame(DATA_EDIT$GC)
GC_PLOT_DATA$GC_LABEL<-DATA_EDIT$GC_LABEL
Gc_p<-na.omit(GC_PLOT_DATA) %>% ggplot(GC_PLOT_DATA, mapping = aes(DATA_EDIT.GC,group=GC_LABEL,fill=GC_LABEL))+
  geom_histogram(position="identity",alpha=0.5,binwidth=0.5)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 aspect.ratio = 0.2) +
  xlab(label="GC content (%)")+
  ylab(label="Number of genes")+
  scale_fill_manual(values = c("Nostocales" = "#3690c0",
                                 "Pleurocapsales"="#016c59",
                                 "Ascomycota"="#ef3b2c"
                                 ))+ scale_x_continuous(breaks=seq(0,100,5))