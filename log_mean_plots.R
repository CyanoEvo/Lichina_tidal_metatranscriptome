library(ggplot2)
library(dplyr)
library(cowplot)

dat<-data.frame(DATA)

dat <- dat %>% filter(HT1A_01_04_NumReads != 0)
dat <- dat %>% filter(HT2A_01_04_NumReads != 0)
dat <- dat %>% filter(HT3A_01_04_NumReads != 0)
dat <- dat %>% filter(LT1A_01_04_NumReads != 0)
dat <- dat %>% filter(LT2A_01_04_NumReads != 0)
dat <- dat %>% filter(LT3A_01_04_NumReads != 0)
dat <- dat %>% filter(HT1A_18_05_NumReads != 0)
dat <- dat %>% filter(HT2A_18_05_NumReads != 0)
dat <- dat %>% filter(HT3A_18_05_NumReads != 0)
dat <- dat %>% filter(LT1A_18_05_NumReads != 0)
dat <- dat %>% filter(LT2A_18_05_NumReads != 0)
dat <- dat %>% filter(LT3A_18_05_NumReads != 0)


dat <-mutate(dat,Taxon = factor(Taxon, levels=c("Cyanobacteria", "Ascomycota", "Chlorophyta"))) 
ggplot(dat, aes(x=HT_logMEAN_01_04, y=LT_logMEAN_01_04)) +
  geom_point(aes(colour=Taxon, alpha=0.5)) +
  theme_bw() +
  xlab("")


ggplot(dat, aes(x=(log(HT1A_01_04_TPM+1)+log(HT2A_01_04_TPM+1)+log(HT3A_01_04_TPM+1))/3, y=(log(LT1A_01_04_TPM+1)+log(LT2A_01_04_TPM+1)+log(LT3A_01_04_TPM+1))/3)) + 
  geom_point(aes(colour=order), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
 # scale_colour_manual(values=cols) + 
  theme(legend.position = c(0.2, 0.85))


cols<-c("grey","#3690c0","#016c59")
Cyanos<-subset(dat,order=="Nostocales" | order=="Pleurocapsales")

#Create objects focusing on groups of interest
Nostocs <- as_tibble(dat) %>% mutate(order2 = forcats::fct_other(order, keep = c("Nostocales")))
Nostocs$order2 <- as.character(Nostocs$order2) %>% replace_na('Other')
Nostocs$order2 <- factor(Nostocs$order2, level = c("Other","Nostocales"))

Pleuros <- as_tibble(dat) %>% mutate(order2 = forcats::fct_other(order, keep = c("Pleurocapsales")))
Pleuros$order2 <- as.character(Pleuros$order2) %>% replace_na('Other')
Pleuros$order2 <- factor(Pleuros$order2, level = c("Other","Pleurocapsales"))

Fungi <- as_tibble(dat) %>% mutate(phylum2 = forcats::fct_other(phylum, keep = c("Ascomycota")))
Fungi$phylum2 <- as.character(Fungi$phylum2) %>% replace_na('Other')
Fungi$phylum2 <- factor(Fungi$phylum2, level = c("Other","Ascomycota"))

#plot Dry days
Nostocs_LT_dry<-ggplot(Nostocs %>% arrange(order2), aes(x=(log(HT1A_01_04_TPM+1)+log(HT2A_01_04_TPM+1)+log(HT3A_01_04_TPM+1))/3, y=(log(LT1A_01_04_TPM+1)+log(LT2A_01_04_TPM+1)+log(LT3A_01_04_TPM+1))/3)) + 
  geom_point(aes(colour=order2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("Nostocales" = "#008837",
                                 "Other"="lightgrey")) + 
  theme(legend.position = "NONE") +
  xlab("Nostocales TPM: High-tide (dry)") +
  ylab("Nostocales TPM: Low-tide (dry)")

Pleuros_LT_dry<-ggplot(Pleuros %>% arrange(order2), aes(x=(log(HT1A_01_04_TPM+1)+log(HT2A_01_04_TPM+1)+log(HT3A_01_04_TPM+1))/3, y=(log(LT1A_01_04_TPM+1)+log(LT2A_01_04_TPM+1)+log(LT3A_01_04_TPM+1))/3)) + 
  geom_point(aes(colour=order2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("Pleurocapsales"="#0A50A1",
                                 "Other"="lightgrey")) + 
  theme(legend.position = "NONE") +
  xlab("Pleurocapsales TPM: High-tide (dry)") +
  ylab("Pleurocapsales TPM: Low-tide (dry)")

Fungi_LT_dry<-ggplot(Fungi %>% arrange(phylum2), aes(x=(log(HT1A_01_04_TPM+1)+log(HT2A_01_04_TPM+1)+log(HT3A_01_04_TPM+1))/3, y=(log(LT1A_01_04_TPM+1)+log(LT2A_01_04_TPM+1)+log(LT3A_01_04_TPM+1))/3)) + 
  geom_point(aes(colour=phylum2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values=c("Ascomycota"="#D60066",
                               "Other"="lightgrey"))  + 
  theme(legend.position = "NONE") +
  xlab("Ascomycota TPM: High-tide (dry)") +
  ylab("Ascomycota TPM: Low-tide (dry)")

#plot Wet days
Nostocs_LT_wet<-ggplot(Nostocs %>% arrange(order2), aes(x=(log(HT1A_18_05_TPM+1)+log(HT2A_18_05_TPM+1)+log(HT3A_18_05_TPM+1))/3, y=(log(LT1A_18_05_TPM+1)+log(LT2A_18_05_TPM+1)+log(LT3A_18_05_TPM+1))/3)) + 
  geom_point(aes(colour=order2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("Nostocales" = "#008837",
                                 "Other"="lightgrey")) + 
  theme(legend.position = "NONE") +
  xlab("Nostocales TPM: High-tide (wet)") +
  ylab("Nostocales TPM: Low-tide (wet)")

Pleuros_LT_wet<-ggplot(Pleuros %>% arrange(order2), aes(x=(log(HT1A_18_05_TPM+1)+log(HT2A_18_05_TPM+1)+log(HT3A_18_05_TPM+1))/3, y=(log(LT1A_18_05_TPM+1)+log(LT2A_18_05_TPM+1)+log(LT3A_18_05_TPM+1))/3)) + 
  geom_point(aes(colour=order2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("Pleurocapsales"="#0A50A1",
                                 "Other"="lightgrey")) + 
  theme(legend.position = "NONE") +
  xlab("Pleurocapsales TPM: High-tide (wet)") +
  ylab("Pleurocapsales TPM: Low-tide (wet)")

Fungi_LT_wet<-ggplot(Fungi %>% arrange(phylum2), aes(x=(log(HT1A_18_05_TPM+1)+log(HT2A_18_05_TPM+1)+log(HT3A_18_05_TPM+1))/3, y=(log(LT1A_18_05_TPM+1)+log(LT2A_18_05_TPM+1)+log(LT3A_18_05_TPM+1))/3)) + 
  geom_point(aes(colour=phylum2), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values=c("Ascomycota"="#D60066",
                               "Other"="lightgrey"))  + 
  theme(legend.position = "NONE") +
  xlab("Ascomycota TPM: High-tide (wet)") +
  ylab("Ascomycota TPM: Low-tide (wet)")





plot_grid(Nostocs_LT_dry,Nostocs_LT_wet,Pleuros_LT_dry,Pleuros_LT_wet,Fungi_LT_dry,Fungi_LT_wet,ncol=2)


ggplot(dat, aes(x=HT_logMEAN_18_05, y=LT_logMEAN_18_05)) + 
  geom_point(aes(colour=Taxon, alpha=0.5)) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  facet_wrap(~Taxon)+ geom_abline(intercept = 0, slope = 1, linetype = 2)

ggplot(dat %>%
         arrange(Taxon),
       aes(x = HT_logMEAN_18_05, y = LT_logMEAN_18_05, color = Taxon)) +
  geom_point(alpha=0.3)+ geom_abline(intercept = 0, slope = 1, linetype = 2)



#EXPLORING ASCOMYCETE CLASSES
Fungi<-subset(dat,phylum=="Ascomycota")
ggplot(Fungi, aes(x=(log(HT1A_01_04_TPM+1)+log(HT2A_01_04_TPM+1)+log(HT3A_01_04_TPM+1))/3, y=(log(LT1A_01_04_TPM+1)+log(LT2A_01_04_TPM+1)+log(LT3A_01_04_TPM+1))/3)) + 
  geom_point(aes(colour=order), alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "NONE") +
  facet_wrap(facets = "order")
#############################################################################
Fungi<-subset(dat,phylum=="Ascomycota")

ggplot(ple_diff_ex_TPM, aes(x=(log(HT1A_18_05_ple_TPM+1)+log(HT2A_18_05_ple_TPM+1)+log(HT3A_18_05_ple_TPM+1))/3, y=(log(LT1A_18_05_ple_TPM+1)+log(LT2A_18_05_ple_TPM+1)+log(LT3A_18_05_ple_TPM+1))/3)) + 
  geom_point(alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values="#ef3b2c") + 
  theme(legend.position = c(0.2, 0.85))

ggplot(nos_diff_ex_TPM, aes(x=(log(HT1A_18_05_nos_TPM+1)+log(HT2A_18_05_nos_TPM+1)+log(HT3A_18_05_nos_TPM+1))/3, y=(log(LT1A_18_05_nos_TPM+1)+log(LT2A_18_05_nos_TPM+1)+log(LT3A_18_05_nos_TPM+1))/3)) + 
  geom_point(alpha=0.5) +
  xlim(value= c(0,4)) +
  ylim(value = c(0,4)) +
  theme_classic() +
  theme(aspect.ratio=1) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values="#ef3b2c") + 
  theme(legend.position = c(0.2, 0.85))