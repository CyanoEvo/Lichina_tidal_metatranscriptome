asco_genes<-subset(DATA, phylum=="Ascomycota")
asco_genes<-as.vector(asco_genes$Name)
write_csv(file = "~/Projects/LICHINA_SALT_TRANSCRIPTOMES/11a_salmon_symnbiont/ascomycota/ascomycota_UNIGENES.list", data.frame(asco_genes))

ple_genes<-subset(DATA, order=="Pleurocapsales")
ple_genes<-as.vector(ple_genes$Name)
write_csv(file = "~/Projects/LICHINA_SALT_TRANSCRIPTOMES/11a_salmon_symnbiont/pleurocapsales/pleurocapsales_UNIGENES.list", data.frame(ple_genes))

nos_genes<-subset(DATA, order=="Nostocales")
nos_genes<-as.vector(nos_genes$Name)
write_csv(file = "~/Projects/LICHINA_SALT_TRANSCRIPTOMES/11a_salmon_symnbiont/nostocales/nostocales_UNIGENES.list", data.frame(nos_genes))
