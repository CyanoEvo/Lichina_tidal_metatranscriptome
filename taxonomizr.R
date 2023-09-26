library(taxonomizr)
library("data.table")
library("readr")
library(ggplot2)
#note this will require a lot of hard drive space, bandwidth and time to process all the data from NCBI
read.names.sql('names.dmp','accessionTaxa.sql')
read.nodes.sql('nodes.dmp','accessionTaxa.sql')
#read.accession2taxid(list.files('.','accession2taxid.gz$'),'accessionTaxa.sql')

taxaId<-data.frame(diamond_out$V4)
taxaId_s<-separate(data = taxaId, col = diamond_out.V4,  into=c("v1","v2","v3","v4","v5","v6","v7","v8","v9"), sep = "\\;")
taxaId_s1<-as.list(taxaId_s$v1)
taxa<-getTaxonomy(taxaId_s1,'accessionTaxa.sql')
print(taxa)
write_csv(file = "~/Projects/LICHINA_SALT_TRANSCRIPTOMES/10a_wrangling/tax_IDS", data.frame(taxa))

taxa<-data.frame(taxa)

phylum_counts <- data.frame(table(taxa$phylum))
threshold<-500
phylum_counts_m<-phylum_counts %>% count(Var1 = fct_collapse(as.character(Var1),Other = unique(as.character(Var1)[Freq < threshold])),wt = Freq, name = "Freq1")
phylum_counts_m$Var1 <- factor(phylum_counts_m$Var1,levels = c("Other","Proteobacteria","Ascomycota","Cyanobacteria"))
ggplot(phylum_counts_m, aes(x=Var1, y=Freq1)) + 
  geom_bar(stat = "identity")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ coord_flip()

cyanobacteria_counts<-subset(taxa,phylum=="Cyanobacteria")
cyanobacteria_order_counts <- as.table(table(cyanobacteria_counts$order))
cyanobacteria_order_prop <- data.frame(prop.table(cyanobacteria_order_counts))
threshold<-500
cyanobacteria_order_prop_m<-cyanobacteria_order_prop %>% count(Var1 = fct_collapse(as.character(Var1),Other = unique(as.character(Var1)[Freq < threshold])),wt = Freq, name = "Freq1")
#cyanobacteria_order_prop_m$Var1 <- factor(phylum_prop_m$Var1,levels = c("Other","Proteobacteria","Ascomycota","Cyanobacteria"))
ggplot(cyanobacteria_order_prop_m, aes(x=Var1, y=Freq1)) + 
  geom_bar(stat = "identity")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ coord_flip()