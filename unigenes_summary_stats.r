#this script will count all of the occurences each taxon in the data_new file 
#and output some summary graphs
library(tidyverse)
library(cowplot)
library(RColorBrewer)

col_palette<-brewer.pal(6,"YlGn")

summary_data<-DATA
summary_data$phylum <- as.character(summary_data$phylum) %>% replace_na('Unknown')
summary_data$superkingdom <- as.character(summary_data$superkingdom) %>% replace_na('Unknown')
other_bacteria_taxa<- filter(summary_data, superkingdom =='Bacteria')
other_bacteria_taxa<- filter(other_bacteria_taxa, phylum !='Cyanobacteria')
other_bacteria_taxa<- filter(other_bacteria_taxa, phylum !='Unknown')
other_bacteria_counts<-data.frame(table(other_bacteria_taxa$superkingdom))%>%with(sum(Freq[Var1 == 'Bacteria']))
other_euk_taxa<- filter(summary_data, superkingdom =='Eukaryota')
other_euk_taxa<- filter(other_euk_taxa, phylum !='Ascomycota')
other_euk_taxa<- filter(other_euk_taxa, phylum !='Unknown')
other_euk_counts<-data.frame(table(other_euk_taxa$superkingdom))%>%with(sum(Freq[Var1 == 'Eukaryota']))
archaea_taxa<- filter(summary_data, superkingdom =='Archaea')
archaea_taxa<- filter(summary_data, phylum !='Unknown')
archaea_counts<-data.frame(table(archaea_taxa$superkingdom))%>%with(sum(Freq[Var1 == 'Archaea']))




phyla_counts<-data.frame(table(summary_data$phylum))
main_taxa<- filter(phyla_counts, Freq > 1000)
others_count<-with(phyla_counts, sum(Freq[Var1 != 'Cyanobacteria' & Var1 != 'Ascomycota']))
main_taxa<- add_row(main_taxa, Var1 = "Other bacteria", Freq = other_bacteria_counts)
main_taxa<- add_row(main_taxa, Var1 = "Other eukaryotes", Freq = other_euk_counts)
main_taxa<- add_row(main_taxa, Var1 = "Archaea", Freq = archaea_counts)

phylum_colours<-c("#6a51a3","#e31a1c","#3690c0","#ccece6","#feb24c","#f6eff7")

phylum_plot<-ggplot(main_taxa, aes(x = factor(Var1, level = c('Unknown', 'Archaea',  'Other eukaryotes', 'Other bacteria','Ascomycota', 'Cyanobacteria')), y=Freq)) + 
  geom_bar(aes(fill=Var1), stat = "identity", colour="black") +
  scale_fill_manual(values = phylum_colours) +
  xlab("") +
  ylab("Total Unigenes") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 1) + coord_flip()

cyanos_only<- filter(DATA, phylum == 'Cyanobacteria')
cyanos_only$order <- as.character(cyanos_only$order) %>% replace_na('Unknown')
cyano_order_counts<-data.frame(table(cyanos_only$order))
main_cyanos<- filter(cyano_order_counts, Freq > 100)
cyano_others_count<-with(cyano_order_counts, sum(Freq[Freq < 100]))
main_cyanos<- add_row(main_cyanos, Var1 = "Others", Freq = cyano_others_count)

cyano_colours<-c("#08589e","#2b8cbe","#4eb3d3","#f0f9e8","#a8ddb5","#ccebc5","#7bccc4","grey")

cyano_order_plot<-ggplot(main_cyanos, aes(x=factor(Var1, level = c('Unknown','Others', 'Chroococcales', 'Pseudanabaenales', 'Synechococcales', 'Oscillatoriales', 'Nostocales','Pleurocapsales')), y=Freq)) + 
  geom_bar(aes(fill=Var1), stat = "identity", colour="black") +
  scale_fill_manual(values = cyano_colours) +
  xlab("") +
  ylab("Total Unigenes") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 1)+ coord_flip()

plot_grid(phylum_plot, cyano_order_plot, ncol=1, align="v")





ascomycota_only<- filter(taxa, phylum == 'Ascomycota')
ascomycota_only$order <- ascomycota_only$order %>% replace_na('Unknown')
ascomycota_order_counts<-data.frame(table(ascomycota_only$order))
main_ascomycota<- filter(ascomycota_order_counts, Freq > 100)
ascomycota_others_count<-with(ascomycota_order_counts, sum(Freq[Freq < 100]))
main_ascomycota<- add_row(main_ascomycota, Var1 = "Others", Freq = ascomycota_others_count)

#ascomycota_colours<-c("#08589e","#2b8cbe","#4eb3d3","#f0f9e8","#a8ddb5","#ccebc5","#7bccc4")
factor(Var1, level = c('Unknown','Others', 'Saccharomycetales', 'Mytilinidiales', 'Verrucariales', 'Xylariales', 'Mycosphaerellales', 'Caliciales', 'Glomerales', 'Botryosphaeriales', 'Dothideales', 'Lecanorales', 'Hypocreales', 'Onygenales', 'Pleosporales', 'Eurotiales', 'Umbilicariales', 'Xylonales', 'Helotiales','Geoglossales'))

ascomycota_order_plot<-ggplot(main_ascomycota, aes(x=factor(Var1, level = c('Unknown','Others',   'Verrucariales', 'Saccharomycetales','Xylariales', 'Mytilinidiales','Mycosphaerellales', 'Caliciales', 'Glomerellales', 'Botryosphaeriales', 'Dothideales', 'Lecanorales', 'Hypocreales', 'Onygenales', 'Pleosporales', 'Eurotiales', 'Umbilicariales', 'Xylonales', 'Helotiales','Geoglossales')), y=Freq)) + 
  geom_bar(aes(fill=Var1), stat = "identity", colour="black") +
  #  scale_fill_manual(values = cyano_colours) +
  xlab("Order") +
  ylab("Total Unigenes") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 1)+ coord_flip()




# Libraries
library(packcircles)
library(ggplot2)

# Create data
data <- data.frame(group=paste("Group", letters[1:20]), value=sample(seq(1,100),20)) 

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(data$value, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(data, packing)

# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)

# Make the plot
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size=value, label = group)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

