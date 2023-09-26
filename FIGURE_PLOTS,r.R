#Figure 1
library(cowplot)
#5.5*10
plot_grid(pro_abundance,euk_abundance,ncol=1)


plot_grid(Rivularia_01_04p,Rivularia_18_05p,
          Pleurocapsa_01_04p,Pleurocapsa_18_05p,
          Ascomycota_01_04p,Ascomycota_18_05p,ncol=6)

pp<-plot_grid(phylum_plot, cyano_order_plot, ncol=2, align="v")
Gc_p
lm_p<-plot_grid(Nostocs_LT_dry,Nostocs_LT_wet,Pleuros_LT_dry,Pleuros_LT_wet,Fungi_LT_dry,Fungi_LT_wet,ncol=2)

plot_grid(pp,Gc_p,lm_p, ncol=1)

plot_grid(nos_CS_plot,ple_CS_plot,nos_ST_plot,ple_ST_plot,nos_EPS_plot,ple_EPS_plot,ncol=2, align="hv")


asc_top<-plot_grid(asc_tre_plot,asc_ST_plot, align="h",ncol=2)

asc_GH_plot

plot_grid(asc_top,asc_GH_plot,ncol=1)

