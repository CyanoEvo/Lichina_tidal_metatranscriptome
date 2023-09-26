#ASCO GENE PLOTS

vec<-c("HT","HT","HT","LT","LT","LT","HT","HT","HT","LT","LT","LT")
vec1<-c("DRY","DRY","DRY","DRY","DRY","DRY","WET","WET","WET","WET","WET","WET")
vec2<-c("HT_DRY","HT_DRY","HT_DRY","LT_DRY","LT_DRY","LT_DRY","HT_WET","HT_WET","HT_WET","LT_WET","LT_WET","LT_WET")

asc_diff_ex_READS_df <- tibble::rownames_to_column(asc_diff_ex_TPM, "NAME")

#ASCOMYCETE GENERAL
asc_ST<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1"|
                                           NAME == "NODE_14384_length_3379_cov_78.120387_g7996_i1.p1"|
                                           NAME == "NODE_4567_length_6455_cov_422.087120_g3060_i0.p1"|
                                           NAME == "NODE_13676_length_3494_cov_536.718211_g9088_i0.p1"|
                                           NAME == "NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1"|
                                           NAME == "NODE_16183_length_3114_cov_396.395594_g8006_i2.p1"|
                                           NAME == "NODE_4962_length_6202_cov_2127.062163_g814_i2.p1"|
                                           NAME == "NODE_1044_length_12085_cov_179.367716_g730_i0.p2"))
                                           #NAME == "NODE_10303_length_4194_cov_265.035671_g6823_i0.p1"|
                                           #NAME == "NODE_16489_length_3071_cov_218.059039_g11031_i0.p1"|
                                           #NAME == "NODE_18046_length_2881_cov_105.940527_g11297_i1.p1"|
                                           #NAME == "NODE_5851_length_5712_cov_344.059940_g3319_i2.p1"|
                                           #NAME == "NODE_9287_length_4466_cov_124.225131_g6160_i0.p1"|
                                           #NAME == "NODE_24063_length_2331_cov_863.463685_g15307_i1.p1"|
                                           #NAME == "NODE_7692_length_4969_cov_416.356209_g5101_i0.p1"|
                                           #NAME == "NODE_6333_length_5497_cov_401.459440_g4201_i0.p1"|
                                           #NAME == "NODE_4815_length_6302_cov_223.018141_g3217_i0.p1"|
                                           #NAME == "NODE_5199_length_6066_cov_338.280828_g3004_i1.p1"|
                                           #NAME == "NODE_6016_length_5641_cov_295.210489_g3993_i0.p1"|
                                           #NAME == "NODE_4787_length_6317_cov_625.977578_g3200_i0.p1"|
                                           #NAME == "NODE_13641_length_3500_cov_130.069157_g9063_i0.p1"|
                                           #NAME == "NODE_7461_length_5052_cov_108.702752_g4955_i0.p1"|
                                           #NAME == "NODE_1092_length_11881_cov_92.157859_g761_i0.p4"|
                                           #NAME == "NODE_32715_length_1841_cov_194.979638_g22800_i1.p1"|
                                           #NAME == "NODE_75155_length_933_cov_183.155814_g1593_i4.p1"|
                                           #NAME == "NODE_29556_length_1994_cov_216.781884_g20722_i0.p1"))

asc_ST$NAME[asc_ST$NAME == "NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1"] <- "ST.1"
asc_ST$NAME[asc_ST$NAME == "NODE_14384_length_3379_cov_78.120387_g7996_i1.p1"] <- "ST.2"
asc_ST$NAME[asc_ST$NAME == "NODE_4567_length_6455_cov_422.087120_g3060_i0.p1"] <- "ST.3"
asc_ST$NAME[asc_ST$NAME == "NODE_13676_length_3494_cov_536.718211_g9088_i0.p1"] <- "ST.4"
asc_ST$NAME[asc_ST$NAME == "NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1"] <- "ST.5"
asc_ST$NAME[asc_ST$NAME == "NODE_16183_length_3114_cov_396.395594_g8006_i2.p1"] <- "ST.6"
asc_ST$NAME[asc_ST$NAME == "NODE_4962_length_6202_cov_2127.062163_g814_i2.p1"] <- "ST.7"
asc_ST$NAME[asc_ST$NAME == "NODE_1044_length_12085_cov_179.367716_g730_i0.p2"] <- "ST.8"
#asc_ST$NAME[asc_ST$NAME == "NODE_10303_length_4194_cov_265.035671_g6823_i0.p1"] <- "UreA"
#asc_ST$NAME[asc_ST$NAME == "NODE_16489_length_3071_cov_218.059039_g11031_i0.p1"] <- "AAP.1"
#asc_ST$NAME[asc_ST$NAME == "NODE_18046_length_2881_cov_105.940527_g11297_i1.p1"] <- "AAP.2"
#asc_ST$NAME[asc_ST$NAME == "NODE_5851_length_5712_cov_344.059940_g3319_i2.p1"] <- "AAP.3"
#asc_ST$NAME[asc_ST$NAME == "NODE_9287_length_4466_cov_124.225131_g6160_i0.p1"] <- "AAP.4"
#asc_ST$NAME[asc_ST$NAME == "NODE_24063_length_2331_cov_863.463685_g15307_i1.p1"] <- "Gln.synt"
#asc_ST$NAME[asc_ST$NAME == "NODE_7692_length_4969_cov_416.356209_g5101_i0.p1"] <- "Glu.synt"
#asc_ST$NAME[asc_ST$NAME == "NODE_6333_length_5497_cov_401.459440_g4201_i0.p1"] <- "ABC.1"
#asc_ST$NAME[asc_ST$NAME == "NODE_4815_length_6302_cov_223.018141_g3217_i0.p1"] <- "ABC.2"
#asc_ST$NAME[asc_ST$NAME == "NODE_5199_length_6066_cov_338.280828_g3004_i1.p1"] <- "ABC.3"
#asc_ST$NAME[asc_ST$NAME == "NODE_6016_length_5641_cov_295.210489_g3993_i0.p1"] <- "ABC.4"
#asc_ST$NAME[asc_ST$NAME == "NODE_4787_length_6317_cov_625.977578_g3200_i0.p1"] <- "ABC.5"
#asc_ST$NAME[asc_ST$NAME == "NODE_13641_length_3500_cov_130.069157_g9063_i0.p1"] <- "ABC.6"
#asc_ST$NAME[asc_ST$NAME == "NODE_7461_length_5052_cov_108.702752_g4955_i0.p1"] <- "ABC.7"
#asc_ST$NAME[asc_ST$NAME == "NODE_1092_length_11881_cov_92.157859_g761_i0.p4"] <- "ABC.8"
#asc_ST$NAME[asc_ST$NAME == "NODE_32715_length_1841_cov_194.979638_g22800_i1.p1"] <- "ABC.9"
#asc_ST$NAME[asc_ST$NAME == "NODE_75155_length_933_cov_183.155814_g1593_i4.p1"] <- "MFS_2.1"
#asc_ST$NAME[asc_ST$NAME == "NODE_29556_length_1994_cov_216.781884_g20722_i0.p1"] <- "MFS_2.2"

asc_ST <- tibble::column_to_rownames(asc_ST, "NAME")
asc_ST_t<-t(asc_ST)
asc_ST_t<-data.frame(asc_ST_t)
asc_ST_t$condition <- vec
asc_ST_t$condition1 <- vec1
asc_ST_t$condition2 <- vec2

asc_ST_tg<-asc_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_ST_tg.summary <- asc_ST_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))
asc_ST_tg.summary_DRY_ONLY<-subset(asc_ST_tg.summary,condition1=="DRY")
asc_ST_plot<-ggplot(asc_ST_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
  #           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  #ylim(0,100) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("ST.1","ST.2","ST.3","ST.4","ST.5","ST.6","ST.7","ST.8"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")+ geom_hline(yintercept=100, lty = 2)

#ASCOMYCETE GH
asc_GH<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_127486_length_619_cov_4086.091575_g110534_i0.p1"|
                                           NAME == "NODE_9354_length_4444_cov_298.459620_g5071_i1.p1"|
                                           NAME == "NODE_8047_length_4854_cov_448.158544_g5339_i0.p1"|
                                           NAME == "NODE_129517_length_612_cov_4088.222635_g112540_i0.p1"|
                                           NAME == "NODE_17039_length_3001_cov_122.702869_g11427_i0.p1"|
                                           NAME == "NODE_34231_length_1775_cov_582.856052_g24353_i0.p1"|
                                           NAME == "NODE_25983_length_2201_cov_313.836466_g17973_i0.p1"|
                                           NAME == "NODE_52713_length_1250_cov_347.120646_g33508_i1.p1"|
                                           NAME == "NODE_24123_length_2327_cov_249.172582_g16557_i0.p1"|
                                           NAME == "NODE_4753_length_6336_cov_116.700942_g3177_i0.p2"|
                                           NAME == "NODE_22485_length_2455_cov_218.863560_g15349_i0.p1"|
                                           NAME == "NODE_9284_length_4468_cov_284.925370_g6157_i0.p1"|
                                           NAME == "NODE_1872_length_9527_cov_163.303998_g1309_i0.p3"|
                                           NAME == "NODE_9878_length_4299_cov_812.515381_g6117_i1.p1"|
                                           NAME == "NODE_15393_length_3223_cov_152.253968_g9526_i1.p1"|
                                           NAME == "NODE_8607_length_4667_cov_254.324771_g5595_i1.p1"|
                                           NAME == "NODE_16424_length_3081_cov_60.676197_g10985_i0.p2"|
                                           NAME == "NODE_16424_length_3081_cov_60.676197_g10985_i0.p3"|
                                           NAME == "NODE_2790_length_8042_cov_367.728071_g1928_i0.p2"|
                                           NAME == "NODE_3735_length_7073_cov_161.156429_g2175_i2.p2" |
                                           NAME == "NODE_14482_length_3364_cov_316.976907_g9109_i2.p1"|
                                           NAME == "NODE_16252_length_3104_cov_415.364236_g10871_i0.p5"|
                                           NAME == "NODE_10636_length_4117_cov_275.776459_g7046_i0.p1"|
                                           NAME == "NODE_13543_length_3518_cov_461.879536_g2557_i2.p1"|
                                           NAME == "NODE_42939_length_1475_cov_555.584165_g26782_i1.p1"|
                                           NAME == "NODE_7801_length_4932_cov_93.028607_g5171_i0.p5"|
                                           NAME == "NODE_8716_length_4632_cov_175.666155_g5787_i0.p1"|
                                           NAME == "NODE_11578_length_3894_cov_50.567391_g1497_i1.p1"|
                                           NAME == "NODE_15866_length_3158_cov_266.793517_g10583_i0.p1"|
                                           NAME == "NODE_5273_length_6024_cov_1611.300118_g3279_i1.p1"|
                                           NAME == "NODE_48399_length_1340_cov_403.876875_g36078_i0.p1"))

asc_GH$NAME[asc_GH$NAME == "NODE_127486_length_619_cov_4086.091575_g110534_i0.p1"] <- "GH17.1"
asc_GH$NAME[asc_GH$NAME == "NODE_9354_length_4444_cov_298.459620_g5071_i1.p1"] <- "GH17.2"
asc_GH$NAME[asc_GH$NAME == "NODE_8047_length_4854_cov_448.158544_g5339_i0.p1"] <- "GH17.3"
asc_GH$NAME[asc_GH$NAME == "NODE_129517_length_612_cov_4088.222635_g112540_i0.p1"] <- "GH17.4"
asc_GH$NAME[asc_GH$NAME == "NODE_17039_length_3001_cov_122.702869_g11427_i0.p1"] <- "GH5.1"
asc_GH$NAME[asc_GH$NAME == "NODE_34231_length_1775_cov_582.856052_g24353_i0.p1"] <- "GH5.2"
asc_GH$NAME[asc_GH$NAME == "NODE_25983_length_2201_cov_313.836466_g17973_i0.p1"] <- "GH16.1"
asc_GH$NAME[asc_GH$NAME == "NODE_52713_length_1250_cov_347.120646_g33508_i1.p1"] <- "GH16.2"
asc_GH$NAME[asc_GH$NAME == "NODE_24123_length_2327_cov_249.172582_g16557_i0.p1"] <- "GH16.3"
asc_GH$NAME[asc_GH$NAME == "NODE_4753_length_6336_cov_116.700942_g3177_i0.p2"] <- "GH16.4"
asc_GH$NAME[asc_GH$NAME == "NODE_22485_length_2455_cov_218.863560_g15349_i0.p1"] <- "GH3.1"
asc_GH$NAME[asc_GH$NAME == "NODE_9284_length_4468_cov_284.925370_g6157_i0.p1"] <- "GH3.2"
asc_GH$NAME[asc_GH$NAME == "NODE_1872_length_9527_cov_163.303998_g1309_i0.p3"] <- "GH3.3"
asc_GH$NAME[asc_GH$NAME == "NODE_9878_length_4299_cov_812.515381_g6117_i1.p1"] <- "GH47.1"
asc_GH$NAME[asc_GH$NAME == "NODE_15393_length_3223_cov_152.253968_g9526_i1.p1"] <- "GH47.2"
asc_GH$NAME[asc_GH$NAME == "NODE_8607_length_4667_cov_254.324771_g5595_i1.p1"] <- "GH47.3"
asc_GH$NAME[asc_GH$NAME == "NODE_16424_length_3081_cov_60.676197_g10985_i0.p2"] <- "GH47.4"
asc_GH$NAME[asc_GH$NAME == "NODE_16424_length_3081_cov_60.676197_g10985_i0.p3"] <- "GH47.5"
asc_GH$NAME[asc_GH$NAME == "NODE_2790_length_8042_cov_367.728071_g1928_i0.p2"] <- "GH18.1"
asc_GH$NAME[asc_GH$NAME == "NODE_3735_length_7073_cov_161.156429_g2175_i2.p2"] <- "GH18.2"
asc_GH$NAME[asc_GH$NAME == "NODE_14482_length_3364_cov_316.976907_g9109_i2.p1"] <- "GH37"
asc_GH$NAME[asc_GH$NAME == "NODE_16252_length_3104_cov_415.364236_g10871_i0.p5"] <- "GH32"
asc_GH$NAME[asc_GH$NAME == "NODE_10636_length_4117_cov_275.776459_g7046_i0.p1"] <- "GH31.1"
asc_GH$NAME[asc_GH$NAME == "NODE_13543_length_3518_cov_461.879536_g2557_i2.p1"] <- "GH31.2"
asc_GH$NAME[asc_GH$NAME == "NODE_42939_length_1475_cov_555.584165_g26782_i1.p1"] <- "GH20.1"
asc_GH$NAME[asc_GH$NAME == "NODE_7801_length_4932_cov_93.028607_g5171_i0.p5"] <- "GH20.2"
asc_GH$NAME[asc_GH$NAME == "NODE_8716_length_4632_cov_175.666155_g5787_i0.p1"] <- "GH63"
asc_GH$NAME[asc_GH$NAME == "NODE_11578_length_3894_cov_50.567391_g1497_i1.p1"] <- "GH81"
asc_GH$NAME[asc_GH$NAME == "NODE_15866_length_3158_cov_266.793517_g10583_i0.p1"] <- "GH128.1"
asc_GH$NAME[asc_GH$NAME == "NODE_5273_length_6024_cov_1611.300118_g3279_i1.p1"] <- "GH65"
asc_GH$NAME[asc_GH$NAME == "NODE_48399_length_1340_cov_403.876875_g36078_i0.p1"] <- "GH115"
#asc_GH$NAME[asc_GH$NAME == "NODE_39896_length_1565_cov_1941.174263_g24942_i1.p1"] <- "HK_tubB"


asc_GH <- tibble::column_to_rownames(asc_GH, "NAME")
asc_GH_t<-t(asc_GH)
asc_GH_t<-data.frame(asc_GH_t)
asc_GH_t$condition <- vec
asc_GH_t$condition1 <- vec1
asc_GH_t$condition2 <- vec2

asc_GH_tg<-asc_GH_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_GH_tg.summary <- asc_GH_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))
asc_GH_tg.summary_DRY_ONLY<-subset(asc_GH_tg.summary,condition1=="DRY")

asc_GH_plot<-ggplot(asc_GH_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
  #           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  #ylim(0,100) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.25) +
  scale_x_discrete(limits = c("GH3.1","GH3.2","GH3.3","","GH5.1","GH5.2","","GH16.1","GH16.2","GH16.3","GH16.4","","GH17.1","GH17.2","GH17.3","GH17.4","","GH18.1","GH18.2","","GH20.1","GH20.2","","GH31.1","GH31.2","","GH32","","GH37","","GH47.1","GH47.2","GH47.3","GH47.4","GH47.5","","GH63","","GH65","","GH81","","GH115","","GH128.1"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")+ geom_hline(yintercept=100, lty = 2)

tbl <- tableGrob(asc_GH_t, rows=NULL)

#ASCOMYCETE Trehalose
asc_tre<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_35370_length_1729_cov_433.846014_g1044_i3.p1"|
                                            NAME == "NODE_48677_length_1334_cov_124.367962_g36314_i0.p1"|
                                            NAME == "NODE_1489_length_10513_cov_660.778640_g1044_i0.p2"|
                                            NAME == "NODE_11285_length_3958_cov_200.838095_g7478_i0.p2"|
                                            NAME == "NODE_1068_length_11988_cov_1377.443559_g746_i0.p3"))



asc_tre$NAME[asc_tre$NAME == "NODE_35370_length_1729_cov_433.846014_g1044_i3.p1"] <- "GT20.1"
asc_tre$NAME[asc_tre$NAME == "NODE_48677_length_1334_cov_124.367962_g36314_i0.p1"] <- "GT20.2"
asc_tre$NAME[asc_tre$NAME == "NODE_1489_length_10513_cov_660.778640_g1044_i0.p2"] <- "GT20.3"
asc_tre$NAME[asc_tre$NAME == "NODE_11285_length_3958_cov_200.838095_g7478_i0.p2"] <- "Trehalose_phosphatase"
asc_tre$NAME[asc_tre$NAME == "NODE_1068_length_11988_cov_1377.443559_g746_i0.p3"] <- "Mannitol_dehydrogenase"

asc_tre <- tibble::column_to_rownames(asc_tre, "NAME")
asc_tre_t<-t(asc_tre)
asc_tre_t<-data.frame(asc_tre_t)
asc_tre_t$condition <- vec
asc_tre_t$condition1 <- vec1
asc_tre_t$condition2 <- vec2

asc_tre_tg<-asc_tre_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_tre_tg.summary <- asc_tre_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))
asc_tre_tg.summary_DRY_ONLY<-subset(asc_tre_tg.summary,condition1=="DRY")

asc_tre_plot<-ggplot(asc_tre_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
  #           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  #ylim(0,100) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("GT20.1","GT20.2","GT20.3","Trehalose_phosphatase","Mannitol_dehydrogenase"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")+ geom_hline(yintercept=100, lty = 2)



library(tidyverse)
library(rstatix)
library(ggpubr)
asc_GH_tg_DRY<-subset(asc_GH_tg, condition1=="DRY")
asc_GH_tg_DRY$logTPM<-log(asc_GH_tg_DRY$TPM+1)
stat.test <- asc_GH_tg_DRY %>%
  group_by(Gene) %>%
  t_test(logTPM ~ condition2, paired=FALSE) %>%
  #  adjust_pvalue(method = "BH") %>%
  add_significance()
stat_test_df<-stat.test
# Create the plot
myplot  <-ggboxplot(
  asc_ST_tg_WET, x = "condition2", y = "logTPM",
  fill = "condition1", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~Gene, ncol=11)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "condition2")
myplot + stat_pvalue_manual(stat.test, label = "p.signif")


