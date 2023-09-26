library(tidyverse)
library(janitor)

vec<-c("HT","HT","HT","LT","LT","LT","HT","HT","HT","LT","LT","LT")
vec1<-c("DRY","DRY","DRY","DRY","DRY","DRY","WET","WET","WET","WET","WET","WET")
vec2<-c("HT_DRY","HT_DRY","HT_DRY","LT_DRY","LT_DRY","LT_DRY","HT_WET","HT_WET","HT_WET","LT_WET","LT_WET","LT_WET")

nos_diff_ex_READS_df <- tibble::rownames_to_column(nos_diff_ex_TPM, "NAME")
#NOSTOCALES
#HK
#selected from Pinto et al 2012
nos_HK<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"|
                                           NAME == "NODE_889_length_12823_cov_960.285569_g623_i0.p1"|
                                           NAME == "NODE_1912_length_9446_cov_1628.686546_g1338_i0.p1"|
                                           NAME == "NODE_15916_length_3151_cov_3045.708252_g10621_i0.p3"|
                                           NAME == "NODE_106_length_25244_cov_2665.611577_g79_i0.p7"|
                                           NAME == "NODE_185_length_21819_cov_1942.958107_g139_i0.p2"|
                                           NAME == "NODE_5235_length_6040_cov_573.266130_g3493_i0.p2"|
                                           NAME == "NODE_13822_length_3468_cov_15743.773785_g9198_i0.p2"|
                                           NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"|
                                           NAME == "NODE_106234_length_709_cov_355.850629_g89712_i0.p1"))

nos_HK$NAME[nos_HK$NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"] <- "petB"
nos_HK$NAME[nos_HK$NAME == "NODE_889_length_12823_cov_960.285569_g623_i0.p1"] <- "ppc"
nos_HK$NAME[nos_HK$NAME == "NODE_1912_length_9446_cov_1628.686546_g1338_i0.p1"] <- "pyrG"
nos_HK$NAME[nos_HK$NAME == "NODE_15916_length_3151_cov_3045.708252_g10621_i0.p3"] <- "rnpA"
nos_HK$NAME[nos_HK$NAME == "NODE_106_length_25244_cov_2665.611577_g79_i0.p7"] <- "rpoA"
nos_HK$NAME[nos_HK$NAME == "NODE_185_length_21819_cov_1942.958107_g139_i0.p2"] <- "rpoB"
nos_HK$NAME[nos_HK$NAME == "NODE_5235_length_6040_cov_573.266130_g3493_i0.p2"] <- "rps1b.1"
nos_HK$NAME[nos_HK$NAME == "NODE_13822_length_3468_cov_15743.773785_g9198_i0.p2"] <- "rps1b.2"
nos_HK$NAME[nos_HK$NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"] <- "ilvD.1"
nos_HK$NAME[nos_HK$NAME == "NODE_106234_length_709_cov_355.850629_g89712_i0.p1"] <- "ilvD.2"

nos_HK <- tibble::column_to_rownames(nos_HK, "NAME")
nos_HK_t<-t(nos_HK)
nos_HK_t<-data.frame(nos_HK_t)
nos_HK_t$condition <- vec
nos_HK_t$condition1 <- vec1
nos_HK_t$condition2 <- vec2

nos_HK_tg<-nos_HK_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

nos_HK_tg.summary <- nos_HK_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_HK_tg.summary_DRY_ONLY<-subset(nos_HK_tg.summary,condition1=="DRY")
nos_HK_plot<-ggplot(nos_HK_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
 # scale_x_discrete(limits = c("petB","ppc","rps1B","rpoA","ilvD","secA"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Houskeeping (Nostocales)")

#Photosynthesis
nos_PS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p1"|
                                         NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p5"|
                                        NAME == "NODE_81960_length_870_cov_125019.715182_g66510_i0.p1"|
                                        NAME == "NODE_95777_length_768_cov_93695.579856_g79650_i0.p1"|
                                        NAME == "NODE_79543_length_891_cov_12799.143032_g64218_i0.p1"|
                                        #NAME == "NODE_148345_length_556_cov_34865.022774_g131151_i0.p1"|
                                        #NAME == "NODE_42184_length_1497_cov_13336.356742_g27908_i1.p2"|
                                          NAME == "NODE_7966_length_4882_cov_15047.541485_g5282_i0.p5"|
                                          NAME == "NODE_2677_length_8191_cov_2194.075758_g111_i7.p4"|
                                          NAME == "NODE_3079_length_7723_cov_9104.390327_g2104_i0.p2"))

nos_PS$NAME[nos_PS$NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p1"] <- "psaA.1"
nos_PS$NAME[nos_PS$NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p5"] <- "psaA.2"
nos_PS$NAME[nos_PS$NAME == "NODE_81960_length_870_cov_125019.715182_g66510_i0.p1"] <- "psbA.1"
nos_PS$NAME[nos_PS$NAME == "NODE_95777_length_768_cov_93695.579856_g79650_i0.p1"] <- "psbA.2"
nos_PS$NAME[nos_PS$NAME == "NODE_79543_length_891_cov_12799.143032_g64218_i0.p1"] <- "psbA.3"
#nos_PS$NAME[nos_PS$NAME == "NODE_148345_length_556_cov_34865.022774_g131151_i0.p1"] <- "PsbA.4"
#nos_PS$NAME[nos_PS$NAME == "NODE_42184_length_1497_cov_13336.356742_g27908_i1.p2"] <- "PsbA.5"
nos_PS$NAME[nos_PS$NAME == "NODE_7966_length_4882_cov_15047.541485_g5282_i0.p5"] <- "cpcA.1"
nos_PS$NAME[nos_PS$NAME == "NODE_2677_length_8191_cov_2194.075758_g111_i7.p4"] <- "cpcA.2"
nos_PS$NAME[nos_PS$NAME == "NODE_3079_length_7723_cov_9104.390327_g2104_i0.p2"] <- "rbcL"


nos_PS <- tibble::column_to_rownames(nos_PS, "NAME")
nos_PS_t<-t(nos_PS)
nos_PS_t<-data.frame(nos_PS_t)
nos_PS_t$condition <- vec
nos_PS_t$condition1 <- vec1
nos_PS_t$condition2 <- vec2

nos_PS_tg<-nos_PS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

nos_PS_tg.summary <- nos_PS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_PS_tg.summary_DRY_ONLY<-subset(nos_PS_tg.summary,condition1=="DRY")

nos_PS_plot<-ggplot(nos_PS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("psaA.1","psaA.2","psbA.1","psbA.2","psbA.3","cpcA.1","cpcA.2","rbcL"))+
   xlab("") +
  ylab("Nostocales transcripts per million (TPM)")+
  ggtitle("Photosynthesis (Nostocales)")+ geom_hline(yintercept=20000, lty = 2)

#Compatible solutes
nos_CS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_12033_length_3790_cov_118.915254_g7984_i0.p3"|
                                           NAME == "NODE_2612_length_8277_cov_218.485617_g1802_i0.p3"|
                                           NAME == "NODE_31055_length_1917_cov_838.574295_g21875_i0.p2"|
                                           NAME == "NODE_10804_length_4072_cov_1164.011003_g7147_i0.p1"|
                                           NAME == "NODE_9636_length_4362_cov_924.362789_g6391_i0.p2"|
                                           NAME == "NODE_14416_length_3374_cov_112.850954_g9589_i0.p1"|
                                           NAME == "NODE_8127_length_4825_cov_879.277778_g5394_i0.p1"|
                                           NAME == "NODE_2706_length_8156_cov_549.081777_g1867_i0.p1"|
                                           NAME == "NODE_172604_length_501_cov_401.943925_g155254_i0.p1"|
                                           NAME == "NODE_252277_length_393_cov_351.981250_g234788_i0.p1"|
                                           NAME == "NODE_140_length_23980_cov_1386.421634_g107_i0.p6"))

nos_CS$NAME[nos_CS$NAME == "NODE_12033_length_3790_cov_118.915254_g7984_i0.p3"] <- "spsA.1"
nos_CS$NAME[nos_CS$NAME == "NODE_2612_length_8277_cov_218.485617_g1802_i0.p3"] <- "spsA.2"
nos_CS$NAME[nos_CS$NAME == "NODE_31055_length_1917_cov_838.574295_g21875_i0.p2"] <- "spp"
nos_CS$NAME[nos_CS$NAME == "NODE_10804_length_4072_cov_1164.011003_g7147_i0.p1"] <- "treS.1"
nos_CS$NAME[nos_CS$NAME == "NODE_9636_length_4362_cov_924.362789_g6391_i0.p2"] <- "treS.2"
nos_CS$NAME[nos_CS$NAME == "NODE_14416_length_3374_cov_112.850954_g9589_i0.p1"] <- "treS.3"
nos_CS$NAME[nos_CS$NAME == "NODE_8127_length_4825_cov_879.277778_g5394_i0.p1"] <- "treY"
nos_CS$NAME[nos_CS$NAME == "NODE_2706_length_8156_cov_549.081777_g1867_i0.p1"] <- "treZ"
nos_CS$NAME[nos_CS$NAME == "NODE_172604_length_501_cov_401.943925_g155254_i0.p1"] <- "mtdh.1"
nos_CS$NAME[nos_CS$NAME == "NODE_252277_length_393_cov_351.981250_g234788_i0.p1"] <- "mtdh.2"
nos_CS$NAME[nos_CS$NAME == "NODE_140_length_23980_cov_1386.421634_g107_i0.p6"] <- "DDG"

nos_CS <- tibble::column_to_rownames(nos_CS, "NAME")
nos_CS_t<-t(nos_CS)
nos_CS_t<-data.frame(nos_CS_t)
nos_CS_t$condition <- vec
nos_CS_t$condition1 <- vec1
nos_CS_t$condition2 <- vec2

nos_CS_tg<-nos_CS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
nos_CS_tg$logTPM<-log(nos_CS_tg$TPM)

nos_CS_tg.summary <- nos_CS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_CS_tg.summary_DRY_ONLY<-subset(nos_CS_tg.summary,condition1=="DRY")
nos_CS_plot<-ggplot(nos_CS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("spsA.1","spsA.2","spp","treS.1","treS.2","treS.3","treY","treZ","mtdh.1","mtdh.2","DDG"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Compatible solutes (Nostocales)")+ geom_hline(yintercept=100, lty = 2)

#SUGAR TRANSPORT
nos_ST<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_1164_length_11617_cov_1128.037509_g816_i0.p3"|
                                           NAME == "NODE_54162_length_1222_cov_147.726719_g41060_i0.p1"|
                                           NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p14"|
                                          # NAME == "NODE_15345_length_3231_cov_364.330906_g10223_i0.p1"|
                                          # NAME == "NODE_2840_length_7989_cov_439.225366_g1958_i0.p2"|
                                          # NAME == "NODE_9347_length_4447_cov_83.253086_g6203_i0.p3"|
                                          # NAME == "NODE_57884_length_1157_cov_325.278598_g44336_i0.p1"|
                                          # NAME == "NODE_17864_length_2902_cov_62.169671_g12004_i0.p1"|
                                          # NAME == "NODE_2511_length_8405_cov_151.375060_g1733_i0.p5"|
                                          # NAME == "NODE_4033_length_6854_cov_257.324436_g2719_i0.p2"|
                                          # NAME == "NODE_99632_length_744_cov_163.602086_g83347_i0.p1"|
                                          # NAME == "NODE_80125_length_886_cov_48.207872_g64770_i0.p1"|
                                          # NAME == "NODE_6945_length_5236_cov_517.852411_g4620_i0.p2"|
                                          # NAME == "NODE_88983_length_814_cov_430.758435_g73175_i0.p1"|
                                          # NAME == "NODE_1206_length_11448_cov_120.996308_g848_i0.p8"|
                                          # NAME == "NODE_311_length_18330_cov_105.961659_g224_i0.p11"|
                                           NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p16"|
                                           NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p8"|
                                           NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p6"|
                                           NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p1"))


nos_ST$NAME[nos_ST$NAME == "NODE_1164_length_11617_cov_1128.037509_g816_i0.p3"] <- "alr4781.1"
nos_ST$NAME[nos_ST$NAME == "NODE_54162_length_1222_cov_147.726719_g41060_i0.p1"] <- "alr4781.2"
nos_ST$NAME[nos_ST$NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p14"] <- "all0261"
#nos_ST$NAME[nos_ST$NAME == "NODE_15345_length_3231_cov_364.330906_g10223_i0.p1"] <- "all1823.1"
#nos_ST$NAME[nos_ST$NAME == "NODE_2840_length_7989_cov_439.225366_g1958_i0.p2"] <- "all1823.2"
#nos_ST$NAME[nos_ST$NAME == "NODE_9347_length_4447_cov_83.253086_g6203_i0.p3"] <- "all1823.3"
#nos_ST$NAME[nos_ST$NAME == "NODE_57884_length_1157_cov_325.278598_g44336_i0.p1"] <- "all1823.4"
#nos_ST$NAME[nos_ST$NAME == "NODE_17864_length_2902_cov_62.169671_g12004_i0.p1"] <- "all1823.5"
#nos_ST$NAME[nos_ST$NAME == "NODE_2511_length_8405_cov_151.375060_g1733_i0.p5"] <- "all1823.6"
#nos_ST$NAME[nos_ST$NAME == "NODE_4033_length_6854_cov_257.324436_g2719_i0.p2"] <- "all1823.7"
#nos_ST$NAME[nos_ST$NAME == "NODE_99632_length_744_cov_163.602086_g83347_i0.p1"] <- "all1823.8"
#nos_ST$NAME[nos_ST$NAME == "NODE_80125_length_886_cov_48.207872_g64770_i0.p1"] <- "all1823.9"
#nos_ST$NAME[nos_ST$NAME == "NODE_6945_length_5236_cov_517.852411_g4620_i0.p2"] <- "all1823.10"
#nos_ST$NAME[nos_ST$NAME == "NODE_88983_length_814_cov_430.758435_g73175_i0.p1"] <- "all1823.11"
#nos_ST$NAME[nos_ST$NAME == "NODE_1206_length_11448_cov_120.996308_g848_i0.p8"] <- "all1823.12"
#nos_ST$NAME[nos_ST$NAME == "NODE_311_length_18330_cov_105.961659_g224_i0.p11"] <- "all1823.13"
nos_ST$NAME[nos_ST$NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p16"] <- "alr2532.1"
nos_ST$NAME[nos_ST$NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p8"] <- "alr2532.2"
nos_ST$NAME[nos_ST$NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p6"] <- "alr2532.3"
nos_ST$NAME[nos_ST$NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p1"] <- "alr2532.4"

nos_ST <- tibble::column_to_rownames(nos_ST, "NAME")
nos_ST_t<-t(nos_ST)
nos_ST_t<-data.frame(nos_ST_t)
nos_ST_t$condition <- vec
nos_ST_t$condition1 <- vec1
nos_ST_t$condition2 <- vec2

nos_ST_tg<-nos_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

nos_ST_tg.summary <- nos_ST_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_ST_tg.summary_DRY_ONLY<-subset(nos_ST_tg.summary,condition1=="DRY")
nos_ST_plot<-ggplot(nos_ST_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("alr4781.1","alr4781.2","all0261","alr2532.1","alr2532.2","alr2532.3","alr2532.4"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Sugar transport (Nostocales)")+ geom_hline(yintercept=100, lty = 2)

#EPS
nos_EPS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p1"|
                                          NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p2"|
                                          NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p3"|
                                          NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p5"|
                                          NAME == "NODE_441_length_16193_cov_75.542184_g323_i0.p3"|
                                          NAME == "NODE_14243_length_3402_cov_434.991589_g9475_i0.p1"|
                                          NAME == "NODE_13298_length_3560_cov_641.725839_g8837_i0.p1"|
                                          NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p1"|
                                          NAME == "NODE_65157_length_1047_cov_216.448665_g50892_i0.p1"|
                                          NAME == "NODE_17_length_36777_cov_40.073507_g12_i0.p7"|
                                          NAME == "NODE_872_length_12899_cov_410.938562_g612_i0.p3"|
                                          NAME == "NODE_11762_length_3853_cov_285.506085_g7799_i0.p4"|
                                          NAME == "NODE_11724_length_3862_cov_129.924518_g7773_i0.p3"|
                                          NAME == "NODE_138615_length_583_cov_151.256863_g121521_i0.p1"|
                                          NAME == "NODE_4220_length_6703_cov_57.689894_g2834_i0.p6"|
                                          NAME == "NODE_1851_length_9582_cov_58.061521_g1292_i0.p1"))



nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p1"] <- "wzc"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p2"] <- "wzy"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p3"] <- "wzx.1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p5"] <- "wzx.2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_441_length_16193_cov_75.542184_g323_i0.p3"] <- "wzx.3"
nos_EPS$NAME[nos_EPS$NAME == "NODE_14243_length_3402_cov_434.991589_g9475_i0.p1"] <- "kpsD"
nos_EPS$NAME[nos_EPS$NAME == "NODE_13298_length_3560_cov_641.725839_g8837_i0.p1"] <- "kpsE.1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p1"] <- "kpsE.2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_65157_length_1047_cov_216.448665_g50892_i0.p1"] <- "kpsT.1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_17_length_36777_cov_40.073507_g12_i0.p7"] <- "kpsT.2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_872_length_12899_cov_410.938562_g612_i0.p3"] <- "kpsT.3"
nos_EPS$NAME[nos_EPS$NAME == "NODE_11762_length_3853_cov_285.506085_g7799_i0.p4"] <- "kpsT.4"
nos_EPS$NAME[nos_EPS$NAME == "NODE_11724_length_3862_cov_129.924518_g7773_i0.p3"] <- "kpsT.5"
nos_EPS$NAME[nos_EPS$NAME == "NODE_138615_length_583_cov_151.256863_g121521_i0.p1"] <- "kpsM.1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4220_length_6703_cov_57.689894_g2834_i0.p6"] <- "kpsM.2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_1851_length_9582_cov_58.061521_g1292_i0.p1"] <- "bcsA"
nos_EPS <- tibble::column_to_rownames(nos_EPS, "NAME")
nos_EPS_t<-t(nos_EPS)
nos_EPS_t<-data.frame(nos_EPS_t)
nos_EPS_t$condition <- vec
nos_EPS_t$condition1 <- vec1
nos_EPS_t$condition2 <- vec2

nos_EPS_tg<-nos_EPS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

nos_EPS_tg.summary <- nos_EPS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_EPS_tg.summary_DRY_ONLY<-subset(nos_EPS_tg.summary,condition1=="DRY")
nos_EPS_plot<-ggplot(nos_EPS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("wzc","wzy","wzx.1","wzx.2","wzx.3", "kpsD","kpsE.1","kpsE.2","kpsT.1","kpsT.2","kpsT.3","kpsT.4","kpsT.5","kpsM.1","kpsM.2","bcsA"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("EPS export (Nostocales)")+ geom_hline(yintercept=100, lty = 2)





#Nitrogen
nos_NT<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_10694_length_4102_cov_4262.831472_g7079_i0.p3"|
                                         NAME == "NODE_20304_length_2646_cov_2351.947532_g9211_i1.p1"|
                                         NAME == "NODE_13842_length_3463_cov_1751.647493_g9211_i0.p2"|
                                         NAME == "NODE_117058_length_659_cov_6419.416382_g100276_i0.p1"|
                                         NAME == "NODE_13842_length_3463_cov_1751.647493_g9211_i0.p4"|
                                         NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p2"|
                                         NAME == "NODE_3972_length_6895_cov_886.344767_g2681_i0.p1"|
                                         NAME == "NODE_361_length_17230_cov_65.716559_g264_i0.p5"|
                                         NAME == "NODE_917_length_12716_cov_949.146089_g643_i0.p3"))

nos_NT$NAME[nos_NT$NAME == "NODE_10694_length_4102_cov_4262.831472_g7079_i0.p3"] <- "NifH"
nos_NT$NAME[nos_NT$NAME == "NODE_20304_length_2646_cov_2351.947532_g9211_i1.p1"] <- "NifD"
nos_NT$NAME[nos_NT$NAME == "NODE_13842_length_3463_cov_1751.647493_g9211_i0.p2"] <- "NifK.1"
nos_NT$NAME[nos_NT$NAME == "NODE_117058_length_659_cov_6419.416382_g100276_i0.p1"] <- "NifK.2"
nos_NT$NAME[nos_NT$NAME == "NODE_13842_length_3463_cov_1751.647493_g9211_i0.p4"] <- "NifK.3"
nos_NT$NAME[nos_NT$NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p2"] <- "HetR"
nos_NT$NAME[nos_NT$NAME == "NODE_3972_length_6895_cov_886.344767_g2681_i0.p1"] <- "alr0991"
nos_NT$NAME[nos_NT$NAME == "NODE_361_length_17230_cov_65.716559_g264_i0.p5"] <- "alr0992"
nos_NT$NAME[nos_NT$NAME == "NODE_917_length_12716_cov_949.146089_g643_i0.p3"] <- "NrtP"

nos_NT <- tibble::column_to_rownames(nos_NT, "NAME")
nos_NT_t<-t(nos_NT)
nos_NT_t<-data.frame(nos_NT_t)
nos_NT_t$condition <- vec
nos_NT_t$condition1 <- vec1
nos_NT_t$condition2 <- vec2

nos_NT_tg<-nos_NT_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
nos_NT_tg$logTPM<-log(nos_NT_tg$TPM)

nos_NT_tg.summary <- nos_NT_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

nos_NT_tg.summary_DRY_ONLY<-subset(nos_NT_tg.summary,condition1=="DRY")
nos_NT_plot<-ggplot(nos_NT_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#37ABC8", "#005544")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("NifH","NifD","NifK.1","NifK.2","NifK.3","HetR","alr0991","alr0992","NrtP"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Nitrogen (Nostocales)")+ geom_hline(yintercept=100, lty = 2)

plot_grid(nos_PS_plot,ple_PS_plot,nos_CS_plot,ple_CS_plot,nos_ST_plot,ple_ST_plot,nos_EPS_plot,ple_EPS_plot,ncol=2, align="hv")

#PLeurocapsales
#HK
#selected from Pinto et al 2012
ple_HK<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_1856_length_9577_cov_102.255787_g1297_i0.p1"|
                                         NAME == "NODE_461_length_15967_cov_59.193784_g335_i0.p7"|
                                         NAME == "NODE_864_length_12967_cov_50.308128_g605_i0.p1"|
                                         NAME == "NODE_1208_length_11446_cov_85.420206_g850_i0.p1"|
                                         NAME == "NODE_105_length_25260_cov_67.393100_g78_i0.p1"|
                                         NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p1"|
                                         NAME == "NODE_1640_length_10104_cov_99.995813_g1142_i0.p2"|
                                         NAME == "NODE_12921_length_3621_cov_55.577508_g8558_i0.p1"|
                                         NAME == "NODE_206_length_21172_cov_92.048960_g154_i0.p7"|
                                         NAME == "NODE_378_length_16938_cov_78.142425_g276_i0.p4"|
                                         NAME == "NODE_119_length_24550_cov_65.461454_g91_i0.p2"|
                                         NAME == "NODE_46349_length_1390_cov_5.285497_g34316_i0.p1"|
                                         NAME == "NODE_181181_length_485_cov_3.973301_g163799_i0.p1"|
                                         NAME == "NODE_97_length_25645_cov_423.605702_g72_i0.p6"|
                                         NAME == "NODE_728_length_13753_cov_53.621345_g422_i1.p4"|
                                         NAME == "NODE_46375_length_1389_cov_16.399696_g34339_i0.p1"|
                                         NAME == "NODE_21_length_35495_cov_1118.934645_g15_i0.p10"|
                                         NAME == "NODE_9688_length_4350_cov_54.675707_g6431_i0.p1"|
                                         NAME == "NODE_1914_length_9439_cov_267.062460_g772_i1.p2"|
                                         NAME == "NODE_1856_length_9577_cov_102.255787_g1297_i0.p1"|
                                         NAME == "NODE_16996_length_3006_cov_14.371292_g11065_i1.p1"|
                                         NAME == "NODE_145898_length_563_cov_3.181633_g128723_i0.p1"))

ple_HK$NAME[ple_HK$NAME == "NODE_1856_length_9577_cov_102.255787_g1297_i0.p1"] <- "ilvD"
ple_HK$NAME[ple_HK$NAME == "NODE_461_length_15967_cov_59.193784_g335_i0.p7"] <- "petB"
ple_HK$NAME[ple_HK$NAME == "NODE_864_length_12967_cov_50.308128_g605_i0.p1"] <- "ppC.1"
ple_HK$NAME[ple_HK$NAME == "NODE_1208_length_11446_cov_85.420206_g850_i0.p1"] <- "ppC.2"
ple_HK$NAME[ple_HK$NAME == "NODE_105_length_25260_cov_67.393100_g78_i0.p1"] <- "ppC.3"
ple_HK$NAME[ple_HK$NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p1"] <- "ppC.4"
ple_HK$NAME[ple_HK$NAME == "NODE_1640_length_10104_cov_99.995813_g1142_i0.p2"] <- "pyrG"
ple_HK$NAME[ple_HK$NAME == "NODE_12921_length_3621_cov_55.577508_g8558_i0.p1"] <- "rpoA1"
ple_HK$NAME[ple_HK$NAME == "NODE_206_length_21172_cov_92.048960_g154_i0.p7"] <- "rpoA2"
ple_HK$NAME[ple_HK$NAME == "NODE_378_length_16938_cov_78.142425_g276_i0.p4"] <- "rpoA3"
ple_HK$NAME[ple_HK$NAME == "NODE_119_length_24550_cov_65.461454_g91_i0.p2"] <- "rpoB1"
ple_HK$NAME[ple_HK$NAME == "NODE_46349_length_1390_cov_5.285497_g34316_i0.p1"] <- "rpoB2"
ple_HK$NAME[ple_HK$NAME == "NODE_181181_length_485_cov_3.973301_g163799_i0.p1"] <- "rpoB3"
ple_HK$NAME[ple_HK$NAME == "NODE_97_length_25645_cov_423.605702_g72_i0.p6"] <- "rps1B.1"
ple_HK$NAME[ple_HK$NAME == "NODE_728_length_13753_cov_53.621345_g422_i1.p4"] <- "rps1B.2"
ple_HK$NAME[ple_HK$NAME == "NODE_46375_length_1389_cov_16.399696_g34339_i0.p1"] <- "rps1B.3"
ple_HK$NAME[ple_HK$NAME == "NODE_21_length_35495_cov_1118.934645_g15_i0.p10"] <- "rps1B.4"
ple_HK$NAME[ple_HK$NAME == "NODE_9688_length_4350_cov_54.675707_g6431_i0.p1"] <- "rps1B.5"
ple_HK$NAME[ple_HK$NAME == "NODE_1914_length_9439_cov_267.062460_g772_i1.p2"] <- "rps1B.6"
ple_HK$NAME[ple_HK$NAME == "NODE_1856_length_9577_cov_102.255787_g1297_i0.p1"] <- "ilvD"
ple_HK$NAME[ple_HK$NAME == "NODE_16996_length_3006_cov_14.371292_g11065_i1.p1"] <- "secA.1"
ple_HK$NAME[ple_HK$NAME == "NODE_145898_length_563_cov_3.181633_g128723_i0.p1"] <- "secA.2"

ple_HK <- tibble::column_to_rownames(ple_HK, "NAME")
ple_HK_t<-t(ple_HK)
ple_HK_t<-data.frame(ple_HK_t)
ple_HK_t$condition <- vec
ple_HK_t$condition1 <- vec1
ple_HK_t$condition2 <- vec2

ple_HK_tg<-ple_HK_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_HK_tg.summary <- ple_HK_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_HK_tg.summary_DRY_ONLY<-subset(ple_HK_tg.summary,condition1=="DRY")
ple_HK_plot<-ggplot(ple_HK_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#37ABC8", "#005544")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("ilvD","PetB"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Houskeeping (Pleurocapsales)")

#PLeurocapsales
#PS
ple_PS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_19751_length_2701_cov_4954.788813_g13333_i0.p1"|
                                         NAME == "NODE_1511_length_10440_cov_4349.966818_g1058_i0.p1"|
                                         NAME == "NODE_1375_length_10838_cov_1438.095866_g969_i0.p1"|
                                         NAME == "NODE_128818_length_615_cov_3.321033_g111851_i0.p1"|
                                         NAME == "NODE_372759_length_317_cov_1.176230_g355263_i0.p1"|
                                         NAME == "NODE_70_length_27121_cov_248.993382_g53_i0.p6"|
                                         NAME == "NODE_23791_length_2350_cov_26283.451910_g802_i3.p1"|
                                         NAME == "NODE_174_length_22195_cov_217.718380_g132_i0.p8"|
                                         NAME == "NODE_2433_length_8539_cov_9166.021970_g1679_i0.p2"|
                                         NAME == "NODE_1935_length_9397_cov_2067.955062_g657_i2.p3"|
                                         NAME == "NODE_103641_length_722_cov_870.023112_g87213_i0.p1"|
                                         NAME == "NODE_41_length_30295_cov_377.508338_g29_i0.p20"|
                                         NAME == "NODE_2632_length_8241_cov_291.273384_g1815_i0.p7"|
                                         NAME == "NODE_14124_length_3420_cov_96.905886_g9397_i0.p4"|
                                         NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p5"
                                         #NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p6"|
                                         #NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p19"|
                                         #NAME == "NODE_6321_length_5503_cov_104.721179_g4191_i0.p8"|
                                         #NAME == "NODE_7497_length_5038_cov_463.560524_g4982_i0.p1"|
                                         #NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p1"|
                                         #NAME == "NODE_7638_length_4993_cov_64.689228_g5068_i0.p1"|
                                         #NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p1"|
                                         #NAME == "NODE_6744_length_5315_cov_351.019649_g464_i2.p1"|
                                         #NAME == "NODE_6967_length_5229_cov_37.146431_g4191_i1.p1"|
                                         #NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p2"
                                         ))

ple_PS$NAME[ple_PS$NAME == "NODE_19751_length_2701_cov_4954.788813_g13333_i0.p1"] <- "psaA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_1511_length_10440_cov_4349.966818_g1058_i0.p1"] <- "psaA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_1375_length_10838_cov_1438.095866_g969_i0.p1"] <- "psaA.3"
ple_PS$NAME[ple_PS$NAME == "NODE_128818_length_615_cov_3.321033_g111851_i0.p1"] <- "psaA.4"
ple_PS$NAME[ple_PS$NAME == "NODE_372759_length_317_cov_1.176230_g355263_i0.p1"] <- "psaA.5"
ple_PS$NAME[ple_PS$NAME == "NODE_70_length_27121_cov_248.993382_g53_i0.p6"] <- "psbA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_23791_length_2350_cov_26283.451910_g802_i3.p1"] <- "psbA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_174_length_22195_cov_217.718380_g132_i0.p8"] <- "psbA.3"
ple_PS$NAME[ple_PS$NAME == "NODE_2433_length_8539_cov_9166.021970_g1679_i0.p2"] <- "psbA.4"
ple_PS$NAME[ple_PS$NAME == "NODE_1935_length_9397_cov_2067.955062_g657_i2.p3"] <- "psbA.5"
ple_PS$NAME[ple_PS$NAME == "NODE_103641_length_722_cov_870.023112_g87213_i0.p1"] <- "psbA.6"
ple_PS$NAME[ple_PS$NAME == "NODE_41_length_30295_cov_377.508338_g29_i0.p20"] <- "cpcA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_2632_length_8241_cov_291.273384_g1815_i0.p7"] <- "cpcA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_14124_length_3420_cov_96.905886_g9397_i0.p4"] <- "cpeA"
ple_PS$NAME[ple_PS$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p5"] <- "rbcL"
#ple_PS$NAME[ple_PS$NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p6"] <- "ccmL.1"
#ple_PS$NAME[ple_PS$NAME == "NODE_6321_length_5503_cov_104.721179_g4191_i0.p8"] <- "ccmL.3"
#ple_PS$NAME[ple_PS$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p19"] <- "ccmL.2"
#ple_PS$NAME[ple_PS$NAME == "NODE_7497_length_5038_cov_463.560524_g4982_i0.p1"] <- "ccmM.1"
#ple_PS$NAME[ple_PS$NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p1"] <- "ccmM.3"
#ple_PS$NAME[ple_PS$NAME == "NODE_7638_length_4993_cov_64.689228_g5068_i0.p1"] <- "ccmM.4"
#ple_PS$NAME[ple_PS$NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p1"] <- "ccmM.5"
#ple_PS$NAME[ple_PS$NAME == "NODE_6744_length_5315_cov_351.019649_g464_i2.p1"] <- "ccmM.6"
#ple_PS$NAME[ple_PS$NAME == "NODE_6967_length_5229_cov_37.146431_g4191_i1.p1"] <- "ccmM.7"
#ple_PS$NAME[ple_PS$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p2"] <- "ccmM.8"


ple_PS <- tibble::column_to_rownames(ple_PS, "NAME")
ple_PS_t<-t(ple_PS)
ple_PS_t<-data.frame(ple_PS_t)
ple_PS_t$condition <- vec
ple_PS_t$condition1 <- vec1
ple_PS_t$condition2 <- vec2

ple_PS_tg<-ple_PS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_PS_tg.summary <- ple_PS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_PS_tg.summary_DRY_ONLY<-subset(ple_PS_tg.summary,condition1=="DRY")
ple_PS_plot<-ggplot(ple_PS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("psaA.1","psaA.2","psaA.3","psaA.4","psaA.5","psbA.1","psbA.2","psbA.3","psbA.4","psbA.5","psbA.6","cpcA.1","cpcA.2","cpeA","rbcL"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Photosynthesis (Pleurocapsales)")+ geom_hline(yintercept=20000, lty = 2)

#PLeurocapsales
#compatible solutes
ple_CS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_3722_length_7081_cov_19.808647_g2518_i0.p1"|
                                           NAME == "NODE_250_length_19831_cov_73.726237_g179_i0.p4"|
                                           NAME == "NODE_7834_length_4922_cov_286.902042_g5196_i0.p1"|
                                           NAME == "NODE_19716_length_2705_cov_64.979103_g13309_i0.p1"|
                                           NAME == "NODE_10088_length_4250_cov_27.823558_g6695_i0.p1"|
                                           NAME == "NODE_392_length_16787_cov_100.502573_g290_i0.p5"|
                                           NAME == "NODE_1835_length_9612_cov_306.439145_g1280_i0.p3"|
                                           NAME == "NODE_202_length_21322_cov_531.130971_g150_i0.p10"|
                                           NAME == "NODE_390_length_16796_cov_39.811278_g288_i0.p2"|
                                           NAME == "NODE_10320_length_4190_cov_25.909400_g6835_i0.p1"|
                                           NAME == "NODE_4743_length_6345_cov_35.767219_g3171_i0.p1"|
                                           NAME == "NODE_11434_length_3927_cov_15.738194_g5474_i1.p1"|
                                           NAME == "NODE_9217_length_4489_cov_12.658514_g6109_i0.p1"|
                                           NAME == "NODE_5658_length_5812_cov_106.701342_g3768_i0.p2"|
                                           NAME == "NODE_360_length_17275_cov_695.772585_g263_i0.p4"|
                                           NAME == "NODE_380_length_16917_cov_590.818867_g278_i0.p3"|
                                           NAME == "NODE_4957_length_6205_cov_42.477495_g3315_i0.p1"))

ple_CS$NAME[ple_CS$NAME == "NODE_3722_length_7081_cov_19.808647_g2518_i0.p1"] <- "spsA.1"
ple_CS$NAME[ple_CS$NAME == "NODE_250_length_19831_cov_73.726237_g179_i0.p4"] <- "spsA.2"
ple_CS$NAME[ple_CS$NAME == "NODE_7834_length_4922_cov_286.902042_g5196_i0.p1"] <- "spsA.3"
ple_CS$NAME[ple_CS$NAME == "NODE_19716_length_2705_cov_64.979103_g13309_i0.p1"] <- "spsA.4"
ple_CS$NAME[ple_CS$NAME == "NODE_10088_length_4250_cov_27.823558_g6695_i0.p1"] <- "spsA.5"
ple_CS$NAME[ple_CS$NAME == "NODE_392_length_16787_cov_100.502573_g290_i0.p5"] <- "spsA.6"
ple_CS$NAME[ple_CS$NAME == "NODE_1835_length_9612_cov_306.439145_g1280_i0.p3"] <- "spsA.7"
ple_CS$NAME[ple_CS$NAME == "NODE_202_length_21322_cov_531.130971_g150_i0.p10"] <- "spp"
ple_CS$NAME[ple_CS$NAME == "NODE_390_length_16796_cov_39.811278_g288_i0.p2"] <- "treS.1"
ple_CS$NAME[ple_CS$NAME == "NODE_10320_length_4190_cov_25.909400_g6835_i0.p1"] <- "treS.2"
ple_CS$NAME[ple_CS$NAME == "NODE_4743_length_6345_cov_35.767219_g3171_i0.p1"] <- "treS.3"
ple_CS$NAME[ple_CS$NAME == "NODE_11434_length_3927_cov_15.738194_g5474_i1.p1"] <- "treZ"
ple_CS$NAME[ple_CS$NAME == "NODE_9217_length_4489_cov_12.658514_g6109_i0.p1"] <- "ggpS.1"
ple_CS$NAME[ple_CS$NAME == "NODE_5658_length_5812_cov_106.701342_g3768_i0.p2"] <- "ggpS.2"
ple_CS$NAME[ple_CS$NAME == "NODE_360_length_17275_cov_695.772585_g263_i0.p4"] <- "ggpS.3"
ple_CS$NAME[ple_CS$NAME == "NODE_380_length_16917_cov_590.818867_g278_i0.p3"] <- "ggpP.1"
ple_CS$NAME[ple_CS$NAME == "NODE_4957_length_6205_cov_42.477495_g3315_i0.p1"] <- "ggpP.2"


ple_CS <- tibble::column_to_rownames(ple_CS, "NAME")
ple_CS_t<-t(ple_CS)
ple_CS_t<-data.frame(ple_CS_t)
ple_CS_t$condition <- vec
ple_CS_t$condition1 <- vec1
ple_CS_t$condition2 <- vec2

ple_CS_tg<-ple_CS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_CS_tg.summary <- ple_CS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_CS_tg.summary_DRY_ONLY<-subset(ple_CS_tg.summary,condition1=="DRY")
ple_CS_plot<-ggplot(ple_CS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("spsA.1","spsA.2","spsA.3","spsA.4","spsA.5","spsA.6","spsA.7","spp","treS.1","treS.2","treS.3","treZ","ggpP.1","ggpP.2","ggpS.1","ggpS.2","ggpS.3"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Compatible solutes (Pleurocapsales)")+ geom_hline(yintercept=100, lty = 2)

#PLeurocapsales
#EPS
ple_EPS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p9"|
                                           NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p4"|
                                            NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p6"|
                                            NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p6"|
                                            NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p2"|
                                            NAME == "NODE_5511_length_5888_cov_74.971109_g3672_i0.p2"|
                                            NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p2"|
                                            NAME == "NODE_113_length_24930_cov_55.630325_g86_i0.p5"|
                                            NAME == "NODE_5270_length_6026_cov_32.412901_g3515_i0.p2"|
                                            NAME == "NODE_7419_length_5065_cov_29.705729_g4926_i0.p2"|
                                          NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p4"|
                                            NAME =="NODE_842_length_13068_cov_93.375298_g10_i3.p1"|
                                            NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p1"|
                                            NAME == "NODE_265_length_19485_cov_195.053111_g188_i0.p2"|
                                            NAME == "NODE_1556_length_10299_cov_2219.429493_g1060_i1.p1"|
                                            NAME == "NODE_147_length_23217_cov_41.362124_g112_i0.p2"))

ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p9"] <- "wza"
ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p4"] <- "wzc"
ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p6"] <- "wzy"
ple_EPS$NAME[ple_EPS$NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p6"] <- "kpsD.1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p2"] <- "kpsD.2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_5511_length_5888_cov_74.971109_g3672_i0.p2"] <- "kpsD.3"
ple_EPS$NAME[ple_EPS$NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p2"] <- "kpsD.4"
ple_EPS$NAME[ple_EPS$NAME == "NODE_113_length_24930_cov_55.630325_g86_i0.p5"] <- "kpsD.5"
ple_EPS$NAME[ple_EPS$NAME == "NODE_5270_length_6026_cov_32.412901_g3515_i0.p2"] <- "kpsD.6"
ple_EPS$NAME[ple_EPS$NAME == "NODE_7419_length_5065_cov_29.705729_g4926_i0.p2"] <- "kpsD.7"
ple_EPS$NAME[ple_EPS$NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p4"] <- "kpsE.1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p1"] <- "kpsE.2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p1"] <- "kpsE.3"
ple_EPS$NAME[ple_EPS$NAME == "NODE_265_length_19485_cov_195.053111_g188_i0.p2"] <- "bcsA.1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_1556_length_10299_cov_2219.429493_g1060_i1.p1"] <- "bcsA.2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_147_length_23217_cov_41.362124_g112_i0.p2"] <- "bcsA.3"

ple_EPS <- tibble::column_to_rownames(ple_EPS, "NAME")
ple_EPS_t<-t(ple_EPS)
ple_EPS_t<-data.frame(ple_EPS_t)
ple_EPS_t$condition <- vec
ple_EPS_t$condition1 <- vec1
ple_EPS_t$condition2 <- vec2

ple_EPS_tg<-ple_EPS_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_EPS_tg.summary <- ple_EPS_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_EPS_tg.summary_DRY_ONLY<-subset(ple_EPS_tg.summary,condition1=="DRY")
ple_EPS_plot<-ggplot(ple_EPS_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("wza","wzc","wzy","kpsD.1","kpsD.2","kpsD.3","kpsD.4","kpsD.5","kpsD.6","kpsD.7","kpsE.1","kpsE.2","kpsE.3","bcsA.1","bcsA.2","bcsA.3"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("EPS export (Pleurocapsales)")+ geom_hline(yintercept=100, lty = 2)

#PLeurocapsales
#ST
ple_ST<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_356_length_17472_cov_223.329444_g259_i0.p3"|
                                         NAME == "NODE_561_length_14978_cov_299.386582_g399_i0.p6"|   
                                           NAME == "NODE_5592_length_5844_cov_58.440998_g3721_i0.p3"|
                                           NAME == "NODE_2371_length_8622_cov_428.524857_g1633_i0.p4"|
                                           NAME == "NODE_1112_length_11791_cov_391.154463_g777_i0.p3"|
                                           NAME == "NODE_240_length_20051_cov_536.532736_g170_i0.p5"|
                                           NAME == "NODE_203_length_21295_cov_555.963199_g151_i0.p8"|
                                           #NAME == "NODE_2779_length_8058_cov_11.868378_g1921_i0.p5"|
                                           NAME == "NODE_1024_length_12170_cov_58.469703_g717_i0.p5"|
                                           NAME == "NODE_637_length_14308_cov_149.386231_g451_i0.p6"|
                                           NAME == "NODE_538_length_15291_cov_285.033710_g380_i0.p7"|
                                           NAME == "NODE_8405_length_4735_cov_43.500644_g5582_i0.p3"|
                                           NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p3"|
                                           #NAME == "NODE_919_length_12712_cov_173.039164_g644_i0.p6"|
                                           NAME == "NODE_7155_length_5157_cov_103.254524_g1280_i3.p1"|
                                           NAME == "NODE_33276_length_1817_cov_191.809060_g23580_i0.p1"|
                                           NAME == "NODE_40265_length_1554_cov_175.994598_g29270_i0.p1"|
                                           NAME == "NODE_671_length_14102_cov_489.907406_g478_i0.p4"|
                                           NAME == "NODE_219_length_20772_cov_357.297164_g23_i2.p10"|
                                           #NAME == "NODE_18588_length_2821_cov_110.438501_g12506_i0.p1"|
                                           #NAME == "NODE_27764_length_2091_cov_220.791378_g19352_i0.p1"|
                                           NAME == "NODE_874_length_12891_cov_268.735372_g614_i0.p5"))

ple_ST$NAME[ple_ST$NAME == "NODE_356_length_17472_cov_223.329444_g259_i0.p3"] <- "alr4781"
ple_ST$NAME[ple_ST$NAME == "NODE_561_length_14978_cov_299.386582_g399_i0.p6"] <- "all0261.01"
ple_ST$NAME[ple_ST$NAME == "NODE_5592_length_5844_cov_58.440998_g3721_i0.p3"] <- "all0261.02"
ple_ST$NAME[ple_ST$NAME == "NODE_2371_length_8622_cov_428.524857_g1633_i0.p4"] <- "all0261.03"
ple_ST$NAME[ple_ST$NAME == "NODE_1112_length_11791_cov_391.154463_g777_i0.p3"] <- "all1823.01"
ple_ST$NAME[ple_ST$NAME == "NODE_240_length_20051_cov_536.532736_g170_i0.p5"] <- "all1823.02"
ple_ST$NAME[ple_ST$NAME == "NODE_203_length_21295_cov_555.963199_g151_i0.p8"] <- "all1823.03"
#ple_ST$NAME[ple_ST$NAME == "NODE_2779_length_8058_cov_11.868378_g1921_i0.p5"] <- "all1823.4"
ple_ST$NAME[ple_ST$NAME == "NODE_1024_length_12170_cov_58.469703_g717_i0.p5"] <- "all1823.05"
ple_ST$NAME[ple_ST$NAME == "NODE_637_length_14308_cov_149.386231_g451_i0.p6"] <- "all1823.06"
ple_ST$NAME[ple_ST$NAME == "NODE_538_length_15291_cov_285.033710_g380_i0.p7"] <- "all1823.07"
ple_ST$NAME[ple_ST$NAME == "NODE_8405_length_4735_cov_43.500644_g5582_i0.p3"] <- "all1823.08"
ple_ST$NAME[ple_ST$NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p3"] <- "all1823.09"
#ple_ST$NAME[ple_ST$NAME == "NODE_919_length_12712_cov_173.039164_g644_i0.p6"] <- "all1823.10"
ple_ST$NAME[ple_ST$NAME == "NODE_7155_length_5157_cov_103.254524_g1280_i3.p1"] <- "all1823.11"
ple_ST$NAME[ple_ST$NAME == "NODE_33276_length_1817_cov_191.809060_g23580_i0.p1"] <- "all1823.12"
ple_ST$NAME[ple_ST$NAME == "NODE_40265_length_1554_cov_175.994598_g29270_i0.p1"] <- "all1823.13"
ple_ST$NAME[ple_ST$NAME == "NODE_671_length_14102_cov_489.907406_g478_i0.p4"] <- "all1823.14"
ple_ST$NAME[ple_ST$NAME == "NODE_219_length_20772_cov_357.297164_g23_i2.p10"] <- "all1823.15"
#ple_ST$NAME[ple_ST$NAME == "NODE_18588_length_2821_cov_110.438501_g12506_i0.p1"] <- "all1823.16"
#ple_ST$NAME[ple_ST$NAME == "NODE_27764_length_2091_cov_220.791378_g19352_i0.p1"] <- "all1823.17"
ple_ST$NAME[ple_ST$NAME == "NODE_874_length_12891_cov_268.735372_g614_i0.p5"] <- "alr2532"

ple_ST <- tibble::column_to_rownames(ple_ST, "NAME")
ple_ST_t<-t(ple_ST)
ple_ST_t<-data.frame(ple_ST_t)
ple_ST_t$condition <- vec
ple_ST_t$condition1 <- vec1
ple_ST_t$condition2 <- vec2

ple_ST_tg<-ple_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_ST_tg.summary <- ple_ST_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_ST_tg.summary_DRY_ONLY<-subset(ple_ST_tg.summary,condition1=="DRY")
ple_ST_plot<-ggplot(ple_ST_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("psaA.1"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Sugar transport (Pleurocapsales)")+ geom_hline(yintercept=100, lty = 2)




#PLeurocapsales
#Nitrogen
ple_NT<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_16766_length_3034_cov_278.412023_g11235_i0.p1"|
                                         NAME == "NODE_35_length_31865_cov_1182.712947_g24_i0.p5"|
                                           NAME == "NODE_36114_length_1700_cov_60.636755_g25872_i0.p1"|
                                           NAME == "NODE_18172_length_2867_cov_38.164281_g12221_i0.p1"|
                                           NAME == "NODE_690_length_13975_cov_291.506690_g485_i0.p1"|
                                           NAME == "NODE_6790_length_5299_cov_25.533678_g4530_i0.p2"|
                                           NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p2"|
                                           NAME == "NODE_6571_length_5392_cov_27.662907_g4379_i0.p1"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p1"|
                                           NAME == "NODE_4951_length_6207_cov_67.991686_g3311_i0.p1"|
                                           NAME == "NODE_193475_length_465_cov_2.278061_g176045_i0.p1"|
                                           NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p4"|
                                           NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p2"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p2"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p3"|
                                           NAME == "NODE_1508_length_10443_cov_71.441273_g1057_i0.p6"|
                                           NAME == "NODE_3011_length_7804_cov_26.634459_g2060_i0.p5"|
                                           NAME == "NODE_2611_length_8278_cov_64.763681_g1801_i0.p5"))

ple_NT$NAME[ple_NT$NAME == "NODE_16766_length_3034_cov_278.412023_g11235_i0.p1"] <- "alr0990.1"
ple_NT$NAME[ple_NT$NAME == "NODE_35_length_31865_cov_1182.712947_g24_i0.p5"] <- "alr0990.2"
ple_NT$NAME[ple_NT$NAME == "NODE_36114_length_1700_cov_60.636755_g25872_i0.p1"] <- "alr0990.3"
ple_NT$NAME[ple_NT$NAME == "NODE_18172_length_2867_cov_38.164281_g12221_i0.p1"] <- "alr0990.4"
ple_NT$NAME[ple_NT$NAME == "NODE_690_length_13975_cov_291.506690_g485_i0.p1"] <- "alr0991"
ple_NT$NAME[ple_NT$NAME == "NODE_6790_length_5299_cov_25.533678_g4530_i0.p2"] <- "nrtP.1"
ple_NT$NAME[ple_NT$NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p2"] <- "nrtP.2"
ple_NT$NAME[ple_NT$NAME == "NODE_6571_length_5392_cov_27.662907_g4379_i0.p1"] <- "urtA.1"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p1"] <- "urtA.2"
ple_NT$NAME[ple_NT$NAME == "NODE_4951_length_6207_cov_67.991686_g3311_i0.p1"] <- "urtA.3"
ple_NT$NAME[ple_NT$NAME == "NODE_193475_length_465_cov_2.278061_g176045_i0.p1"] <- "urtA.4"
ple_NT$NAME[ple_NT$NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p4"] <- "urtB"
ple_NT$NAME[ple_NT$NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p2"] <- "urtC.1"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p2"] <- "urtC.2"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p3"] <- "urtD.1"
ple_NT$NAME[ple_NT$NAME == "NODE_1508_length_10443_cov_71.441273_g1057_i0.p6"] <- "urtD.2"
ple_NT$NAME[ple_NT$NAME == "NODE_3011_length_7804_cov_26.634459_g2060_i0.p5"] <- "urtE.1"
ple_NT$NAME[ple_NT$NAME == "NODE_2611_length_8278_cov_64.763681_g1801_i0.p5"] <- "urtE.2"

ple_NT <- tibble::column_to_rownames(ple_NT, "NAME")
ple_NT_t<-t(ple_NT)
ple_NT_t<-data.frame(ple_NT_t)
ple_NT_t$condition <- vec
ple_NT_t$condition1 <- vec1
ple_NT_t$condition2 <- vec2

ple_NT_tg<-ple_NT_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)

ple_NT_tg.summary <- ple_NT_tg %>%
  group_by(condition2,condition1,condition,Gene) %>%
  summarise(
    sd = sd(TPM, na.rm = TRUE),
    len = mean(TPM))

ple_NT_tg.summary_DRY_ONLY<-subset(ple_NT_tg.summary,condition1=="DRY")
ple_NT_plot<-ggplot(ple_NT_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  # ylim(0,65000) +
  scale_fill_manual(values = c("#37ABC8", "#005544")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("psaA.1"))+
  xlab("") +
  ylab("Transcripts per million (TPM)")+
  ggtitle("Nitrogen (Pleurocapsales)")+ geom_hline(yintercept=100, lty = 2)

plot_grid(ple_HK_plot,ple_PS_plot,ple_NT_plot,ple_EPS_plot,ple_ST_plot,ple_CS_plot,ncol=2, align="hv")

library(tidyverse)
library(rstatix)
library(ggpubr)
nos_CS_tg_DRY<-subset(nos_CS_tg, condition1=="DRY")
nos_CS_tg_DRY$logTPM<-log(nos_CS_tg_DRY$TPM+1)
stat.test <- nos_CS_tg_DRY %>%
  group_by(Gene) %>%
  t_test(logTPM ~ condition2, paired=FALSE) %>%
  #  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
# Create the plot
myplot  <-ggboxplot(
  ple_PS_tg_DRY, x = "condition2", y = "TPM",
  fill = "condition1", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~Gene, ncol=11)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "condition2")
myplot + stat_pvalue_manual(stat.test, label = "p.signif")

plot_grid(nos_PS_plot,nos_EPS_plot,nos_ST_plot,nos_CS_plot,ncol=2, align = "hv")

#heatmap
library(pheatmap)
nos_PS_hm<-nos_PS
nos_PS_hm$type<-"PS"
nos_CS_hm<-nos_CS
nos_CS_hm$type<-"CS"
nos_ST_hm<-nos_ST
nos_ST_hm$type<-"ST"
nos_EPS_hm<-nos_EPS
nos_EPS_hm$type<-"EPS"
nos_hm <- rbind(nos_PS_hm, nos_CS_hm,nos_ST_hm,nos_EPS_hm)
nos_hm_mat <- as.matrix(nos_CS_hm[,1:6])
nos_hm_mat_ln<-log(nos_hm_mat)
nos_hm_mat_scale <- scale(nos_hm_mat)
pheatmap(nos_hm_mat_ln,main = "pheatmap default")
