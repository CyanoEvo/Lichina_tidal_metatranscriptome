library(tidyverse)
library(janitor)

vec<-c("HT","HT","HT","LT","LT","LT","HT","HT","HT","LT","LT","LT")
vec1<-c("DRY","DRY","DRY","DRY","DRY","DRY","WET","WET","WET","WET","WET","WET")
vec2<-c("HT_DRY","HT_DRY","HT_DRY","LT_DRY","LT_DRY","LT_DRY","HT_WET","HT_WET","HT_WET","LT_WET","LT_WET","LT_WET")

nos_diff_ex_READS_df <- tibble::rownames_to_column(nos_diff_ex_TPM, "NAME")
ple_diff_ex_READS_df <- tibble::rownames_to_column(ple_diff_ex_TPM, "NAME")
asc_diff_ex_READS_df <- tibble::rownames_to_column(asc_diff_ex_TPM, "NAME")


#NOSTOCALES HK GENES
#selected from Pinto et al 2012
nos_HK<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"|
                                         NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"|
                                           NAME == "NODE_889_length_12823_cov_960.285569_g623_i0.p1"|
                                           NAME == "NODE_106_length_25244_cov_2665.611577_g79_i0.p7"|
                                           NAME == "NODE_5235_length_6040_cov_573.266130_g3493_i0.p2"|
                                           NAME == "NODE_1089_length_11885_cov_831.030393_g759_i0.p1"))

nos_HK$NAME[nos_HK$NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"] <- "petB"
nos_HK$NAME[nos_HK$NAME == "NODE_889_length_12823_cov_960.285569_g623_i0.p1"] <- "ppc"
nos_HK$NAME[nos_HK$NAME == "NODE_5235_length_6040_cov_573.266130_g3493_i0.p2"] <- "rps1B"
nos_HK$NAME[nos_HK$NAME == "NODE_106_length_25244_cov_2665.611577_g79_i0.p7"] <- "rpoA"
nos_HK$NAME[nos_HK$NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"] <- "ilvD"
nos_HK$NAME[nos_HK$NAME == "NODE_1089_length_11885_cov_831.030393_g759_i0.p1"] <- "secA"

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
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("petB","ppc","rps1B","rpoA","ilvD","secA"))+
  xlab("") +
  ylab("Nostocales transcripts per million (TPM)")

ple_HK<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_461_length_15967_cov_59.193784_g335_i0.p7"|
                                           NAME == "NODE_864_length_12967_cov_50.308128_g605_i0.p1"|
                                           NAME == "NODE_1208_length_11446_cov_85.420206_g850_i0.p1"|
                                           NAME == "NODE_105_length_25260_cov_67.393100_g78_i0.p1"|
                                           NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p1"|
                                           NAME == "NODE_97_length_25645_cov_423.605702_g72_i0.p6"|
                                           NAME == "NODE_728_length_13753_cov_53.621345_g422_i1.p4"|
                                           NAME == "NODE_46375_length_1389_cov_16.399696_g34339_i0.p1"|
                                           NAME == "NODE_2090_length_9130_cov_20.736889_g1452_i0.p2"|
                                           NAME == "NODE_12921_length_3621_cov_55.577508_g8558_i0.p1"))

ple_HK$NAME[ple_HK$NAME == "NODE_461_length_15967_cov_59.193784_g335_i0.p7"] <- "petB"
ple_HK$NAME[ple_HK$NAME == "NODE_864_length_12967_cov_50.308128_g605_i0.p1"] <- "ppc.1"
ple_HK$NAME[ple_HK$NAME == "NODE_1208_length_11446_cov_85.420206_g850_i0.p1"] <- "ppc.2"
ple_HK$NAME[ple_HK$NAME == "NODE_105_length_25260_cov_67.393100_g78_i0.p1"] <- "ppc.3"
ple_HK$NAME[ple_HK$NAME == "NODE_207_length_21109_cov_427.554478_g155_i0.p1"] <- "ppc.4"
ple_HK$NAME[ple_HK$NAME == "NODE_97_length_25645_cov_423.605702_g72_i0.p6"] <- "rps1B.1"
ple_HK$NAME[ple_HK$NAME == "NODE_728_length_13753_cov_53.621345_g422_i1.p4"] <- "rps1B.2"
ple_HK$NAME[ple_HK$NAME == "NODE_46375_length_1389_cov_16.399696_g34339_i0.p1"] <- "rps1B.3"
ple_HK$NAME[ple_HK$NAME == "NODE_2090_length_9130_cov_20.736889_g1452_i0.p2"] <- "rps1B.4"
ple_HK$NAME[ple_HK$NAME == "NODE_12921_length_3621_cov_55.577508_g8558_i0.p1"] <- "rpoA.1"



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
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("petB","ppc.1","ppc.2","ppc.3","ppc.4","rps1B.1","rps1B.2","rps1B.3","rps1B.4","rpoA.1"))+
  xlab("") +
  ylab("Pleurocapsales transcripts per million (TPM)")


#NOSTOCALES PS GENES
nos_PS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p1"|
#                                         NAME == "NODE_33200_length_1820_cov_185.679450_g23515_i0.p1"|
                                        NAME == "NODE_81960_length_870_cov_125019.715182_g66510_i0.p1"|
                                          NAME == "NODE_7966_length_4882_cov_15047.541485_g5282_i0.p5"|
                                          NAME == "NODE_3079_length_7723_cov_9104.390327_g2104_i0.p2"|
                                          NAME == "NODE_216_length_20811_cov_1216.247613_g161_i0.p3"))

nos_PS$NAME[nos_PS$NAME == "NODE_2427_length_8550_cov_29497.134128_g1675_i0.p1"] <- "PsaA"
#nos_PS$NAME[nos_PS$NAME == "NODE_33200_length_1820_cov_185.679450_g23515_i0.p1"] <- "PsbA.1"
nos_PS$NAME[nos_PS$NAME == "NODE_81960_length_870_cov_125019.715182_g66510_i0.p1"] <- "PsbA.2"
nos_PS$NAME[nos_PS$NAME == "NODE_7966_length_4882_cov_15047.541485_g5282_i0.p5"] <- "CpcA"
nos_PS$NAME[nos_PS$NAME == "NODE_3079_length_7723_cov_9104.390327_g2104_i0.p2"] <- "RbcL"
nos_PS$NAME[nos_PS$NAME == "NODE_216_length_20811_cov_1216.247613_g161_i0.p3"] <- "CccM"

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
 # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
#           alpha = .1,fill = "black")+
#  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 0, ymax = Inf,
#           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("PsaA","PsbA1","PsbA2","CpcA","CpeA","RbcL","CccM"))+
  scale_x_discrete(limits = c("PsaA","PsbA.2","CpcA","","RbcL","CccM"))+
   xlab("") +
  ylab("Nostocales transcripts per million (TPM)")

#Pleurocapsales PS GENES
ple_PS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_19751_length_2701_cov_4954.788813_g13333_i0.p1"|
                                            NAME == "NODE_1511_length_10440_cov_4349.966818_g1058_i0.p1"|
                                            NAME == "NODE_1375_length_10838_cov_1438.095866_g969_i0.p1"|
                                              NAME == "NODE_4757_length_6333_cov_56.899042_g3181_i0.p1"|
                                              NAME == "NODE_70_length_27121_cov_248.993382_g53_i0.p6"|
                                           NAME == "NODE_23791_length_2350_cov_26283.451910_g802_i3.p1"|
                                              NAME == "NODE_174_length_22195_cov_217.718380_g132_i0.p8"|
                                           NAME == "NODE_41_length_30295_cov_377.508338_g29_i0.p20"|
                                              NAME == "NODE_2632_length_8241_cov_291.273384_g1815_i0.p7"|
                                           NAME == "NODE_14124_length_3420_cov_96.905886_g9397_i0.p4"|
                                              NAME == "NODE_47921_length_1351_cov_28.586072_g35681_i0.p2"|
                                           NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p5"|
                                              NAME == "NODE_7497_length_5038_cov_463.560524_g4982_i0.p1"|
                                              NAME == "NODE_9293_length_4464_cov_51.033478_g6166_i0.p1"|
                                           NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p1"|
                                            NAME == "NODE_7638_length_4993_cov_64.689228_g5068_i0.p1"|
                                            NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p1"|
                                            NAME == "NODE_6744_length_5315_cov_351.019649_g464_i2.p1"|
                                            NAME == "NODE_6967_length_5229_cov_37.146431_g4191_i1.p1"|
                                            NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p2"
))

ple_PS$NAME[ple_PS$NAME == "NODE_19751_length_2701_cov_4954.788813_g13333_i0.p1"] <- "PsaA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_1511_length_10440_cov_4349.966818_g1058_i0.p1"] <- "PsaA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_1375_length_10838_cov_1438.095866_g969_i0.p1"] <- "PsaA.3"
ple_PS$NAME[ple_PS$NAME == "NODE_4757_length_6333_cov_56.899042_g3181_i0.p1"] <- "PsaA.4"
ple_PS$NAME[ple_PS$NAME == "NODE_70_length_27121_cov_248.993382_g53_i0.p6"] <- "PsbA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_23791_length_2350_cov_26283.451910_g802_i3.p1"] <- "PsbA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_174_length_22195_cov_217.718380_g132_i0.p8"] <- "PsbA.3"
ple_PS$NAME[ple_PS$NAME == "NODE_41_length_30295_cov_377.508338_g29_i0.p20"] <- "CpcA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_2632_length_8241_cov_291.273384_g1815_i0.p7"] <- "CpcA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_14124_length_3420_cov_96.905886_g9397_i0.p4"] <- "CpeA.1"
ple_PS$NAME[ple_PS$NAME == "NODE_47921_length_1351_cov_28.586072_g35681_i0.p2"] <- "CpeA.2"
ple_PS$NAME[ple_PS$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p5"] <- "RbcL"
ple_PS$NAME[ple_PS$NAME == "NODE_7497_length_5038_cov_463.560524_g4982_i0.p1"] <- "CcmM.1"
ple_PS$NAME[ple_PS$NAME == "NODE_9293_length_4464_cov_51.033478_g6166_i0.p1"] <- "CcmM.2"
ple_PS$NAME[ple_PS$NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p1"] <- "CcmM.3"
ple_PS$NAME[ple_PS$NAME == "NODE_7638_length_4993_cov_64.689228_g5068_i0.p1"] <- "CcmM.4"
ple_PS$NAME[ple_PS$NAME == "NODE_6273_length_5527_cov_89.561239_g4156_i0.p1"] <- "CcmM.5"
ple_PS$NAME[ple_PS$NAME == "NODE_6744_length_5315_cov_351.019649_g464_i2.p1"] <- "CcmM.6"
ple_PS$NAME[ple_PS$NAME == "NODE_6967_length_5229_cov_37.146431_g4191_i1.p1"] <- "CcmM.7"
ple_PS$NAME[ple_PS$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p2"] <- "CcmM.8"

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


ple_PS_plot<-ggplot(ple_PS_tg.summary, aes(x = Gene, y = len)) +
 # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
#           alpha = .1,fill = "black")+
 # annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 0, ymax = Inf,
  #         alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  ylim(0,65000) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("PsaA.1","PsaA.2","PsaA.3","PsaA.4","PsbA.1","PsbA.2","PsbA.3","CpcA.1","CpcA.2","CpeA.1","CpeA.2","RbcL","CcmM.1","CcmM.2","CcmM.3","CcmM.4","CcmM.5","CcmM.6","CcmM.7","CcmM.8"))+
  xlab("") +
  ylab("Pleurocapsales transcripts per million (TPM)")

plot_grid(nos_PS_plot,ple_PS_plot,ncol=2,align = "hv")

#NOSTOCALES EPS GENES
nos_EPS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p1"|
                                            NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p2"|
                                            NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p3"|
                                            NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p5"|
                                            NAME == "NODE_14243_length_3402_cov_434.991589_g9475_i0.p1"|
                                            NAME == "NODE_13298_length_3560_cov_641.725839_g8837_i0.p1"|
                                           NAME == "NODE_13762_length_3480_cov_54.459935_g9151_i0.p2"|
                                          NAME == "NODE_24263_length_2317_cov_24.334670_g16656_i0.p2" |
                                          NAME == "NODE_12244_length_3756_cov_50.163182_g8122_i0.p1"|
                                           NAME == "NODE_65157_length_1047_cov_216.448665_g50892_i0.p1"|
                                           NAME == "NODE_138615_length_583_cov_151.256863_g121521_i0.p1"|
                                           NAME == "NODE_4220_length_6703_cov_57.689894_g2834_i0.p6" |
                                      NAME == "NODE_1851_length_9582_cov_58.061521_g1292_i0.p1"))

nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p1"] <- "Wzc"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p2"] <- "Wzy"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p3"] <- "Wzx"
nos_EPS$NAME[nos_EPS$NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p5"] <- "Wzx1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_14243_length_3402_cov_434.991589_g9475_i0.p1"] <- "KpsD"
nos_EPS$NAME[nos_EPS$NAME == "NODE_13298_length_3560_cov_641.725839_g8837_i0.p1"] <- "KpsE"
nos_EPS$NAME[nos_EPS$NAME == "NODE_13762_length_3480_cov_54.459935_g9151_i0.p2"] <- "KpsT1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_24263_length_2317_cov_24.334670_g16656_i0.p2"] <- "KpsT2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_12244_length_3756_cov_50.163182_g8122_i0.p1"] <- "KpsT3"
nos_EPS$NAME[nos_EPS$NAME == "NODE_65157_length_1047_cov_216.448665_g50892_i0.p1"] <- "KpsT4"
nos_EPS$NAME[nos_EPS$NAME == "NODE_138615_length_583_cov_151.256863_g121521_i0.p1"] <- "KpsM1"
nos_EPS$NAME[nos_EPS$NAME == "NODE_4220_length_6703_cov_57.689894_g2834_i0.p6"] <- "KpsM2"
nos_EPS$NAME[nos_EPS$NAME == "NODE_1851_length_9582_cov_58.061521_g1292_i0.p1"] <- "BcsA"

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
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  annotate("rect", xmin = 5.5, xmax = 9.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  annotate("rect", xmin = 10.5, xmax = 12.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  ylim(0,200) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),      
        aspect.ratio = 0.5) +
        scale_x_discrete(limits = c("Wza","Wzc","Wzy","Wzx","","KpsD","KpsE","KpsT4","KpsM1","","BcsA","BcsC"))+
        xlab("") +
        ylab("Nostocales transcripts per million (TPM)")


#PLEUROCAPSALES EPS GENES

ple_EPS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p9"|
                                            NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p4"|
                                            NAME == "NODE_355_length_17525_cov_279.467683_g258_i0.p1"|
                                            NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p6"|
                                            NAME == "NODE_7914_length_4897_cov_80.121683_g5247_i0.p1" |
                                            NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p6"|
                                            NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p2"|
                                            NAME == "NODE_5511_length_5888_cov_74.971109_g3672_i0.p2"|
                                            NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p2"|
                                            NAME == "NODE_113_length_24930_cov_55.630325_g86_i0.p5"|
                                            NAME == "NODE_5270_length_6026_cov_32.412901_g3515_i0.p2"|
                                            NAME == "NODE_7419_length_5065_cov_29.705729_g4926_i0.p2"|
                                            NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p4"|
                                            NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p1"|
                                            NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p1"|
                                            NAME == "NODE_15936_length_3149_cov_53.289012_g10633_i0.p1"|
                                            NAME == "NODE_425_length_16344_cov_102.803085_g315_i0.p3"|
                                            NAME == "NODE_558_length_14984_cov_218.878613_g119_i1.p1"|
                                            NAME == "NODE_107_length_25205_cov_75.313982_g80_i0.p1"|
                                            NAME == "NODE_115_length_24677_cov_151.241058_g88_i0.p12"|
                                            NAME == "NODE_17_length_36777_cov_40.073507_g12_i0.p22"))

ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p9"] <- "Wza"
ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p4"] <- "Wzc1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_355_length_17525_cov_279.467683_g258_i0.p1"] <- "Wzc2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p6"] <- "Wzy"
ple_EPS$NAME[ple_EPS$NAME == "NODE_7914_length_4897_cov_80.121683_g5247_i0.p1"] <- "Wzx"
ple_EPS$NAME[ple_EPS$NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p6"] <- "KpsD1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p2"] <- "KpsD2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_5511_length_5888_cov_74.971109_g3672_i0.p2"] <- "KpsD3"
ple_EPS$NAME[ple_EPS$NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p2"] <- "KpsD4"
ple_EPS$NAME[ple_EPS$NAME == "NODE_113_length_24930_cov_55.630325_g86_i0.p5"] <- "KpsD5"
ple_EPS$NAME[ple_EPS$NAME == "NODE_5270_length_6026_cov_32.412901_g3515_i0.p2"] <- "KpsD6"
ple_EPS$NAME[ple_EPS$NAME == "NODE_7419_length_5065_cov_29.705729_g4926_i0.p2"] <- "KpsD7"
ple_EPS$NAME[ple_EPS$NAME == "NODE_8_length_39956_cov_673.668856_g6_i0.p4"] <- "KpsE1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_842_length_13068_cov_93.375298_g10_i3.p1"] <- "KpsE2"
ple_EPS$NAME[ple_EPS$NAME == "NODE_6409_length_5463_cov_259.661967_g4257_i0.p1"] <- "KpsE3"
ple_EPS$NAME[ple_EPS$NAME == "NODE_15936_length_3149_cov_53.289012_g10633_i0.p1"] <- "KpsE4"
ple_EPS$NAME[ple_EPS$NAME == "NODE_425_length_16344_cov_102.803085_g315_i0.p3"] <- "KpsE5"
ple_EPS$NAME[ple_EPS$NAME == "NODE_558_length_14984_cov_218.878613_g119_i1.p1"] <- "KpsE6"
ple_EPS$NAME[ple_EPS$NAME == "NODE_107_length_25205_cov_75.313982_g80_i0.p1"] <- "KpsE7"
ple_EPS$NAME[ple_EPS$NAME == "NODE_115_length_24677_cov_151.241058_g88_i0.p12"] <- "KpsM1"
ple_EPS$NAME[ple_EPS$NAME == "NODE_17_length_36777_cov_40.073507_g12_i0.p22"] <- "KpsM2"

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

ple_EPS_plot<-ggplot(ple_EPS_tg.summary, aes(x = Gene, y = len))+
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  annotate("rect", xmin = 5.5, xmax = 9.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  annotate("rect", xmin = 10.5, xmax = 12.5, ymin = 0, ymax = Inf,
           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  ylim(0,500) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),
        aspect.ratio = 0.5) +
  scale_x_discrete(limits = c("Wza","Wzc2","Wzy","Wzx","","KpsD1","KpsE1","KpsT","KpsM2","","BcsA","BcsC"))+
  xlab("") +
  ylab("Pleurocapsales transcripts per million (TPM)")



plot_grid(nos_EPS_plot,ple_EPS_plot,ncol=2,align = "hv")

#NOSTOCALES SUGAR GENES
nos_ST<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_1164_length_11617_cov_1128.037509_g816_i0.p3"|
                                           NAME == "NODE_15345_length_3231_cov_364.330906_g10223_i0.p1"|
                                           NAME == "NODE_54162_length_1222_cov_147.726719_g41060_i0.p1"|
                                           NAME == "NODE_2840_length_7989_cov_439.225366_g1958_i0.p2"|
                                           NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p7"|
                                           NAME == "NODE_2181_length_8992_cov_169.673954_g1508_i0.p3"|
                                           NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p3"|
                                           NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p8"|
                                           NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p6"|
                                           NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p12"|
                                           NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p7"|
                                           NAME == "NODE_25957_length_2203_cov_278.194836_g17955_i0.p1"|
                                           NAME == "NODE_8066_length_4847_cov_55.883117_g5352_i0.p1"|
                                           NAME == "NODE_20031_length_2674_cov_42.771242_g13532_i0.p1"|
                                           NAME == "NODE_3315_length_7470_cov_2218.795728_g2261_i0.p2"|
                                           NAME == "NODE_2792_length_8042_cov_45.165014_g1929_i0.p2"|
                                           NAME == "NODE_36_length_31522_cov_335.764221_g25_i0.p2"|
                                           NAME == "NODE_1251_length_11290_cov_1903.872515_g881_i0.p1"|
                                           NAME == "NODE_216_length_20811_cov_1216.247613_g161_i0.p2"|
                                           NAME == "NODE_18896_length_2787_cov_1020.099853_g12727_i0.p1"|
                                           NAME == "NODE_87_length_26090_cov_1811.419303_g66_i0.p6"|
                                           NAME == "NODE_5932_length_5676_cov_437.333214_g3936_i0.p1"|
                                           NAME == "NODE_11718_length_3863_cov_73.964908_g7768_i0.p1"|
                                           NAME == "NODE_5640_length_5823_cov_43.254087_g3757_i0.p2"|
                                           NAME == "NODE_10804_length_4072_cov_1164.011003_g7147_i0.p1"|
                                           NAME == "NODE_8127_length_4825_cov_879.277778_g5394_i0.p1"|
                                           NAME == "NODE_1658_length_10060_cov_419.707019_g1154_i0.p2"|
                                           NAME == "NODE_1912_length_9446_cov_1628.686546_g1338_i0.p1"|
                                           NAME == "NODE_185_length_21819_cov_1942.958107_g139_i0.p2"|
                                           NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"|
                                           NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"))

nos_ST$NAME[nos_ST$NAME == "NODE_1164_length_11617_cov_1128.037509_g816_i0.p3"] <- "GlsCD1"
nos_ST$NAME[nos_ST$NAME == "NODE_15345_length_3231_cov_364.330906_g10223_i0.p1"] <- "GlsCD2"
nos_ST$NAME[nos_ST$NAME == "NODE_54162_length_1222_cov_147.726719_g41060_i0.p1"] <- "GlsCD3"
nos_ST$NAME[nos_ST$NAME == "NODE_2840_length_7989_cov_439.225366_g1958_i0.p2"] <- "GlsCD4"
nos_ST$NAME[nos_ST$NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p7"] <- "GlsP"
nos_ST$NAME[nos_ST$NAME == "NODE_2181_length_8992_cov_169.673954_g1508_i0.p3"] <- "GlsR1"
nos_ST$NAME[nos_ST$NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p3"] <- "GlsR2"
nos_ST$NAME[nos_ST$NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p8"] <- "GlsQ1"
nos_ST$NAME[nos_ST$NAME == "NODE_1034_length_12131_cov_160.349643_g722_i0.p6"] <- "GlsQ2"
nos_ST$NAME[nos_ST$NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p12"] <- "FraC"
nos_ST$NAME[nos_ST$NAME == "NODE_258_length_19637_cov_209.469485_g184_i0.p7"] <- "FraD"
nos_ST$NAME[nos_ST$NAME == "NODE_25957_length_2203_cov_278.194836_g17955_i0.p1"] <- "SbtA1"
nos_ST$NAME[nos_ST$NAME == "NODE_8066_length_4847_cov_55.883117_g5352_i0.p1"] <- "SbtA2"
nos_ST$NAME[nos_ST$NAME == "NODE_20031_length_2674_cov_42.771242_g13532_i0.p1"] <- "SbtA3"
nos_ST$NAME[nos_ST$NAME == "NODE_3315_length_7470_cov_2218.795728_g2261_i0.p2"] <- "BicA1"
nos_ST$NAME[nos_ST$NAME == "NODE_2792_length_8042_cov_45.165014_g1929_i0.p2"] <- "BicA2"
nos_ST$NAME[nos_ST$NAME == "NODE_36_length_31522_cov_335.764221_g25_i0.p2"] <- "BicA3"
nos_ST$NAME[nos_ST$NAME == "NODE_1251_length_11290_cov_1903.872515_g881_i0.p1"] <- "NdhF1"
nos_ST$NAME[nos_ST$NAME == "NODE_216_length_20811_cov_1216.247613_g161_i0.p2"] <- "NdhF2"
nos_ST$NAME[nos_ST$NAME == "NODE_18896_length_2787_cov_1020.099853_g12727_i0.p1"] <- "NdhF3"
nos_ST$NAME[nos_ST$NAME == "NODE_87_length_26090_cov_1811.419303_g66_i0.p6"] <- "Alr3705"
nos_ST$NAME[nos_ST$NAME == "NODE_5932_length_5676_cov_437.333214_g3936_i0.p1"] <- "InvAB"
nos_ST$NAME[nos_ST$NAME == "NODE_11718_length_3863_cov_73.964908_g7768_i0.p1"] <- "SBP"
nos_ST$NAME[nos_ST$NAME == "NODE_5640_length_5823_cov_43.254087_g3757_i0.p2"] <- "TMDL"
nos_ST$NAME[nos_ST$NAME == "NODE_10804_length_4072_cov_1164.011003_g7147_i0.p1"] <- "TreS"
nos_ST$NAME[nos_ST$NAME == "NODE_8127_length_4825_cov_879.277778_g5394_i0.p1"] <- "TreY"
nos_ST$NAME[nos_ST$NAME == "NODE_1658_length_10060_cov_419.707019_g1154_i0.p2"] <- "Ggps"
nos_ST$NAME[nos_ST$NAME == "NODE_1912_length_9446_cov_1628.686546_g1338_i0.p1"] <- "HK_pyrG"
nos_ST$NAME[nos_ST$NAME == "NODE_185_length_21819_cov_1942.958107_g139_i0.p2"] <- "HK_rpoB"
nos_ST$NAME[nos_ST$NAME == "NODE_89_length_25912_cov_1373.786215_g67_i0.p12"] <- "HK_pet1"
nos_ST$NAME[nos_ST$NAME == "NODE_33747_length_1795_cov_127.213124_g23960_i0.p1"] <- "HK_ilvD"


nos_ST <- tibble::column_to_rownames(nos_ST, "NAME")
nos_ST_t<-t(nos_ST)
nos_ST_t<-data.frame(nos_ST_t)
nos_ST_t$condition <- vec
nos_ST_t$condition1 <- vec1
nos_ST_t$condition2 <- vec2

nos_ST_tg<-nos_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
nos_ST_tg.summary <- nos_ST_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))
nos_ST_tg.summary_DRY_ONLY<-subset(nos_ST_tg.summary,condition1=="DRY")
nos_ST_plot<-ggplot(nos_ST_tg.summary_DRY_ONLY, aes(x = Gene, y = len)) +
 # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
#           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
 # ylim(0,175) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("GlsCD1","GlsCD2","GlsCD3","GlsCD4","GlsP","GlsR1","GlsR2","GlsQ1","GlsQ2","","FraC","FraD","","SbtA1","SbtA2","SbtA3","","BicA1","BicA2","BicA3","","NdhF1","NdhF2","NdhF3","","Alr3705","","InvAB","","SBP","","TMDL"))+
  scale_x_discrete(limits = c("GlsCD1","GlsCD2","GlsCD3","GlsCD4","GlsP","GlsR1","GlsR2","GlsQ1","GlsQ2","","TreS","TreY","Ggps","HK_pyrG","HK_rpoB","HK_pet1","HK_ilvD"))+
    xlab("") +
  ylab("Nostocales transcripts per million (TPM)")

#PLEUROCAPSALES SUGAR GENES
ple_ST<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_356_length_17472_cov_223.329444_g259_i0.p3"|
                                           NAME == "NODE_1112_length_11791_cov_391.154463_g777_i0.p3"|
                                           NAME == "NODE_240_length_20051_cov_536.532736_g170_i0.p5"|
                                           NAME == "NODE_561_length_14978_cov_299.386582_g399_i0.p6"|
                                           NAME == "NODE_5592_length_5844_cov_58.440998_g3721_i0.p3"|
                                           NAME == "NODE_2371_length_8622_cov_428.524857_g1633_i0.p4"|
                                           NAME == "NODE_874_length_12891_cov_268.735372_g614_i0.p5"|
                                           NAME == "NODE_12148_length_3772_cov_90.048391_g8052_i0.p1"|
                                           NAME == "NODE_83_length_26285_cov_124.652335_g63_i0.p2"|
                                           NAME == "NODE_1258_length_11267_cov_26.880918_g886_i0.p1"|
                                           NAME == "NODE_23095_length_2406_cov_9.263609_g15805_i0.p1"|
                                           NAME == "NODE_4_length_44673_cov_539.121996_g3_i0.p3"|
                                           NAME == "NODE_3235_length_7553_cov_85.053877_g2210_i0.p1"|
                                           NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p1"|
                                           NAME == "NODE_831_length_13112_cov_408.719610_g583_i0.p1"|
                                           NAME == "NODE_21553_length_2535_cov_63.720146_g14668_i0.p1"|
                                           NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p3"|
                                           NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p2"|
                                           NAME == "NODE_921_length_12697_cov_270.632129_g478_i1.p3"|
                                           NAME == "NODE_3303_length_7477_cov_63.941113_g2252_i0.p5"|
                                           NAME == "NODE_1672_length_10023_cov_110.507035_g1165_i0.p4"|
                                           NAME == "NODE_2027_length_9227_cov_211.794625_g1411_i0.p3"|
                                           NAME == "NODE_941_length_12600_cov_131.586493_g658_i0.p3"|
                                           NAME == "NODE_2100_length_9113_cov_67.578650_g1457_i0.p4"|
                                           NAME == "NODE_1640_length_10104_cov_99.995813_g1142_i0.p2"|
                                           NAME == "NODE_119_length_24550_cov_65.461454_g91_i0.p2"))


ple_ST$NAME[ple_ST$NAME == "NODE_356_length_17472_cov_223.329444_g259_i0.p3"] <- "GlsCD1"
ple_ST$NAME[ple_ST$NAME == "NODE_1112_length_11791_cov_391.154463_g777_i0.p3"] <- "GlsCD2"
ple_ST$NAME[ple_ST$NAME == "NODE_240_length_20051_cov_536.532736_g170_i0.p5"] <- "GlsCD3"
ple_ST$NAME[ple_ST$NAME == "NODE_561_length_14978_cov_299.386582_g399_i0.p6"] <- "GlsP1"
ple_ST$NAME[ple_ST$NAME == "NODE_5592_length_5844_cov_58.440998_g3721_i0.p3"] <- "GlsP2"
ple_ST$NAME[ple_ST$NAME == "NODE_2371_length_8622_cov_428.524857_g1633_i0.p4"] <- "GlsP3"
ple_ST$NAME[ple_ST$NAME == "NODE_874_length_12891_cov_268.735372_g614_i0.p5"] <- "GlsQ"
ple_ST$NAME[ple_ST$NAME == "NODE_12148_length_3772_cov_90.048391_g8052_i0.p1"] <- "SbtA"
ple_ST$NAME[ple_ST$NAME == "NODE_83_length_26285_cov_124.652335_g63_i0.p2"] <- "BicA1"
ple_ST$NAME[ple_ST$NAME == "NODE_1258_length_11267_cov_26.880918_g886_i0.p1"] <- "BicA2"
ple_ST$NAME[ple_ST$NAME == "NODE_23095_length_2406_cov_9.263609_g15805_i0.p1"] <- "BicA3"
ple_ST$NAME[ple_ST$NAME == "NODE_4_length_44673_cov_539.121996_g3_i0.p3"] <- "BicA4"
ple_ST$NAME[ple_ST$NAME == "NODE_3235_length_7553_cov_85.053877_g2210_i0.p1"] <- "BicA5"
ple_ST$NAME[ple_ST$NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p1"] <- "BicA6"
ple_ST$NAME[ple_ST$NAME == "NODE_831_length_13112_cov_408.719610_g583_i0.p1"] <- "NdhF1"
ple_ST$NAME[ple_ST$NAME == "NODE_21553_length_2535_cov_63.720146_g14668_i0.p1"] <- "NdhF2"
ple_ST$NAME[ple_ST$NAME == "NODE_190_length_21684_cov_465.159410_g142_i0.p3"] <- "NdhF3"
ple_ST$NAME[ple_ST$NAME == "NODE_423_length_16370_cov_3813.623182_g313_i0.p2"] <- "NdhF4"
ple_ST$NAME[ple_ST$NAME == "NODE_921_length_12697_cov_270.632129_g478_i1.p3"] <- "SBP"
ple_ST$NAME[ple_ST$NAME == "NODE_3303_length_7477_cov_63.941113_g2252_i0.p5"] <- "TMDS"
ple_ST$NAME[ple_ST$NAME == "NODE_1672_length_10023_cov_110.507035_g1165_i0.p4"] <- "TMDL"
ple_ST$NAME[ple_ST$NAME == "NODE_2027_length_9227_cov_211.794625_g1411_i0.p3"] <- "Ggps.1"
ple_ST$NAME[ple_ST$NAME == "NODE_941_length_12600_cov_131.586493_g658_i0.p3"] <- "Ggps.2"
ple_ST$NAME[ple_ST$NAME == "NODE_2100_length_9113_cov_67.578650_g1457_i0.p4"] <- "Ggps.3"
ple_ST$NAME[ple_ST$NAME == "NODE_1640_length_10104_cov_99.995813_g1142_i0.p2"] <- "HK_pyrG"
ple_ST$NAME[ple_ST$NAME == "NODE_119_length_24550_cov_65.461454_g91_i0.p2"] <- "HK_rpoB"

ple_ST <- tibble::column_to_rownames(ple_ST, "NAME")
ple_ST_t<-t(ple_ST)
ple_ST_t<-data.frame(ple_ST_t)
ple_ST_t$condition <- vec
ple_ST_t$condition1 <- vec1
ple_ST_t$condition2 <- vec2

ple_ST_tg<-ple_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
ple_ST_tg.summary <- ple_ST_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

ple_ST_plot<-ggplot(ple_ST_tg.summary, aes(x = Gene, y = len)) +
  # annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf,
  #           alpha = .1,fill = "black")+
  geom_col(aes(fill = condition), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = len-sd, ymax = len+sd, group = condition),
    width = 0.2, position = position_dodge(0.8))+
   ylim(0,100) +
  scale_fill_manual(values = c("#696969", "#C0C0C0")) +
  facet_wrap(facets = "condition1",ncol = 2) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank(),  
        aspect.ratio = 0.5) +
#  scale_x_discrete(limits = c("GlsCD1","GlsCD2","GlsCD3","GlsP1","GlsP2","GlsP3","GlsR","GlsQ","","SbtA","","BicA1","BicA2","BicA3","BicA4","BicA5","BicA6","","NdhF1","NdhF2","NdhF3","NdhF4","","SBP","TMDS","TMDL"))+
  scale_x_discrete(limits = c("GlsCD1","GlsCD2","GlsCD3","GlsP1","GlsP2","GlsP3","GlsR","GlsQ","","Ggps.1","Ggps.2","Ggps.3","HK_pyrG","HK_rpoB"))+
   xlab("") +
  ylab("Pleurocapsales transcripts per million (TPM)")


#NOSTOCALES NITROGEN GENES
nos_NT<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_10694_length_4102_cov_4262.831472_g7079_i0.p3"|
                                           NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p2"|
                                           NAME == "NODE_3972_length_6895_cov_886.344767_g2681_i0.p1"|
                                           NAME == "NODE_361_length_17230_cov_65.716559_g264_i0.p5"|
                                           NAME == "NODE_917_length_12716_cov_949.146089_g643_i0.p3"))

nos_NT$NAME[nos_NT$NAME == "NODE_10694_length_4102_cov_4262.831472_g7079_i0.p3"] <- "NifH"
nos_NT$NAME[nos_NT$NAME == "NODE_8701_length_4636_cov_426.174227_g5774_i0.p2"] <- "HetR"
nos_NT$NAME[nos_NT$NAME == "NODE_3972_length_6895_cov_886.344767_g2681_i0.p1"] <- "Amt4.1"
nos_NT$NAME[nos_NT$NAME == "NODE_361_length_17230_cov_65.716559_g264_i0.p5"] <- "Amt4.2"
nos_NT$NAME[nos_NT$NAME == "NODE_917_length_12716_cov_949.146089_g643_i0.p3"] <- "NrtP"

nos_NT <- tibble::column_to_rownames(nos_NT, "NAME")
nos_NT_t<-t(nos_NT)
nos_NT_t<-data.frame(nos_NT_t)
nos_NT_t$condition <- vec
nos_NT_t$condition1 <- vec1
nos_NT_t$condition2 <- vec2

nos_NT_tg<-nos_NT_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
nos_NT_tg.summary <- nos_NT_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

nos_NT_plot<-ggplot(nos_NT_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("NifH","HetR","","Amt4.1","Amt4.2","NrtP"))+
  xlab("") +
  ylab("Nostocales transcripts per million (TPM)")

#PLEUROCAPSALES NITROGEN GENES
ple_NT<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_6571_length_5392_cov_27.662907_g4379_i0.p1"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p1"|
                                           NAME == "NODE_4951_length_6207_cov_67.991686_g3311_i0.p1"|
                                           NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p4"|
                                           NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p2"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p2"|
                                           NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p3"|
                                           NAME == "NODE_3011_length_7804_cov_26.634459_g2060_i0.p5"|
                                           NAME == "NODE_16766_length_3034_cov_278.412023_g11235_i0.p1"|
                                           NAME == "NODE_35_length_31865_cov_1182.712947_g24_i0.p5"|
                                           NAME == "NODE_36114_length_1700_cov_60.636755_g25872_i0.p1"|
                                           NAME == "NODE_18172_length_2867_cov_38.164281_g12221_i0.p1"|
                                           NAME == "NODE_6790_length_5299_cov_25.533678_g4530_i0.p2"|
                                           NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p2"))

ple_NT$NAME[ple_NT$NAME == "NODE_6571_length_5392_cov_27.662907_g4379_i0.p1"] <- "UrtA.1"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p1"] <- "UrtA.2"
ple_NT$NAME[ple_NT$NAME == "NODE_4951_length_6207_cov_67.991686_g3311_i0.p1"] <- "UrtA.3"
ple_NT$NAME[ple_NT$NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p4"] <- "UrtB"
ple_NT$NAME[ple_NT$NAME == "NODE_1745_length_9825_cov_77.397252_g1218_i0.p2"] <- "UrtC.1"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p2"] <- "UrtC.2"
ple_NT$NAME[ple_NT$NAME == "NODE_3694_length_7113_cov_290.764773_g2500_i0.p3"] <- "UrtD"
ple_NT$NAME[ple_NT$NAME == "NODE_3011_length_7804_cov_26.634459_g2060_i0.p5"] <- "UrtE"
ple_NT$NAME[ple_NT$NAME == "NODE_16766_length_3034_cov_278.412023_g11235_i0.p1"] <- "Amt4.1"
ple_NT$NAME[ple_NT$NAME == "NODE_35_length_31865_cov_1182.712947_g24_i0.p5"] <- "Amt4.2"
ple_NT$NAME[ple_NT$NAME == "NODE_36114_length_1700_cov_60.636755_g25872_i0.p1"] <- "Amt4.3"
ple_NT$NAME[ple_NT$NAME == "NODE_18172_length_2867_cov_38.164281_g12221_i0.p1"] <- "Amt4.4"
ple_NT$NAME[ple_NT$NAME == "NODE_6790_length_5299_cov_25.533678_g4530_i0.p2"] <- "NrtP.1"
ple_NT$NAME[ple_NT$NAME == "NODE_1783_length_9722_cov_371.517567_g1243_i0.p2"] <- "NrtP.2"

ple_NT <- tibble::column_to_rownames(ple_NT, "NAME")
ple_NT_t<-t(ple_NT)
ple_NT_t<-data.frame(ple_NT_t)
ple_NT_t$condition <- vec
ple_NT_t$condition1 <- vec1
ple_NT_t$condition2 <- vec2

ple_NT_tg<-ple_NT_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
ple_NT_tg.summary <- ple_NT_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

ple_NT_plot<-ggplot(ple_NT_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("UrtA.1","UrtA.2","UrtA.3","UrtB","UrtC.1","UrtC.2","UrtD","UrtE","","Amt4.1","Amt4.2","Amt4.3","Amt4.4","","NrtP.1","NrtP.2"))+
  xlab("") +
  ylab("Pleurocapsales transcripts per million (TPM)")
plot_grid(nos_PS_plot,ple_PS_plot,nos_EPS_plot,ple_EPS_plot,nos_ST_plot,ple_ST_plot,ncol=2,align = "hv")



#ASCOMYCETE GENERAL
asc_ST<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1"|
                                           NAME == "NODE_14384_length_3379_cov_78.120387_g7996_i1.p1"|
                                           NAME == "NODE_4567_length_6455_cov_422.087120_g3060_i0.p1"|
                                           NAME == "NODE_13676_length_3494_cov_536.718211_g9088_i0.p1"|
                                           NAME == "NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1"|
                                           NAME == "NODE_16183_length_3114_cov_396.395594_g8006_i2.p1"|
                                           NAME == "NODE_4962_length_6202_cov_2127.062163_g814_i2.p1"|
                                           NAME == "NODE_1044_length_12085_cov_179.367716_g730_i0.p2"|
                                           NAME == "NODE_10303_length_4194_cov_265.035671_g6823_i0.p1"|
                                           NAME == "NODE_16489_length_3071_cov_218.059039_g11031_i0.p1"|
                                           NAME == "NODE_18046_length_2881_cov_105.940527_g11297_i1.p1"|
                                           NAME == "NODE_5851_length_5712_cov_344.059940_g3319_i2.p1"|
                                           NAME == "NODE_9287_length_4466_cov_124.225131_g6160_i0.p1"|
                                           NAME == "NODE_24063_length_2331_cov_863.463685_g15307_i1.p1"|
                                           NAME == "NODE_7692_length_4969_cov_416.356209_g5101_i0.p1"|
                                           NAME == "NODE_6333_length_5497_cov_401.459440_g4201_i0.p1"|
                                           NAME == "NODE_4815_length_6302_cov_223.018141_g3217_i0.p1"|
                                           NAME == "NODE_5199_length_6066_cov_338.280828_g3004_i1.p1"|
                                           NAME == "NODE_6016_length_5641_cov_295.210489_g3993_i0.p1"|
                                           NAME == "NODE_4787_length_6317_cov_625.977578_g3200_i0.p1"|
                                           NAME == "NODE_13641_length_3500_cov_130.069157_g9063_i0.p1"|
                                           NAME == "NODE_7461_length_5052_cov_108.702752_g4955_i0.p1"|
                                           NAME == "NODE_1092_length_11881_cov_92.157859_g761_i0.p4"|
                                           NAME == "NODE_32715_length_1841_cov_194.979638_g22800_i1.p1"|
                                           NAME == "NODE_75155_length_933_cov_183.155814_g1593_i4.p1"|
                                           NAME == "NODE_29556_length_1994_cov_216.781884_g20722_i0.p1"))

asc_ST$NAME[asc_ST$NAME == "NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1"] <- "ST.1"
asc_ST$NAME[asc_ST$NAME == "NODE_14384_length_3379_cov_78.120387_g7996_i1.p1"] <- "ST.2"
asc_ST$NAME[asc_ST$NAME == "NODE_4567_length_6455_cov_422.087120_g3060_i0.p1"] <- "ST.3"
asc_ST$NAME[asc_ST$NAME == "NODE_13676_length_3494_cov_536.718211_g9088_i0.p1"] <- "ST.4"
asc_ST$NAME[asc_ST$NAME == "NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1"] <- "ST.5"
asc_ST$NAME[asc_ST$NAME == "NODE_16183_length_3114_cov_396.395594_g8006_i2.p1"] <- "ST.6"
asc_ST$NAME[asc_ST$NAME == "NODE_4962_length_6202_cov_2127.062163_g814_i2.p1"] <- "ST.7"
asc_ST$NAME[asc_ST$NAME == "NODE_1044_length_12085_cov_179.367716_g730_i0.p2"] <- "ST.8"
asc_ST$NAME[asc_ST$NAME == "NODE_10303_length_4194_cov_265.035671_g6823_i0.p1"] <- "UreA"
asc_ST$NAME[asc_ST$NAME == "NODE_16489_length_3071_cov_218.059039_g11031_i0.p1"] <- "AAP.1"
asc_ST$NAME[asc_ST$NAME == "NODE_18046_length_2881_cov_105.940527_g11297_i1.p1"] <- "AAP.2"
asc_ST$NAME[asc_ST$NAME == "NODE_5851_length_5712_cov_344.059940_g3319_i2.p1"] <- "AAP.3"
asc_ST$NAME[asc_ST$NAME == "NODE_9287_length_4466_cov_124.225131_g6160_i0.p1"] <- "AAP.4"
asc_ST$NAME[asc_ST$NAME == "NODE_24063_length_2331_cov_863.463685_g15307_i1.p1"] <- "Gln.synt"
asc_ST$NAME[asc_ST$NAME == "NODE_7692_length_4969_cov_416.356209_g5101_i0.p1"] <- "Glu.synt"
asc_ST$NAME[asc_ST$NAME == "NODE_6333_length_5497_cov_401.459440_g4201_i0.p1"] <- "ABC.1"
asc_ST$NAME[asc_ST$NAME == "NODE_4815_length_6302_cov_223.018141_g3217_i0.p1"] <- "ABC.2"
asc_ST$NAME[asc_ST$NAME == "NODE_5199_length_6066_cov_338.280828_g3004_i1.p1"] <- "ABC.3"
asc_ST$NAME[asc_ST$NAME == "NODE_6016_length_5641_cov_295.210489_g3993_i0.p1"] <- "ABC.4"
asc_ST$NAME[asc_ST$NAME == "NODE_4787_length_6317_cov_625.977578_g3200_i0.p1"] <- "ABC.5"
asc_ST$NAME[asc_ST$NAME == "NODE_13641_length_3500_cov_130.069157_g9063_i0.p1"] <- "ABC.6"
asc_ST$NAME[asc_ST$NAME == "NODE_7461_length_5052_cov_108.702752_g4955_i0.p1"] <- "ABC.7"
asc_ST$NAME[asc_ST$NAME == "NODE_1092_length_11881_cov_92.157859_g761_i0.p4"] <- "ABC.8"
asc_ST$NAME[asc_ST$NAME == "NODE_32715_length_1841_cov_194.979638_g22800_i1.p1"] <- "ABC.9"
asc_ST$NAME[asc_ST$NAME == "NODE_75155_length_933_cov_183.155814_g1593_i4.p1"] <- "MFS_2.1"
asc_ST$NAME[asc_ST$NAME == "NODE_29556_length_1994_cov_216.781884_g20722_i0.p1"] <- "MFS_2.2"

asc_ST <- tibble::column_to_rownames(asc_ST, "NAME")
asc_ST_t<-t(asc_ST)
asc_ST_t<-data.frame(asc_ST_t)
asc_ST_t$condition <- vec
asc_ST_t$condition1 <- vec1
asc_ST_t$condition2 <- vec2

asc_ST_tg<-asc_ST_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_ST_tg.summary <- asc_ST_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

asc_ST_plot<-ggplot(asc_ST_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("ST.1","ST.2","ST.3","ST.4","ST.5","ST.6","ST.7","ST.8","","UreA","","AAP.1","AAP.2","AAP.3","AAP.4","","Gln.synt","Glu.synt","","ABC.1","ABC.2","ABC.3","ABC.4","ABC.5","ABC.6","ABC.7","ABC.8","ABC.9","MFS_2.1","MFS_2.2"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")

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
                                           NAME == "NODE_48399_length_1340_cov_403.876875_g36078_i0.p1"|
                                           NAME == "NODE_39896_length_1565_cov_1941.174263_g24942_i1.p1"))

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
asc_GH$NAME[asc_GH$NAME == "NODE_39896_length_1565_cov_1941.174263_g24942_i1.p1"] <- "HK_tubB"


asc_GH <- tibble::column_to_rownames(asc_GH, "NAME")
asc_GH_t<-t(asc_GH)
asc_GH_t<-data.frame(asc_GH_t)
asc_GH_t$condition <- vec
asc_GH_t$condition1 <- vec1
asc_GH_t$condition2 <- vec2

asc_GH_tg<-asc_GH_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_GH_tg.summary <- asc_GH_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

asc_GH_plot<-ggplot(asc_GH_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("GH3.1","GH3.2","GH3.3","","GH5.1","GH5.2","","GH16.1","GH16.2","GH16.3","GH16.4","","GH17.1","GH17.2","GH17.3","GH17.4","","GH18.1","GH18.2","","GH20.1","GH20.2","","GH31.1","GH31.2","","GH32","","GH37","","GH47.1","GH47.2","GH47.3","GH47.4","GH47.5","","GH63","","GH65","","GH81","","GH115","","GH128.1","HK_tubB"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")
#ASCOMYCETE SOD
asc_SOD<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_10663_length_4109_cov_869.988850_g6777_i1.p4"|
                                            NAME ==  "NODE_119255_length_650_cov_1242.746967_g77078_i1.p1"|
                                            NAME == "NODE_34001_length_1784_cov_1363.392168_g24170_i0.p1"|
                                            NAME == "NODE_82198_length_868_cov_2955.817610_g66734_i0.p1"|
                                            NAME == "NODE_3331_length_7450_cov_473.596448_g2273_i0.p3"|
                                            NAME == "NODE_21563_length_2534_cov_2580.070703_g3154_i1.p1"|
                                            NAME == "NODE_11789_length_3846_cov_562.072886_g7815_i0.p1"|
                                            NAME == "NODE_19014_length_2775_cov_682.914508_g10783_i1.p1"|
                                            NAME == "NODE_89258_length_812_cov_4283.902571_g48896_i1.p2"|
                                            NAME == "NODE_2258_length_8831_cov_102.924640_g1559_i0.p1"|
                                            NAME == "NODE_22523_length_2452_cov_246.747373_g15374_i0.p1"|
                                            NAME == "NODE_4714_length_6365_cov_819.396694_g3154_i0.p2"|
                                            NAME == "NODE_26538_length_2165_cov_772.933556_g18411_i0.p1"))

asc_SOD$NAME[asc_SOD$NAME == "NODE_10663_length_4109_cov_869.988850_g6777_i1.p4"] <- "CuSOD"
asc_SOD$NAME[asc_SOD$NAME == "NODE_119255_length_650_cov_1242.746967_g77078_i1.p1"] <- "MnSOD"
asc_SOD$NAME[asc_SOD$NAME == "NODE_34001_length_1784_cov_1363.392168_g24170_i0.p1"] <- "Redoxin1"
asc_SOD$NAME[asc_SOD$NAME == "NODE_82198_length_868_cov_2955.817610_g66734_i0.p1"] <- "Redoxin2"
asc_SOD$NAME[asc_SOD$NAME == "NODE_3331_length_7450_cov_473.596448_g2273_i0.p3"] <- "Thioredoxin1"
asc_SOD$NAME[asc_SOD$NAME == "NODE_21563_length_2534_cov_2580.070703_g3154_i1.p1"] <- "Thioredoxin2"
asc_SOD$NAME[asc_SOD$NAME == "NODE_11789_length_3846_cov_562.072886_g7815_i0.p1"] <- "Thioredoxin3"
asc_SOD$NAME[asc_SOD$NAME == "NODE_19014_length_2775_cov_682.914508_g10783_i1.p1"] <- "Thioredoxin4"
asc_SOD$NAME[asc_SOD$NAME == "NODE_89258_length_812_cov_4283.902571_g48896_i1.p2"] <- "Thioredoxin5"
asc_SOD$NAME[asc_SOD$NAME == "NODE_2258_length_8831_cov_102.924640_g1559_i0.p1"] <- "Thioredoxin6"
asc_SOD$NAME[asc_SOD$NAME == "NODE_22523_length_2452_cov_246.747373_g15374_i0.p1"] <- "Thioredoxin7"
asc_SOD$NAME[asc_SOD$NAME == "NODE_4714_length_6365_cov_819.396694_g3154_i0.p2"] <- "Thioredoxin8"
asc_SOD$NAME[asc_SOD$NAME == "NODE_26538_length_2165_cov_772.933556_g18411_i0.p1"] <- "Hog1"


asc_SOD <- tibble::column_to_rownames(asc_SOD, "NAME")
asc_SOD_t<-t(asc_SOD)
asc_SOD_t<-data.frame(asc_SOD_t)
asc_SOD_t$condition <- vec
asc_SOD_t$condition1 <- vec1
asc_SOD_t$condition2 <- vec2

asc_SOD_tg<-asc_SOD_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_SOD_tg.summary <- asc_SOD_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

asc_SOD_plot<-ggplot(asc_SOD_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("CuSOD","MnSOD","Redoxin1","Redoxin2","Thioredoxin1","Thioredoxin2","Thioredoxin3","Thioredoxin4","Thioredoxin5","Thioredoxin6","Thioredoxin7","Thioredoxin8","Hog1"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")


#ASCOMYCETE Trehalose
asc_tre<-asc_diff_ex_READS_df %>% filter((NAME == "NODE_35370_length_1729_cov_433.846014_g1044_i3.p1"|
                                          NAME == "NODE_48677_length_1334_cov_124.367962_g36314_i0.p1"|
                                          NAME == "NODE_1489_length_10513_cov_660.778640_g1044_i0.p2"|
                                          NAME == "NODE_11285_length_3958_cov_200.838095_g7478_i0.p2"))



asc_tre$NAME[asc_tre$NAME == "NODE_35370_length_1729_cov_433.846014_g1044_i3.p1"] <- "GT20.1"
asc_tre$NAME[asc_tre$NAME == "NODE_48677_length_1334_cov_124.367962_g36314_i0.p1"] <- "GT20.2"
asc_tre$NAME[asc_tre$NAME == "NODE_1489_length_10513_cov_660.778640_g1044_i0.p2"] <- "GT20.3_Trehalose_PPase"
asc_tre$NAME[asc_tre$NAME == "NODE_11285_length_3958_cov_200.838095_g7478_i0.p2"] <- "Trehalose_PPase"

asc_tre <- tibble::column_to_rownames(asc_tre, "NAME")
asc_tre_t<-t(asc_tre)
asc_tre_t<-data.frame(asc_tre_t)
asc_tre_t$condition <- vec
asc_tre_t$condition1 <- vec1
asc_tre_t$condition2 <- vec2

asc_tre_tg<-asc_tre_t%>% gather(Gene, TPM, -condition, -condition1, -condition2)
asc_tre_tg.summary <- asc_tre_tg %>% group_by(condition2,condition1,condition,Gene) %>% summarise( sd = sd(TPM, na.rm = TRUE),len = mean(TPM))

asc_tre_plot<-ggplot(asc_tre_tg.summary, aes(x = Gene, y = len)) +
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
  scale_x_discrete(limits = c("GT20.1","GT20.2","GT20.3_Trehalose_PPase","Trehalose_PPase"))+
  xlab("") +
  ylab("Ascomycota transcripts per million (TPM)")



library(tidyverse)
library(rstatix)
library(ggpubr)
asc_ST_tg_DRY<-subset(asc_ST_tg, condition1=="DRY")
asc_ST_tg_DRY$logTPM<-log(asc_ST_tg_DRY$TPM+1)
stat.test <- asc_ST_tg_DRY %>%
  group_by(Gene) %>%
  t_test(logTPM ~ condition2, paired=FALSE) %>%
#  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
stat_test_df<-data.frame(stat.test)
# Create the plot
myplot  <-ggboxplot(
  asc_ST_tg_DRY, x = "condition2", y = "logTPM",
  fill = "condition1", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~Gene, ncol=11)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "condition2")
myplot + stat_pvalue_manual(stat.test, label = "p.signif")

