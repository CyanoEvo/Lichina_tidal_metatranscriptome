library(tidyverse)
library(janitor)

nos_diff_ex_READS_df <- tibble::rownames_to_column(nos_diff_ex_TPM, "NAME")
nos_EPS<-nos_diff_ex_READS_df %>% filter((NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p1"|
                                          NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p2"|
                                          NAME == "NODE_4760_length_6331_cov_526.782518_g3184_i0.p3"|
                                          NAME == "NODE_636_length_14310_cov_263.887968_g450_i0.p5"|
                                          NAME == "NODE_14243_length_3402_cov_434.991589_g9475_i0.p1"|
                                          NAME == "NODE_13298_length_3560_cov_641.725839_g8837_i0.p1"))

nos_EPS_t<-t(nos_EPS)
nos_EPS_t<-janitor::row_to_names(nos_EPS_t, row_number = 1)
vec<-c("HT","HT","HT","LT","LT","LT","HT","HT","HT","LT","LT","LT")
vec1<-c("DRY","DRY","DRY","DRY","DRY","DRY","WET","WET","WET","WET","WET","WET")

nos_EPS_t<-data.frame(nos_EPS_t)
nos_EPS_t$condition <- vec
nos_EPS_t$condition1 <- vec1

nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p1 = as.numeric(as.character(nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p1))
nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p2 = as.numeric(as.character(nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p2))
nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p3 = as.numeric(as.character(nos_EPS_t$NODE_4760_length_6331_cov_526.782518_g3184_i0.p3))
nos_EPS_t$NODE_636_length_14310_cov_263.887968_g450_i0.p5 = as.numeric(as.character(nos_EPS_t$NODE_636_length_14310_cov_263.887968_g450_i0.p5))
nos_EPS_t$NODE_14243_length_3402_cov_434.991589_g9475_i0.p1 = as.numeric(as.character(nos_EPS_t$NODE_14243_length_3402_cov_434.991589_g9475_i0.p1))
nos_EPS_t$NODE_13298_length_3560_cov_641.725839_g8837_i0.p1 = as.numeric(as.character(nos_EPS_t$NODE_13298_length_3560_cov_641.725839_g8837_i0.p1))



Wza_nos<-ggplot(nos_EPS_t, aes(x=condition, y=)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wza (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzc_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_4760_length_6331_cov_526.782518_g3184_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzc (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzy_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_4760_length_6331_cov_526.782518_g3184_i0.p2)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzy (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzx_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_4760_length_6331_cov_526.782518_g3184_i0.p3)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzx (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzx1_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_636_length_14310_cov_263.887968_g450_i0.p5)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzx1 (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

KpsD_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_14243_length_3402_cov_434.991589_g9475_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("KpsD (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
KpsE_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_13298_length_3560_cov_641.725839_g8837_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("KpsD (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
plot_grid(Wza_nos,Wzc_nos,Wzy_nos,Wzx_nos,Wzx1_nos,KpsD_nos,Kps_nos,ncol=2,align = "hv")

#############PLEUROCAPSALES###########

ple_diff_ex_READS_df <- tibble::rownames_to_column(ple_diff_ex_TPM, "NAME")
ple_EPS<-ple_diff_ex_READS_df %>% filter((NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p9"|
                                            NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p4"|
                                            NAME == "NODE_355_length_17525_cov_279.467683_g258_i0.p1"|
                                            NAME == "NODE_88_length_26018_cov_220.834496_g39_i1.p6"))

NODE_88_length_26018_cov_220.834496_g39_i1.p6

ple_EPS_t<-t(ple_EPS)
ple_EPS_t<-janitor::row_to_names(ple_EPS_t, row_number = 1)
vec<-c("HT","HT","HT","LT","LT","LT","HT","HT","HT","LT","LT","LT")
vec1<-c("DRY","DRY","DRY","DRY","DRY","DRY","WET","WET","WET","WET","WET","WET")

ple_EPS_t<-data.frame(ple_EPS_t)
ple_EPS_t$condition <- vec
ple_EPS_t$condition1 <- vec1

ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p9 = as.numeric(as.character(ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p9))
ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p4 = as.numeric(as.character(ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p4))
ple_EPS_t$NODE_355_length_17525_cov_279.467683_g258_i0.p1 = as.numeric(as.character(ple_EPS_t$NODE_355_length_17525_cov_279.467683_g258_i0.p1))
ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p6 = as.numeric(as.character(ple_EPS_t$NODE_88_length_26018_cov_220.834496_g39_i1.p6))

Wza_ple<-ggplot(ple_EPS_t, aes(x=condition, y=NODE_88_length_26018_cov_220.834496_g39_i1.p9)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wza (Pleurocapsales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzc_ple<-ggplot(ple_EPS_t, aes(x=condition, y=NODE_88_length_26018_cov_220.834496_g39_i1.p4)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzc (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzc1_ple<-ggplot(ple_EPS_t, aes(x=condition, y=NODE_355_length_17525_cov_279.467683_g258_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzc (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzy_ple<-ggplot(ple_EPS_t, aes(x=condition, y=NODE_88_length_26018_cov_220.834496_g39_i1.p6)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzy (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzx_ple<-ggplot(ple_EPS_t, aes(x=condition, y=)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzx (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

Wzx1_nos<-ggplot(nos_EPS_t, aes(x=condition, y=NODE_636_length_14310_cov_263.887968_g450_i0.p5)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("Wzx1 (Nostocales)") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
plot_grid(Wza_nos,Wzc_nos,Wzy_nos,Wzx_nos,Wzx1_nos,ncol=2,align = "hv")
 
#############Ascomycota###########

asc_diff_ex_READS_df <- tibble::rownames_to_column(asc_diff_ex_TPM, "NAME")
asc_SUGTRANS<-asc_diff_ex_READS_df %>% filter(( NAME == "NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1"|
                                                NAME == "NODE_13676_length_3494_cov_536.718211_g9088_i0.p1"|
                                                NAME == "NODE_4567_length_6455_cov_422.087120_g3060_i0.p1"|
                                                NAME == "NODE_16941_length_3012_cov_2397.016672_g11361_i0.p1"|
                                                NAME == "NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1"|
                                                NAME == "NODE_4962_length_6202_cov_2127.062163_g814_i2.p1"|
                                                NAME == "NODE_1044_length_12085_cov_179.367716_g730_i0.p2"|
                                                NAME == "NODE_16183_length_3114_cov_396.395594_g8006_i2.p1"|
                                                NAME == "NODE_12068_length_3785_cov_285.955280_g8006_i0.p1"|
                                                NAME == "NODE_80321_length_884_cov_390.374846_g64964_i0.p1"|
                                                NAME == "NODE_21108_length_2575_cov_123.826938_g13270_i1.p1"|
                                                NAME == "NODE_313074_length_347_cov_231.536496_g295578_i0.p1"|
                                                NAME == "NODE_21371_length_2552_cov_128.323518_g14522_i0.p1"|
                                                NAME == "NODE_76594_length_919_cov_253.602837_g44486_i1.p1"|
                                                NAME == "NODE_10702_length_4099_cov_286.372578_g7083_i0.p1"|
                                                NAME == "NODE_8352_length_4748_cov_225.893904_g5423_i1.p1"|
                                                NAME == "NODE_7660_length_4983_cov_158.582485_g4506_i1.p1"|
                                                NAME == "NODE_22029_length_2494_cov_2721.621231_g10229_i1.p1"|
                                                NAME == "NODE_107641_length_702_cov_2746.181240_g91076_i0.p1"|
                                                NAME == "NODE_162648_length_522_cov_388.020045_g145353_i0.p1"|
                                                NAME == "NODE_96072_length_766_cov_407.738817_g66067_i1.p1"|
                                                NAME == "NODE_3338_length_7444_cov_313.755257_g2279_i0.p1"))
asc_SUGTRANS_t<-t(asc_SUGTRANS)
asc_SUGTRANS_t<-janitor::row_to_names(asc_SUGTRANS_t, row_number = 1)

asc_SUGTRANS_t<-data.frame(asc_SUGTRANS_t)
asc_SUGTRANS_t$condition <- vec
asc_SUGTRANS_t$condition1 <- vec1

asc_SUGTRANS_t$NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1))
asc_SUGTRANS_t$NODE_13676_length_3494_cov_536.718211_g9088_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_13676_length_3494_cov_536.718211_g9088_i0.p1))
asc_SUGTRANS_t$NODE_4567_length_6455_cov_422.087120_g3060_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_4567_length_6455_cov_422.087120_g3060_i0.p1))
asc_SUGTRANS_t$NODE_16941_length_3012_cov_2397.016672_g11361_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_16941_length_3012_cov_2397.016672_g11361_i0.p1))
asc_SUGTRANS_t$NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1))
asc_SUGTRANS_t$NODE_4962_length_6202_cov_2127.062163_g814_i2.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_4962_length_6202_cov_2127.062163_g814_i2.p1))
asc_SUGTRANS_t$NODE_1044_length_12085_cov_179.367716_g730_i0.p2 = as.numeric(as.character(asc_SUGTRANS_t$NODE_1044_length_12085_cov_179.367716_g730_i0.p2))
asc_SUGTRANS_t$NODE_16183_length_3114_cov_396.395594_g8006_i2.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_16183_length_3114_cov_396.395594_g8006_i2.p1))
asc_SUGTRANS_t$NODE_12068_length_3785_cov_285.955280_g8006_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_12068_length_3785_cov_285.955280_g8006_i0.p1))
asc_SUGTRANS_t$NODE_80321_length_884_cov_390.374846_g64964_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_80321_length_884_cov_390.374846_g64964_i0.p1))
asc_SUGTRANS_t$NODE_21108_length_2575_cov_123.826938_g13270_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_21108_length_2575_cov_123.826938_g13270_i1.p1))
asc_SUGTRANS_t$NODE_313074_length_347_cov_231.536496_g295578_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_313074_length_347_cov_231.536496_g295578_i0.p1))
asc_SUGTRANS_t$NODE_21371_length_2552_cov_128.323518_g14522_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_21371_length_2552_cov_128.323518_g14522_i0.p1))
asc_SUGTRANS_t$NODE_76594_length_919_cov_253.602837_g44486_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_76594_length_919_cov_253.602837_g44486_i1.p1))
asc_SUGTRANS_t$NODE_10702_length_4099_cov_286.372578_g7083_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_10702_length_4099_cov_286.372578_g7083_i0.p1))
asc_SUGTRANS_t$NODE_8352_length_4748_cov_225.893904_g5423_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_8352_length_4748_cov_225.893904_g5423_i1.p1))
asc_SUGTRANS_t$NODE_7660_length_4983_cov_158.582485_g4506_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_7660_length_4983_cov_158.582485_g4506_i1.p1))
asc_SUGTRANS_t$NODE_22029_length_2494_cov_2721.621231_g10229_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_22029_length_2494_cov_2721.621231_g10229_i1.p1))
asc_SUGTRANS_t$NODE_107641_length_702_cov_2746.181240_g91076_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_107641_length_702_cov_2746.181240_g91076_i0.p1))
asc_SUGTRANS_t$NODE_162648_length_522_cov_388.020045_g145353_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_162648_length_522_cov_388.020045_g145353_i0.p1))
asc_SUGTRANS_t$NODE_96072_length_766_cov_407.738817_g66067_i1.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_96072_length_766_cov_407.738817_g66067_i1.p1))
asc_SUGTRANS_t$NODE_3338_length_7444_cov_313.755257_g2279_i0.p1 = as.numeric(as.character(asc_SUGTRANS_t$NODE_3338_length_7444_cov_313.755257_g2279_i0.p1))


asc_SUGTRANS_p1<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_16379_length_3087_cov_4368.462840_g10957_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)       
asc_SUGTRANS_p2<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_13676_length_3494_cov_536.718211_g9088_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)   
asc_SUGTRANS_p3<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_4567_length_6455_cov_422.087120_g3060_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1) 

asc_SUGTRANS_p4<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_16941_length_3012_cov_2397.016672_g11361_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p5<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_17453_length_2951_cov_2614.985754_g11361_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p6<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_4962_length_6202_cov_2127.062163_g814_i2.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p7<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_1044_length_12085_cov_179.367716_g730_i0.p2)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p8<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_16183_length_3114_cov_396.395594_g8006_i2.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p9<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_12068_length_3785_cov_285.955280_g8006_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p10<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_80321_length_884_cov_390.374846_g64964_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p11<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_21108_length_2575_cov_123.826938_g13270_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p12<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_313074_length_347_cov_231.536496_g295578_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p13<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_21371_length_2552_cov_128.323518_g14522_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p14<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_76594_length_919_cov_253.602837_g44486_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p15<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_10702_length_4099_cov_286.372578_g7083_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p16<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_8352_length_4748_cov_225.893904_g5423_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p17<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_7660_length_4983_cov_158.582485_g4506_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p18<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_22029_length_2494_cov_2721.621231_g10229_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p19<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_107641_length_702_cov_2746.181240_g91076_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p20<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_162648_length_522_cov_388.020045_g145353_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

asc_SUGTRANS_p21<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_96072_length_766_cov_407.738817_g66067_i1.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
asc_SUGTRANS_p22<-ggplot(asc_SUGTRANS_t, aes(x=condition, y=NODE_3338_length_7444_cov_313.755257_g2279_i0.p1)) + 
  geom_boxplot(fill="red") +
  facet_wrap(~condition1) + ylab("") + xlab(element_blank()) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
plot_grid(asc_SUGTRANS_p1,asc_SUGTRANS_p2,asc_SUGTRANS_p3,asc_SUGTRANS_p4,asc_SUGTRANS_p5,asc_SUGTRANS_p6,asc_SUGTRANS_p7,asc_SUGTRANS_p8,asc_SUGTRANS_p9,asc_SUGTRANS_p10,asc_SUGTRANS_p11,asc_SUGTRANS_p12,asc_SUGTRANS_p13,asc_SUGTRANS_p14,asc_SUGTRANS_p15,asc_SUGTRANS_p16,asc_SUGTRANS_p17,asc_SUGTRANS_p18,asc_SUGTRANS_p19,asc_SUGTRANS_p20,asc_SUGTRANS_p21,asc_SUGTRANS_p22,ncol=4,align = "hv")
