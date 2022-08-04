# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/CHIP/", showWarnings = FALSE, recursive = TRUE)


#LHP1 -----

#Data import
LHP1_output_qPCR <-
  read_excel("data/CHIP/LHP1_output_qPCR.xlsx") %>%
  clean_genotype_code %>%
  mutate(Cp_IGG=as.numeric(Cp_IGG), Cp_LHP1=as.numeric(Cp_LHP1))

#Data processing
LHP1_output_qPCR_PercentInput <- LHP1_output_qPCR %>%
  mutate(Percent_IGG = 2^-(Cp_IGG-Cp_input_IGG)*100,
         Percent_LHP1= 2^-(Cp_LHP1-Cp_input_LHP1)*100,
         LHP1background=Percent_LHP1/Percent_IGG) %>%
  filter(!is.na(LHP1background)) %>%
  mutate(time=as.character(time))


# Statistics and Graph
LHP1_output_qPCR_PercentInput$Position <-
    factor(LHP1_output_qPCR_PercentInput$Position,
           c("CYP705A12.3'","CYP705A12.5'","CYP71A16.3'","CYP71A16.5'",
             "intergenic1","intergenic2","intergenic3","intergenic4", "MHAL.5'","MHAL.3'","MRN1.promotor", "MRN1.5'","MRN1.3'"))
LHP1_output_qPCR_PercentInput$time <-
  factor(LHP1_output_qPCR_PercentInput$time,
         c("0","4 hours"))

#Col t0 vs t4 ----
# Anova
my_graph <- 
LHP1_output_qPCR_PercentInput %>% filter(genotype =="col") %>%
ggbarplot(., x = "Position", y = "LHP1background", add = "mean_se",
          fill = "time", 
          position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("LHP1 background") +
  geom_hline(yintercept=0.5, linetype=2)+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_LHP1_col t0 vs t4.pdf", plot=my_graph, width = 5, height = 3)

#Col vs RNAi t0 ----
# Anova

my_graph <- 
  LHP1_output_qPCR_PercentInput %>% filter(time =="0") %>%
  ggbarplot(., x = "Position", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("LHP1 background") +
  geom_hline(yintercept=0.5,linetype=2) +
  scale_fill_manual(labels=genotype_metadata$genotype_name,
                    breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,)+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_LHP1_col vs RNAi t0.pdf", plot=my_graph, width = 5, height = 3)

#Graph/position ChIP LHP1 ----

LHP1_output_qPCR_PercentInput_intergenic2 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="intergenic2") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic2, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 36) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_intergenic2.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_intergenic4 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="intergenic4") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic4, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 36) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_intergenic4.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1pP <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="MRN1.promotor") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1pP, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 28) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_MRN1pP.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p5 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="MRN1.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p5, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 8) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_MRN1p5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p3 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="MRN1.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p3, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 1.5) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_MRN1p3.pdf", plot=my_graph, width = 2.5, height = 3)
#
#
LHP1_output_qPCR_PercentInput_intergenic1 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="intergenic1") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic1, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 2.5) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_intergenic3 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="intergenic3") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic3, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 27) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_intergenic3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_intergenic4 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="intergenic4") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic4, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 15) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_intergenic4.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MHAL5 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="MHAL.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL5, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 15) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_MHAL5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MHAL3 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="MHAL.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL3, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 12) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_MHAL3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_3 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="CYP71A16.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_3, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 1.2) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_CYP71_3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_5 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="CYP71A16.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_5, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 4) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_CYP71_5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP705_3 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="CYP705A12.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_3, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 4) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_CYP705_3.pdf", plot=my_graph, width = 2.5, height = 3)

#
LHP1_output_qPCR_PercentInput_CYP705_5 <- LHP1_output_qPCR_PercentInput %>%
  filter(Position =="CYP705A12.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_5, x = "time", y = "LHP1background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 2) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_CYP705_5.pdf", plot=my_graph, width = 2.5, height = 3)




#CHIP H3 and H3K27me3 -----
#Data import
H3H3K27me3_output_qPCR <- 
  read_excel("data/CHIP/H3H3K27me3_output_qPCR.xlsx") %>%
  clean_genotype_code %>%
  mutate(Cp=as.numeric(Cp))


#Data processing
qPCR_H3_input <- H3H3K27me3_output_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_H3_IP_input <- H3H3K27me3_output_qPCR %>% 
  filter(condition != "input") %>%
  left_join(qPCR_H3_input)


IP_input_IGG <- qPCR_H3_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_H3 <- qPCR_H3_IP_input %>% 
  filter(condition %in% c("H3")) %>%
  mutate(Percent_input_H3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_H3H3K27me3 <- qPCR_H3_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_H3) %>%
  left_join(IP_input_H3H3K27me3)


IP_RealMethylation_all <- Percent_input_all %>%
  mutate(H3_background=Percent_input_H3/Percent_input_IGG,
         H3K27me3_background=Percent_input_H3K27me3/Percent_input_IGG,
         RealMethylation=H3K27me3_background/H3_background,
         time=as.character(time)) %>%
  dplyr::select(genotype,time,replicate,Position,RealMethylation,H3K27me3_background)



IP_RealMethylation_all$time <- factor(IP_RealMethylation_all$time, c("0","4 hours"))
IP_RealMethylation_all$Position <-
    factor(IP_RealMethylation_all$Position,
           c("CYP705A12.3'","CYP705A12.5'","CYP71A16.3'","CYP71A16.5'",
             "intergenic1","intergenic2", "intergenic3","intergenic4", "MHAL.5'","MHAL.3'","MRN1.promotor",
             "MRN1.5'","MRN1.exon9", "MRN1.3'"))



#Col t0 vs t4 ----

my_graph <- 
    IP_RealMethylation_all %>%
    filter(genotype =="col",
           Position %in% c("CYP705A12.3'","CYP705A12.5'","CYP71A16.3'",
                           "CYP71A16.5'", "intergenic1","intergenic2", "intergenic3","intergenic4", "MHAL.5'",
                           "MHAL.3'","MRN1.5'", "MRN1.3'")) %>% 
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y=0.075) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("(H3K27me3 / H3) \n / background") +theme(legend.position="none") +
  ylim(0,100)
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_H3K27me3_col t0 vs t4.pdf", plot=my_graph, width = 5, height = 3)

#Col vs RNAi t0 ----

my_graph <- 
  IP_RealMethylation_all %>%
  filter(time =="0",
         Position %in% c("CYP705A12.3'","CYP705A12.5'","CYP71A16.3'",
                         "CYP71A16.5'", "intergenic1","intergenic2", "intergenic3","intergenic4", "MHAL.5'",
                         "MHAL.3'","MRN1.5'", "MRN1.3'")) %>% 
  ggbarplot(., x = "Position", y = "RealMethylation", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y=0.09) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("(H3K27me3 / H3) \n / background") +
  scale_fill_manual(labels=genotype_metadata$genotype_name,
                    breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,)+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_H3K27me3_col vs RNAi t0.pdf", plot=my_graph, width = 5, height = 3)


# Chip H3K27me3 on clf mutant -----
H3K27me3_clf_output_qPCR <-
    read_excel("data/CHIP/H3K27me3_clf_output_qPCR.xlsx") %>%
    clean_genotype_code

#Data processing
qPCR_clf_input <- H3K27me3_clf_output_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  select(-Cp, -condition)

qPCR_clf_IP_input <- H3K27me3_clf_output_qPCR %>% 
  filter(condition != "input") %>%
  left_join(qPCR_clf_input)


IP_input_IGG <- qPCR_clf_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  select(-Cp, -Cp_input,-condition)
IP_input_H3K27me3 <- qPCR_clf_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100)%>%
  select(-Cp, -Cp_input,-condition)


Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_H3K27me3)


IP_RealMethylation_all <- Percent_input_all %>%
  mutate(H3K27me3_background=Percent_input_H3K27me3/Percent_input_IGG,
         time=as.character(time)) %>%
  select(genotype,time,replicate,Position,H3K27me3_background) %>%
  mutate(genotype=as.factor(genotype))

IP_RealMethylation_all_t0 <- IP_RealMethylation_all %>%
  filter(time == "0") %>% 
  rename(H3K27me3_t0_background=H3K27me3_background) %>%
  select(-time,-replicate)
IP_RealMethylation_all_t4 <- IP_RealMethylation_all %>%
  filter(time == "4 hours") %>% 
  rename(H3K27me3_t4_background=H3K27me3_background) %>%
  select(-time,-replicate)

IP_RealMethylation_ratio_time <- IP_RealMethylation_all_t0 %>%
  left_join(IP_RealMethylation_all_t4) %>%
  mutate(ratio_time=H3K27me3_t0_background/H3K27me3_t4_background) %>%
  mutate(genotype=as.factor(genotype))
  

# Statistics and Graph
IP_RealMethylation_all$Position <- factor(IP_RealMethylation_all$Position, c("MHAL.5'","MHAL.3'","MRN1.5'","MRN1.3'"))
IP_RealMethylation_ratio_time$Position <- factor(IP_RealMethylation_ratio_time$Position, c("MHAL.5'","MHAL.3'","MRN1.5'","MRN1.3'"))
IP_RealMethylation_all$genotype <- factor(IP_RealMethylation_all$genotype, c("col","clf"))
IP_RealMethylation_ratio_time$genotype <- factor(IP_RealMethylation_ratio_time$genotype, c("col","clf"))
#col t0 vs t4 ----
# Anova
my_graph <- 
  IP_RealMethylation_all %>% filter(genotype =="col") %>%
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  #stat_compare_means(aes(group = time), method="anova", label = "p.signif") +
  theme_bw()+
  xlab("") +
  ylab("H3K27me3 / background")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+ theme(legend.position = "none")
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_clf_H3K27me3_col t0 vs t4.pdf", plot=my_graph, width = 5, height = 5)

my_graph <- 
  IP_RealMethylation_all %>% filter(genotype =="clf") %>%
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.signif") +
  theme_bw()+
  xlab("") +
  ylab("H3K27me3 / background")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+ theme(legend.position = "none")
my_graph




# Save 
ggsave(filename="out/CHIP/CHIP_clf_H3K27me3_clf t0 vs t4.pdf", plot=my_graph, width = 5, height = 5)

my_graph <- 
  IP_RealMethylation_all %>% filter(time =="0") %>% left_join(genotype_metadata) %>%
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype_name", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype_name), method="anova", label = "p.signif") +
  theme_bw()+
  xlab("") +
  ylab("H3K27me3 / background")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  scale_fill_manual(labels=genotype_metadata$genotype_name,
                    breaks=genotype_metadata$genotype_name,
                    values=genotype_metadata$genotype_color,)+ theme(legend.position = "none")
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_clf_H3K27me3_col vs clf t0.pdf", plot=my_graph, width = 5, height = 5)


#Ratio time ----
#Anova
my_graph <- 
  IP_RealMethylation_ratio_time %>%
  ggbarplot(., x = "Position", y = "ratio_time", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.signif", label.y=55) +
  theme_bw()+
  xlab("") +
  ylab("time 0 / time 4 \n (H3K27me3 / background)")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1)) +
  scale_fill_manual(labels=genotype_metadata$genotype_name,
                    breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,)
my_graph

# Save 
ggsave(filename="out/CHIP/CHIP_clf_H3K27me3_ratio t0 t4.pdf", plot=my_graph, width = 5, height = 3)


# ChIP H3K27me3 only ----

#Data import
H3H3K27me3_output_qPCR <- 
  read_excel("data/CHIP/H3H3K27me3_output_qPCR.xlsx", sheet=2) %>%
  clean_genotype_code %>%
  mutate(Cp=as.numeric(Cp))


#Data processing
qPCR_H3K27me3_input <- H3H3K27me3_output_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_H3K27me3_IP_input <- H3H3K27me3_output_qPCR %>% 
  filter(condition != "input") %>%
  left_join(qPCR_H3K27me3_input)


IP_input_IGG <- qPCR_H3K27me3_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3H3K27me3 <- qPCR_H3K27me3_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_H3H3K27me3) 


IP_H3K27me3 <- Percent_input_all %>%
  mutate(H3K27me3_background=Percent_input_H3K27me3/Percent_input_IGG,
         time=as.character(time)) %>%
  dplyr::select(genotype,time,replicate,Position,H3K27me3_background)


IP_H3K27me3 %>% filter(time =="0") %>%
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  xlab("") +
  ylab("H3K27me3 / background")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.25, size=10, hjust=1)) +
  ylim(0,500)


IP_H3K27me3 %>% filter(genotype == "col") %>%
  ggbarplot(., x = "Position", y = "H3K27me3_background", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  #stat_compare_means(aes(group = time), method="anova", label = "p.signif") +
  theme_bw()+
  xlab("") +
  ylab("H3K27me3 / background")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.25, size=10, hjust=1)) 


#Graph/position ChIP H3K27me3 ----

IP_H3K27me3_intergenic2 <- IP_H3K27me3 %>%
  filter(Position =="intergenic2") 

my_graph <- 
  ggbarplot(IP_H3K27me3_intergenic2, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_intergenic2.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_intergenic4 <- IP_H3K27me3 %>%
  filter(Position =="intergenic4") 

my_graph <- 
  ggbarplot(IP_H3K27me3_intergenic4, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 10) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_intergenic4.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_MRN1_pP <- IP_H3K27me3 %>%
  filter(Position =="MRN1.promotor") 

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1_pP, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 140) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_MRN1_pP.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_MRN1_p5 <- IP_H3K27me3 %>%
  filter(Position =="MRN1.5'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1_p5, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 32) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_MRN1_p5.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_MRN1_p3 <- IP_H3K27me3 %>%
  filter(Position =="MRN1.3'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1_p3, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 29) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_MRN1_p3.pdf", plot=my_graph, width = 2.5, height = 3)
#
#
IP_H3K27me3_intergenic1 <- IP_H3K27me3 %>%
  filter(Position =="intergenic1") 

my_graph <- 
  ggbarplot(IP_H3K27me3_intergenic1, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 40) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#
#
IP_H3K27me3_intergenic3 <- IP_H3K27me3 %>%
  filter(Position =="intergenic3") 

my_graph <- 
  ggbarplot(IP_H3K27me3_intergenic3, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_intergenic3.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_MHAL5 <- IP_H3K27me3 %>%
  filter(Position =="MHAL.5'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_MHAL5, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 10) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_MHAL5.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_MHAL3 <- IP_H3K27me3 %>%
  filter(Position =="MHAL.3'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_MHAL3, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 10) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_MHAL3.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_CYP71_3 <- IP_H3K27me3 %>%
  filter(Position =="CYP71A16.3'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_CYP71_3, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 22) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_CYP71_3.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_CYP71_5 <- IP_H3K27me3 %>%
  filter(Position =="CYP71A16.5'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_CYP71_5, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 40) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_CYP71_5.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_CYP705_3 <- IP_H3K27me3 %>%
  filter(Position =="CYP705A12.3'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_CYP705_3, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_CYP705_3.pdf", plot=my_graph, width = 2.5, height = 3)
#
IP_H3K27me3_CYP705_5 <- IP_H3K27me3 %>%
  filter(Position =="CYP705A12.5'") 

my_graph <- 
  ggbarplot(IP_H3K27me3_CYP705_5, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_H3K27me3_CYP705_5.pdf", plot=my_graph, width = 2.5, height = 3)




####  IN PROCESSS do not delete

output_LHP1_replicate <- read_excel("C:/Users/roule/OneDrive/Bureau/Covid19 Work/Paper_mhal/Post-reunion1/ChIP_LHP1/output_LHP1_replicate.xlsx",sheet=5) %>%
  dplyr::select(-Position) %>%
  mutate(Cp=as.numeric(Cp))

qPCR_LHP1_input <- output_LHP1_replicate %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

output_LHP1_ChIP_process <- output_LHP1_replicate %>% 
  filter(condition != "input") %>% 
  left_join(qPCR_LHP1_input)

write.csv(output_LHP1_ChIP_process, file="output_LHP1_ChIP_process")

IP_input_IGG <- output_LHP1_ChIP_process  %>%
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

IP_input_LHP1 <- output_LHP1_ChIP_process %>% 
  filter(condition %in% c("LHP1")) %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

LHP1_output_qPCR_PercentInput_retour <- IP_input_IGG %>%
  left_join(IP_input_LHP1) %>%
  mutate(LHP1background=Percent_input_LHP1/Percent_input_IGG,
         time=as.character(time), 
         replicate=as.character(replicate))


LHP1_output_qPCR_PercentInput_retour %>% filter(genotype =="col") %>%
  ggplot(., aes(time, LHP1background)) +
  geom_point(aes(color=replicate))
  
LHP1_output_qPCR_PercentInput_retour %>% filter(time =="0") %>%
  ggplot(., aes(genotype, LHP1background)) +
  geom_point(aes(color=replicate))




# LHP1 in vitro (the sheet 2 is the one used and published!!! number 5 is the corrected one with the repetition with GFP added) -----
output_qPCR_inVitro <- read_excel("data/CHIP/output_qPCR_inVitro.xlsx", sheet=5)%>%
  filter(Cp>1)


qPCR_LHP1_IP <- output_qPCR_inVitro %>%
  filter(condition %in% c("IGG", "LHP1"))


qPCR_LHP1_input <- output_qPCR_inVitro %>%
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_LHP1_IP_input <- qPCR_LHP1_IP %>%
  left_join(qPCR_LHP1_input)


qPCR_LHP1_IP_input_IGGpercent <- qPCR_LHP1_IP_input %>%
  filter(condition == "IGG") %>%
  mutate(Percent_input = 2^-(Cp-Cp_input)*100)

qPCR_LHP1_IP_input_LHP1percent <- qPCR_LHP1_IP_input %>%
  filter(condition == "LHP1") %>%
  mutate(Percent_input = 2^-(Cp-Cp_input)*100)


qPCR_LHP1_Percent_input <- qPCR_LHP1_IP_input_IGGpercent %>%
  bind_rows(qPCR_LHP1_IP_input_LHP1percent)


qPCR_LHP1_ratio_IGG <- qPCR_LHP1_Percent_input %>%
  dplyr::select(-Cp, -Cp_input) %>%
  filter(condition =="IGG") %>%
  dplyr::select(-condition) %>%
  rename(Percent_IGG=Percent_input)

qPCR_LHP1_ratio_LHP1 <- qPCR_LHP1_Percent_input %>%
  dplyr::select(-Cp, -Cp_input) %>%
  filter(condition =="LHP1") %>%
  dplyr::select(-condition) %>%
  rename(Percent_LHP1=Percent_input)

qPCR_LHP1_ratio <- qPCR_LHP1_ratio_IGG %>%
  left_join(qPCR_LHP1_ratio_LHP1) %>%
  mutate(ratio=Percent_LHP1/Percent_IGG)


ggplot(qPCR_LHP1_ratio, aes(x=concentration, y=ratio)) +
  geom_boxplot(aes(color=position))


stat_CHIP_invitro <- qPCR_LHP1_ratio %>%
  mutate(concentration=as.character(concentration)) %>%
  group_by(concentration, position, RNA) %>%
  summarise(mean=mean(ratio),
            median= median(ratio),
            SD=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillon
            erreur_std=SD/sqrt(n))


#For statistic
my_graph <-   
  qPCR_LHP1_ratio %>% filter(position %in% c("MRN1.promoter")) %>% 
  ggbarplot(., x = "RNA", y = "ratio",add = "mean_se") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("GFP","MARS") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  xlab("") + 
  theme(legend.position="none") +
  scale_fill_grey() +
  facet_wrap(~concentration)

my_graph

ggsave(filename="out/CHIP/statistic_MRN1promoter.pdf", plot=my_graph, width = 3, height = 2)

#



my_graph <- 
  stat_CHIP_invitro %>% filter(position %in% c("MRN1.promoter")) %>%
  ggplot(data=., mapping=aes(x=concentration, y=mean, group=RNA)) +
  geom_bar(stat="identity", position=position, aes(fill=RNA),colour="black") +
  geom_errorbar(mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5),position=position)+
  xlab(label = "")+
  ylab(label = "")+
  geom_hline(yintercept=1, linetype=2) +
  scale_fill_manual(values=selected_gene_invitro$gene_color,)

my_graph



ggsave(filename="out/CHIP/ChIP_invitro_MRN1_promoter.pdf", plot=my_graph, width = 7, height = 4)


position=position_dodge(0.9)


my_graph <- 
  stat_CHIP_invitro %>%
  ggplot(data=., mapping=aes(x=concentration, y=mean, group=position)) +
  geom_bar(stat="identity", position=position, aes(fill=position)) +
  geom_errorbar(mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5),position=position)+
  xlab(label = "")+
  ylab(label = "")+
  facet_wrap(~RNA)+
  geom_hline(yintercept=1, linetype=2)

my_graph
ggsave(filename="out/CHIP/ChIP_invitro_all.pdf", plot=my_graph, width = 6, height = 3)



























# ChIP Line2 H3K27me3 and LHP1 -------

H3K27me3_LHP1_Line2 <- read_excel("data/CHIP/H3K27me3_LHP1_Line2.xlsx")


qPCR_IP <- H3K27me3_LHP1_Line2 %>%
  filter(condition %in% c("IGG", "LHP1", "K27"))


qPCR_input <- H3K27me3_LHP1_Line2 %>%
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_IP_input <- qPCR_IP %>%
  left_join(qPCR_input) %>%
  mutate(time=as.character(time),
         replicate=as.character(replicate))


qPCR_IP_input_IGGpercent <- qPCR_IP_input %>%
  filter(condition == "IGG") %>%
  mutate(Percent_input = 2^-(Cp-Cp_input)*100) %>%
  rename(Percent_input_IGG=Percent_input) %>%
  dplyr::select(-condition, -Cp, -Cp_input)

qPCR_IP_input_LHP1percent <- qPCR_IP_input %>%
  filter(condition %in% c("LHP1")) %>%
  mutate(Percent_input = 2^-(Cp-Cp_input)*100) %>%
  rename(Percent_input_LHP1=Percent_input) %>%
  dplyr::select(-condition, -Cp, -Cp_input)

qPCR_IP_input_K27percent <- qPCR_IP_input %>%
  filter(condition %in% c("K27")) %>%
  mutate(Percent_input = 2^-(Cp-Cp_input)*100) %>%
  rename(Percent_input_K27=Percent_input) %>%
  dplyr::select(-condition, -Cp, -Cp_input)





qPCR_Percent_input <- qPCR_IP_input_IGGpercent %>%
  left_join(qPCR_IP_input_LHP1percent) %>%
  left_join(qPCR_IP_input_K27percent) %>%
  mutate(Enrichment_H3K27me3 = Percent_input_K27 / Percent_input_IGG,
         Enrichment_LHP1 = Percent_input_LHP1 / Percent_input_IGG)



stat_CHIP_H3K27me3 <- qPCR_Percent_input %>%
  group_by(gene, time,genotype) %>%
  summarise(mean=mean(Enrichment_H3K27me3),
            median= median(Enrichment_H3K27me3),
            SD=sd(Enrichment_H3K27me3), #ecart-type
            n=n(), #nombre d'?chantillon
            erreur_std=SD/sqrt(n))

stat_CHIP_LHP1 <- qPCR_Percent_input %>%
  group_by(gene, time,genotype) %>%
  summarise(mean=mean(Enrichment_LHP1),
            median= median(Enrichment_LHP1),
            SD=sd(Enrichment_LHP1), #ecart-type
            n=n(), #nombre d'?chantillon
            erreur_std=SD/sqrt(n))



stat_CHIP_H3K27me3 %>%
  ggplot(., aes(x=time, y=mean, group=genotype)) +
  geom_bar(stat = "identity", position="dodge", aes(fill=genotype)) +
  theme_bw() +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std, group=gene), position = position_dodge(width = 0.9), width=.5)+
  facet_wrap(~gene, scale="free")

stat_CHIP_LHP1 %>%
  ggplot(., aes(x=time, y=mean, group=genotype)) +
  geom_bar(stat = "identity", position="dodge", aes(fill=genotype)) +
  theme_bw() +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std, group=gene), position = position_dodge(width = 0.9), width=.5)+
  facet_wrap(~gene, scale="free")

#Graph/position ChIP K27 RNAi3 ----
qPCR_Percent_input$time <-
  factor(qPCR_Percent_input$time,
         c("0","4 hours"))


LHP1_output_qPCR_PercentInput_intergenic2 <- qPCR_Percent_input %>%
  filter(gene =="intergenic2") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic2, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_intergenic2.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_intergenic1 <- qPCR_Percent_input %>%
  filter(gene =="intergenic1") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic1, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1pP <- qPCR_Percent_input %>%
  filter(gene =="MRN1.promoter") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1pP, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_MRN1pP.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p5 <- qPCR_Percent_input %>%
  filter(gene =="MRN1.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p5, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_MRN1p5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p3 <- qPCR_Percent_input %>%
  filter(gene =="MRN1.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p3, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_MRN1p3.pdf", plot=my_graph, width = 2.5, height = 3)
#
#
LHP1_output_qPCR_PercentInput_intergenic1 <- qPCR_Percent_input %>%
  filter(gene =="intergenic1") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic1, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#

#
LHP1_output_qPCR_PercentInput_MHAL5 <- qPCR_Percent_input %>%
  filter(gene =="MHAL.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL5, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_MHAL5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MHAL3 <- qPCR_Percent_input %>%
  filter(gene =="MHAL.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL3, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_MHAL3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_3 <- qPCR_Percent_input %>%
  filter(gene =="CYP71A16.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_3, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_CYP71_3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_5 <- qPCR_Percent_input %>%
  filter(gene =="CYP71A16.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_5, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_CYP71_5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP705_3 <- qPCR_Percent_input %>%
  filter(gene =="CYP705A12.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_3, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_CYP705_3.pdf", plot=my_graph, width = 2.5, height = 3)

#
LHP1_output_qPCR_PercentInput_CYP705_5 <- qPCR_Percent_input %>%
  filter(gene =="CYP705A12.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_5, x = "time", y = "Enrichment_H3K27me3", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_K27_Line2_CYP705_5.pdf", plot=my_graph, width = 2.5, height = 3)


#Graph/position ChIP LHP1 RNAi3 ----
qPCR_Percent_input$time <-
  factor(qPCR_Percent_input$time,
         c("0","4 hours"))


LHP1_output_qPCR_PercentInput_intergenic2 <- qPCR_Percent_input %>%
  filter(gene =="intergenic2") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic2, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_intergenic2.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_intergenic1 <- qPCR_Percent_input %>%
  filter(gene =="intergenic1") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic1, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,)  + theme(legend.position = "none")
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1pP <- qPCR_Percent_input %>%
  filter(gene =="MRN1.promoter") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1pP, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_MRN1pP.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p5 <- qPCR_Percent_input %>%
  filter(gene =="MRN1.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p5, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_MRN1p5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MRN1p3 <- qPCR_Percent_input %>%
  filter(gene =="MRN1.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MRN1p3, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_MRN1p3.pdf", plot=my_graph, width = 2.5, height = 3)
#
#
LHP1_output_qPCR_PercentInput_intergenic1 <- qPCR_Percent_input %>%
  filter(gene =="intergenic1") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_intergenic1, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_intergenic1.pdf", plot=my_graph, width = 2.5, height = 3)
#

#
LHP1_output_qPCR_PercentInput_MHAL5 <- qPCR_Percent_input %>%
  filter(gene =="MHAL.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL5, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_MHAL5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_MHAL3 <- qPCR_Percent_input %>%
  filter(gene =="MHAL.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_MHAL3, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_MHAL3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_3 <- qPCR_Percent_input %>%
  filter(gene =="CYP71A16.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_3, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_CYP71_3.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP71_5 <- qPCR_Percent_input %>%
  filter(gene =="CYP71A16.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP71_5, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_CYP71_5.pdf", plot=my_graph, width = 2.5, height = 3)
#
LHP1_output_qPCR_PercentInput_CYP705_3 <- qPCR_Percent_input %>%
  filter(gene =="CYP705A12.3'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_3, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_CYP705_3.pdf", plot=my_graph, width = 2.5, height = 3)

#
LHP1_output_qPCR_PercentInput_CYP705_5 <- qPCR_Percent_input %>%
  filter(gene =="CYP705A12.5'") 

my_graph <- 
  ggbarplot(LHP1_output_qPCR_PercentInput_CYP705_5, x = "time", y = "Enrichment_LHP1", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph

ggsave(filename="out/CHIP/ChIP_LHP1_Line2_CYP705_5.pdf", plot=my_graph, width = 2.5, height = 3)



# Ying -------
ying_qpcr_output <- read_excel("data/CHIP/Ying/ying_qpcr_output.xlsx",sheet=1) %>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  unique()%>%
  ungroup()


chip_repet_1 <- read_excel("data/CHIP/Ying/chip repet 1.xlsx",sheet=2)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_2 <- read_excel("data/CHIP/Ying/chip repet 2.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()


chip_repet_3 <- read_excel("data/CHIP/Ying/chip repet 3.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_4 <- read_excel("data/CHIP/Ying/chip repet 4.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_5 <- read_excel("data/CHIP/Ying/chip repet 5.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_6 <- read_excel("data/CHIP/Ying/chip repet 6.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_7 <- read_excel("data/CHIP/Ying/chip repet 7.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_8 <- read_excel("data/CHIP/Ying/chip repet 8.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_9 <- read_excel("data/CHIP/Ying/chip repet 9.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_10 <- read_excel("data/CHIP/Ying/chip repet 10.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_repet_11 <- read_excel("data/CHIP/Ying/chip repet 11.xlsx",sheet=1)%>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na() %>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  select(-Rtech)%>%
  ungroup()%>%
  unique()


chip_repet <- ying_qpcr_output %>%
  bind_rows(chip_repet_1)%>%
  bind_rows(chip_repet_2)%>%
  bind_rows(chip_repet_4)%>%
  bind_rows(chip_repet_5)%>%
  bind_rows(chip_repet_6)%>%
  bind_rows(chip_repet_7)%>%
  bind_rows(chip_repet_8)%>%
  bind_rows(chip_repet_9)%>%
  bind_rows(chip_repet_10)%>%
  bind_rows(chip_repet_11)
  

#Data processing
qPCR_ying_input <- chip_repet %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_ying_input_IP_input <- chip_repet %>% 
  filter(condition != "input") %>%
  left_join(qPCR_ying_input)


IP_input_IGG <- qPCR_ying_input_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3 <- qPCR_ying_input_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_LHP1 <- qPCR_ying_input_IP_input %>% 
  filter(condition %in% c("LHP1")) %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_H3K27me3 %>%
  left_join(IP_input_LHP1) %>%
  left_join(IP_input_IGG)


IP_all <- Percent_input_all %>%
  mutate(H3K27me3_background=Percent_input_H3K27me3/Percent_input_IGG,
         LHP1_background=Percent_input_LHP1/Percent_input_IGG) %>%
  select(genotype,time,replicate,H3K27me3_background,LHP1_background,gene)

IP_all$time <- factor(IP_all$time, c("0", "4 hours"))


#Graph/position ChIP H3K27me3 ----

IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.promoter")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




ggsave(filename="out/CHIP/ChIP_H3K27me3_intergenic2.pdf", plot=my_graph, width = 2.5, height = 3)
#



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="enhancer2")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="enhancer1")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph





IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph






IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP705A12.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP705A12.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph


IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP71A16.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP71A16.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MARS.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MARS.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




# Repetition TR/TB/AC ------------

chip_1 <- read_excel("data/CHIP/Repet_all_TR/chip.1.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_2 <- read_excel("data/CHIP/Repet_all_TR/chip.2.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_3 <- read_excel("data/CHIP/Repet_all_TR/chip.3.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_4 <- read_excel("data/CHIP/Repet_all_TR/chip.4.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_5 <- read_excel("data/CHIP/Repet_all_TR/chip.5.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_6 <- read_excel("data/CHIP/Repet_all_TR/chip.6.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_7 <- read_excel("data/CHIP/Repet_all_TR/chip.7.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_8 <- read_excel("data/CHIP/Repet_all_TR/chip.8.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_9 <- read_excel("data/CHIP/Repet_all_TR/chip.9.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_10 <- read_excel("data/CHIP/Repet_all_TR/chip.10.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()
chip_11 <- read_excel("data/CHIP/Repet_all_TR/chip.11.xlsx")%>%
  dplyr::select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))%>%
  group_by(time,genotype,replicate,condition,gene)%>%
  mutate(Cp=mean(Cp)) %>%
  dplyr::select(-Rtech)%>%
  ungroup()%>%
  unique()

chip_all_TR <- chip_1 %>%
  bind_rows(chip_2)%>%
  bind_rows(chip_3)%>%
  bind_rows(chip_4)%>%
  bind_rows(chip_5)%>%
  bind_rows(chip_6)%>%
  bind_rows(chip_7)%>%
  bind_rows(chip_8)%>%
  bind_rows(chip_9)%>%
  bind_rows(chip_10)%>%
  bind_rows(chip_11)





chip_all_TR %>%
ggplot(., aes(time,Cp,shape=condition,color=genotype))+
  geom_jitter(width = 0.5)+
  facet_wrap(~gene)

chip_all_TR <- read_excel("data/CHIP/Repet_all_TR/chip_all_TR.xlsx")




#Data processing
qPCR_TR_input <- chip_all_TR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_TR_input_IP_input <- chip_all_TR %>% 
  filter(condition != "input") %>%
  left_join(qPCR_TR_input)


IP_input_IGG <- qPCR_TR_input_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3 <- qPCR_TR_input_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_LHP1 <- qPCR_TR_input_IP_input %>% 
  filter(condition %in% c("LHP1")) %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_H3K27me3 %>%
  left_join(IP_input_LHP1) %>%
  left_join(IP_input_IGG)


IP_all <- Percent_input_all %>%
  mutate(H3K27me3_background=Percent_input_H3K27me3/Percent_input_IGG,
         LHP1_background=Percent_input_LHP1/Percent_input_IGG) %>%
  dplyr::select(genotype,time,replicate,H3K27me3_background,LHP1_background,gene)

IP_all$time <- factor(IP_all$time, c("0", "4"))


#some stats
IP_all_stat_H3K27me3 <- IP_all %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(H3K27me3_background),	 		
            median= median(H3K27me3_background),			
            SD=sd(H3K27me3_background),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			
IP_all_stat_LHP1 <- IP_all %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(LHP1_background),	 		
            median= median(LHP1_background),			
            SD=sd(LHP1_background),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


IP_input_H3K27me3_tidy <- IP_input_H3K27me3 %>%
  pivot_longer(Percent_input_H3K27me3)
IP_input_LHP1_tidy <- IP_input_LHP1 %>%
  pivot_longer(Percent_input_LHP1)
IP_input_IGG_tidy <- IP_input_IGG %>%
  pivot_longer(Percent_input_IGG)

IP_tidy <- IP_input_H3K27me3_tidy %>%
  bind_rows(IP_input_LHP1_tidy) %>%
  bind_rows(IP_input_IGG_tidy)
IP_tidy_stat <- IP_tidy %>%			
  group_by(gene, time, genotype,name) %>%	
  summarise(mean=mean(value),	 		
            median= median(value),			
            SD=sd(value),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

IP_tidy_stat %>%
  mutate(time=as.character(time))%>%
ggplot(., aes(time,mean,fill=name))+
  geom_bar(stat="identity", position=position_dodge(.9)) +
  geom_errorbar(mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+
  facet_grid(genotype~gene)






#Graph/position ChIP H3K27me3 ----

IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.promoter")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MRN1.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="intergenic1")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="intergenic2")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP705A12.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP705A12.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph





IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP71A16.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph




IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="CYP71A16.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph



IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MHAL.3'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph


IP_H3K27me3_MRN1pro <- IP_all %>%
  filter(gene=="MHAL.5'")

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "H3K27me3_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

my_graph <- 
  ggbarplot(IP_H3K27me3_MRN1pro, x = "time", y = "LHP1_background", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=1,linetype=2)+
  facet_wrap(~gene)
my_graph

