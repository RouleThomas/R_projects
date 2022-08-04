

# ChIP H3K27me3 only ----

#Data import
H3H3K27me3_output_qPCR <- 
  read_excel("H3K27me3_output_qPCR.xlsx") %>%
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






#LHP1 -----

#Data import
LHP1_output_qPCR <-
  read_excel("data/ChIP/LHP1_output_qPCR.xlsx") %>%
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

my_graph <- 
  LHP1_output_qPCR_PercentInput %>% filter(genotype =="col",
                                           time == 0,
                                           Position %in% c("intergenic1","intergenic2","MRN1.promotor","MRN1.5'","MRN1.3'")) %>%
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



ggsave(filename="out/CHIP/ChIP_LHP1_thesis.pdf", plot=my_graph, width = 5, height = 4)


# LHP1 in vitro -----
output_qPCR_inVitro <- read_excel("output_qPCR_inVitro_.xlsx")%>%
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
  qPCR_LHP1_ratio %>% filter(position %in% c("MRN1.3'")) %>% 
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

H3K27me3_LHP1_Line2 <- read_excel("H3K27me3_LHP1_Line2_.xlsx")


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


# Ying -----
ying_qpcr_output <- read_excel("data/CHIP/Ying/ying_qpcr_output.xlsx")





