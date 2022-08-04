


# ChIP arabido ------------

#Data import
H3K27me3_ChIP_qPCR <- read_excel("ChIP/in/H3K27me3_ChIP_qPCR.xlsx")


#Data processing
qPCR_H3K27me3_input <- H3K27me3_ChIP_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_H3K27me3_IP_input <- H3K27me3_ChIP_qPCR %>% 
  filter(condition != "input") %>%
  left_join(qPCR_H3K27me3_input)


IP_input_IGG <- qPCR_H3K27me3_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3_good <- qPCR_H3K27me3_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3 <- qPCR_H3K27me3_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_H3H3K27me3) 


IP_H3K27me3 <- Percent_input_all %>%
  mutate(Enrichment=Percent_input_H3K27me3/Percent_input_IGG) %>%
  dplyr::select(replicate,region,Enrichment)


IP_H3K27me3 %>% 
  filter(region != "FLC-1") %>%
ggplot(., aes(region, Enrichment))+
  geom_boxplot()


# plot with stat 

IP_H3K27me3$region <- factor(IP_H3K27me3$region, levels=c("housekeeping","PID","MRN1","FLC-2"))

my_graph <- 
  IP_H3K27me3 %>%
  filter(region != "FLC-1") %>%
  ggbarplot(., x = "region", y = "Enrichment",add = "mean_se", fill="darkgrey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("housekeeping", "PID"), c("housekeeping", "MRN1"), c("housekeeping", "FLC-2")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  geom_abline(slope=0, intercept=2,  col = "red",lty=2)+
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/H3K27me3 Arabido.pdf", plot=my_graph, width = 7, height = 4) 


# plot percent input

IP_IGG <- IP_input_IGG %>% 
  pivot_longer(Percent_input_IGG)

IP_H3K27me3 <- IP_input_H3H3K27me3 %>% 
  pivot_longer(Percent_input_H3K27me3)

IP_all <- IP_IGG %>%
  bind_rows(IP_H3K27me3)

IP_all %>%
ggplot(., aes(region, value))+
  geom_boxplot(aes(fill=name))





IP_all$region <- factor(IP_all$region, levels=c("housekeeping","PID","MRN1","FLC-2"))

my_graph <- 
  IP_all %>% filter(region != "FLC-1") %>%
  ggbarplot(., x = "region", y = "value",add = "mean_se", fill="name", position=position_dodge()) +
  stat_compare_means(aes(group=name), label = "p.format", method="t.test") +
  theme_bw() 
my_graph

# Save 
ggsave(filename="ChIP/out/arabido percent input.pdf", plot=my_graph, width = 7, height = 4) 








# ChIP rice FIE ------------



#Data import

X20220309 <- read_excel("ChIP/in/20220309.xls") %>%
  dplyr::select(-Pos)


#Data processing
qPCR_H3K27me3_input <- X20220309 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_H3K27me3_IP_input <- X20220309 %>% 
  filter(condition != "input") %>%
  left_join(qPCR_H3K27me3_input)


IP_input_IGG <- qPCR_H3K27me3_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3H3K27me3 <- qPCR_H3K27me3_IP_input %>% 
  filter(condition %in% c("FIE")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_H3H3K27me3) 


IP_H3K27me3 <- Percent_input_all %>%
  mutate(Enrichment=Percent_input_H3K27me3/Percent_input_IGG) %>%
  dplyr::select(replicate,region,Enrichment)


IP_H3K27me3 %>% 
  ggplot(., aes(region, Enrichment))+
  geom_boxplot()


# plot with stat 

IP_H3K27me3$region <- factor(IP_H3K27me3$region, levels=c("housekeeping","PID","MRN1","FLC-2"))

my_graph <- 
  IP_H3K27me3 %>%
  ggbarplot(., x = "region", y = "Enrichment",add = "mean_se", fill="darkgrey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("housekeeping", "PID"), c("housekeeping", "MRN1"), c("housekeeping", "FLC-2")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  geom_abline(slope=0, intercept=2,  col = "red",lty=2)+
  ylim(0,3)+
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/FIE rice.pdf", plot=my_graph, width = 7, height = 4) 





# ChIP arabido SWI3B GFP RNAse A------------

#Data import
GFP_rChIP_qPCR_pl1 <- read_excel("ChIP/in/GFP_rChIP_qPCR_pl1.xlsx", sheet=2) 

GFP_rChIP_qPCR_pl2 <- read_excel("ChIP/in/GFP_rChIP_qPCR_pl2.xlsx", sheet=2) 

GFP_rChIP_qPCR_pl3 <- read_excel("ChIP/in/GFP_rChIP_qPCR_pl3.xlsx", sheet=2) 

GFP_rChIP_qPCR <- GFP_rChIP_qPCR_pl1 %>%
  bind_rows(GFP_rChIP_qPCR_pl2) %>%
  bind_rows(GFP_rChIP_qPCR_pl3) %>%
  select(-Pos)

#Data processing
GFP_rChIP_qPCR_input <- GFP_rChIP_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition, -RNAse_A)

GFP_rChIP_qPCR_IP_input <- GFP_rChIP_qPCR %>% 
  filter(condition != "input") %>%
  left_join(GFP_rChIP_qPCR_input)


IP_input_IGG <- GFP_rChIP_qPCR_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_GFP <- GFP_rChIP_qPCR_IP_input %>% 
  filter(condition %in% c("GFP")) %>%
  mutate(Percent_input_GFP = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_IGG %>%
  left_join(IP_input_GFP) 


IP_GFP <- Percent_input_all %>%
  mutate(Enrichment=Percent_input_GFP/Percent_input_IGG) %>%
  dplyr::select(replicate,position,Enrichment, RNAse_A)


IP_GFP %>% 
  filter(replicate !=1)%>%
  ggplot(., aes(position, Enrichment))+
  geom_boxplot(aes(fill=RNAse_A))


# plot with stat UNDER WORK -----

IP_H3K27me3$region <- factor(IP_H3K27me3$region, levels=c("housekeeping","PID","MRN1","FLC-2"))

my_graph <- 
  IP_GFP %>%
  ggbarplot(., x = "position", y = "Enrichment",add = "mean_se", fill="RNAse_A") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("housekeeping", "PID"), c("housekeeping", "MRN1"), c("housekeeping", "FLC-2")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  geom_abline(slope=0, intercept=2,  col = "red",lty=2)+
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/H3K27me3 Arabido.pdf", plot=my_graph, width = 7, height = 4) 

# ----
  
# plot percent input

IP_IGG <- IP_input_IGG %>% 
  pivot_longer(Percent_input_IGG)

IP_GFP <- IP_input_GFP %>% 
  pivot_longer(Percent_input_GFP)

IP_all <- IP_IGG %>%
  bind_rows(IP_GFP) 



IP_all %>%
  ggplot(., aes(RNAse_A, value))+
  geom_boxplot(aes(fill=percent_input))




df <- tibble(
  name = c('Percent_input_IGG', 'Percent_input_GFP')
  ,percent_input = c('IGG', 'GFP')
)

IP_all <- IP_IGG %>%
  bind_rows(IP_GFP) %>%
  left_join(df) %>%
  select(-name) %>% 
  unite(position_RNAseA, c(position, RNAse_A), sep = "_", remove = FALSE) %>%
  select(-position)

# statistic ggpubr -----

"AT2G26630_minus", "AT2G26630_plus",
"AT2G04460_minus", "AT2G04460_plus",
"AT5G19015_minus", "AT5G19015_plus"
,"ref4_minus", "ref4_plus",
"MRN1_promoter_minus", "MRN1_promoter_plus",
"Actin2_minus", "Actin2_plus",
"PID_D_minus", "PID_D_plus"
"APOLO_G_minus", "APOLO_G_plus",
"FIL_a_minus", "FIL_a_plus",
"FIL_b_minus", "FIL_b_plus",
"FIL_c_minus", "FIL_c_plus",
"IAMT1_2_minus", "IAMT1_2_plus",
"IAMT1_3_minus", "IAMT1_3_plus",
"IAMT1_4_minus", "IAMT1_4_plus",
"IAMT1_6_minus", "IAMT1_6_plus"


my_graph <- 
IP_all %>%
  ggplot(., aes(x = percent_input, y = value, color = as.character(replicate)))+
  facet_wrap(~position_RNAseA, scale="free")+
  geom_jitter()
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_replicate colored.pdf", plot=my_graph, width = 15, height = 15)


IP_all %>%
  ggbarplot(., x = "percent_input", y = "value",add = c("mean_se","jitter"), fill="percent_input", facet.by = "position_RNAseA", scales="free",ncol=4) +
  theme_bw() +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(face="bold", size=10))





IP_all %>%
  filter(position_RNAseA %in% c("APOLO_G_minus", "APOLO_G_plus"
                                )) %>%
  ggbarplot(., x = "percent_input", y = "value",add = c("mean_se"), fill="percent_input", facet.by = "position_RNAseA", scales="free",ncol=2) +
  theme_bw() +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(face="bold", size=10))



# All
my_graph <- 
IP_all %>%
  ggbarplot(., x = "percent_input", y = "value",add = "mean_se", fill="lightgrey", facet.by = "position_RNAseA", scales="free",ncol=2) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_corr.pdf", plot=my_graph, width = 5, height = 15)


# Positive (affected by RNAsA)
my_graph <- 
  IP_all %>%
  filter(replicate != 1,
         position_RNAseA %in% c("AT2G04460_minus", "AT2G04460_plus",
                                "AT2G26630_minus", "AT2G26630_plus",
                                "IAMT1_4_minus", "IAMT1_4_plus",
                                "FIL_a_minus", "FIL_a_plus")) %>%
  ggbarplot(., x = "percent_input", y = "value",add = c("mean_se","jitter"), fill="lightgrey", facet.by = "position_RNAseA", scales="free",ncol=2) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_positive.pdf", plot=my_graph, width = 5, height = 15)

# Negative not bind (not affected by RNAsA)
my_graph <- 
  IP_all %>%
  filter(replicate != 1,
         position_RNAseA %in% c("FIL_c_minus", "FIL_c_plus",
                                "IAMT1_2_minus", "IAMT1_2_plus",
                                "IAMT1_3_minus", "IAMT1_3_plus",
                                "ref4_minus", "ref4_plus",
                                "MRN1_promoter_minus", "MRN1_promoter_plus"
                                )) %>%
  ggbarplot(., x = "percent_input", y = "value",add = c("mean_se","jitter"), fill="lightgrey", facet.by = "position_RNAseA", scales="free",ncol=2) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_nobinding_notaffected.pdf", plot=my_graph, width = 5, height = 15)

# Negative bind (not affected by RNAsA)
my_graph <- 
  IP_all %>%
  filter(replicate != 1,
         position_RNAseA %in% c("Actin2_minus", "Actin2_plus",
                                "APOLO_G_minus", "APOLO_G_plus",
                                "PID_D_minus", "PID_D_plus"
         )) %>%
  ggbarplot(., x = "percent_input", y = "value",add = "mean_se", fill="lightgrey", facet.by = "position_RNAseA", scales="free",ncol=2) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_binding_notaffected.pdf", plot=my_graph, width = 5, height = 15)


# Negative bind (affected by RNAsA)
my_graph <- 
  IP_all %>%
  filter(replicate != 1,
         position_RNAseA %in% c("FIL_b_minus", "FIL_b_plus",
                                "IAMT1_6_minus", "IAMT1_6_plus"
         )) %>%
  ggbarplot(., x = "percent_input", y = "value",add = "mean_se", fill="lightgrey", facet.by = "position_RNAseA", scales="free",ncol=2) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")
my_graph

# Save 
ggsave(filename="ChIP/out/SWI3B_rChIP_nobinding_affected.pdf", plot=my_graph, width = 5, height = 15)



# ChIP rice H3K27me3 exp1-------
Rice_H3K27me3_qPCR_pl1 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_pl1.xlsx", sheet=2) %>%
  select(-Pos)
Rice_H3K27me3_qPCR_pl2 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_pl2.xlsx", sheet=2) %>%
  select(-Pos)
Rice_H3K27me3_qPCR_pl3 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_exp12.xlsx", sheet=1) %>%
  filter(exp=="exp1")%>%
  select(-Pos, -exp) 
  

Rice_H3K27me3_qPCR <- Rice_H3K27me3_qPCR_pl1 %>%
  bind_rows(Rice_H3K27me3_qPCR_pl2) %>%
  bind_rows(Rice_H3K27me3_qPCR_pl3)


#Data processing
ChIP_qPCR_input <- Rice_H3K27me3_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

ChIP_qPCR_IP_input <- Rice_H3K27me3_qPCR %>% 
  filter(condition != "input") %>%
  left_join(ChIP_qPCR_input)


IP_input_H3 <- ChIP_qPCR_IP_input %>% 
  filter(condition == "H3") %>%
  mutate(Percent_input_H3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3 <- ChIP_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_H3 %>%
  left_join(IP_input_H3K27me3) 

Percent_input_all %>%
  mutate(class = fct_reorder(position, Percent_input_H3K27me3, .fun='median')) %>%
  ggplot(aes(x=reorder(position, Percent_input_H3K27me3), y=Percent_input_H3K27me3)) + 
  geom_boxplot()+
  geom_abline(slope=0, intercept=0.35,  col = "red",lty=2)

Percent_input_all %>%
  mutate(class = fct_reorder(position, Percent_input_H3K27me3, .fun='median')) %>%
  ggplot(aes(x=reorder(position, Percent_input_H3K27me3), y=Percent_input_H3K27me3)) + 
  geom_boxplot()+
  geom_abline(slope=0, intercept=0.35,  col = "red",lty=2)

Percent_input_all %>%
  filter(position %in% c("Actin1", "Pos1", "Pos2"
  )) %>%
  ggbarplot(., x = "position", y = "Percent_input_H3",add = "mean_se", fill="lightgrey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")

#LAB MEETING ----
Percent_input_all %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1", "03g01270-Pwz","MADS6","Actin1","Pos1","Pos2"
  )) %>%
  ggbarplot(., x = "position", y = "Percent_input_H3K27me3",add = "mean_se", fill="position") +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
Percent_input_all %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1","03g01270-Pwz","MADS6","Actin1","Pos1","Pos2"
  )) %>%
  ggbarplot(., x = "position", y = "Percent_input_H3",add = "mean_se", fill="position")+
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
#


IP_H3K27me3 <- Percent_input_all %>%
  mutate(Enrichment=Percent_input_H3K27me3/Percent_input_H3) %>%
  dplyr::select(replicate,position,Enrichment)


IP_H3K27me3 %>%
  ggplot(., aes(position, Enrichment))+
  geom_boxplot()


IP_H3K27me3 %>%
  mutate(class = fct_reorder(position, Enrichment, .fun='median')) %>%
  ggplot(aes(x=reorder(position, Enrichment), y=Enrichment)) + 
  geom_boxplot()+
  geom_abline(slope=0, intercept=0.0156,  col = "red",lty=2)


IP_H3K27me3 %>%
  filter(position %in% c("Actin", "MADS6")) %>%
  ggbarplot(., x = "position", y = "Enrichment",add = "mean_se", fill="lightgrey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")

IP_H3K27me3 %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1","03g01270-Pwz","MADS6","Actin1","Pos1","Pos2")) %>%
  ggbarplot(., x = "position", y = "Enrichment",add = "mean_se", fill="position") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))


IP_H3K27me3 %>%
  ggbarplot(., x = "position", y = "Enrichment",add = "mean_se", fill="lightgrey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")


# ChIP test AB GFP Arabido ------


#Data import
GFP_testAB_qPCR <- read_excel("ChIP/in/GFP_testAB_pl1.xlsx") %>%
  select(-Pos)


#Data processing
GFP_testAB_qPCR_input <- GFP_testAB_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(25)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- GFP_testAB_qPCR %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_IGG <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(IGG, names_to = "Antibody", values_to = "% input")
  
IP_input_CST <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("CST")) %>%
  mutate(CST = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(CST, names_to = "Antibody", values_to = "% input")
IP_input_Abcam <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("Abcam")) %>%
  mutate(Abcam = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(Abcam, names_to = "Antibody", values_to = "% input")
IP_input_Chromtek <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("Chromtek")) %>%
  mutate(Chromtek = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(Chromtek, names_to = "Antibody", values_to = "% input")

IP_all <- IP_input_IGG %>%
  bind_rows(IP_input_CST) %>%
  bind_rows(IP_input_Abcam) %>%
  bind_rows(IP_input_Chromtek)


my_graph <- 
  IP_all %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge()) +
  theme_bw() 
my_graph

# Save 
ggsave(filename="ChIP/out/arabido percent input_test AB all.pdf", plot=my_graph, width = 7, height = 4) 



IP_all$gene <- factor(IP_all$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))
IP_all$Antibody <- factor(IP_all$Antibody, levels=c("IGG","CST","Abcam"))


my_graph <- 
  IP_all %>%
  filter(gene %in% c("MRN1_promoter", "FLC-2", "housekeeping", "actin2", "FT", "YUCCA7"),
         Antibody %in% c("IGG","CST","Abcam")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "red") +
  theme_bw() +
  theme(legend.position = c(0.87, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))
my_graph

# Save 
ggsave(filename="ChIP/out/arabido percent input_test AB enriched and not enriched.pdf", plot=my_graph, width = 10, height = 4) 


# Stat

IP_all$gene <- factor(IP_all$gene, levels=c( "actin2", "housekeeping", "YUCCA7", "FT",  "MRN1_promoter", "FLC-2"))

my_graph <- 
IP_all %>%
  filter(gene %in% c("MRN1_promoter", "FLC-2", "housekeeping", "actin2", "FT", "YUCCA7"),
         Antibody %in% c("Abcam")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="#619CFF", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "red") +
stat_compare_means(method="t.test", 
                   comparisons = list(c("actin2", "housekeeping"),
                                      c("actin2", "YUCCA7"),
                                      c("actin2", "FT"),
                                      c("actin2", "MRN1_promoter"),
                                      c("actin2", "FLC-2")), 
                   label = "p.signif", bracket.size= 0.5)+
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph 

# Save 

ggsave(filename="ChIP/out/GFP Abcam signif.pdf", plot=my_graph, width = 3, height = 4) 






# ChIP rice H3K27me3 exp2-------

Rice_H3K27me3_qPCR_exp2_pl1 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_exp2_pl1.xlsx", sheet=2) %>%
  select(-Pos)
Rice_H3K27me3_qPCR_exp2_pl2 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_exp2_pl2.xlsx", sheet=2) %>%
  select(-Pos)
Rice_H3K27me3_qPCR_pl3 <- read_excel("ChIP/in/Rice_H3K27me3_qPCR_exp12.xlsx", sheet=1) %>%
  filter(exp=="exp2")%>%
  select(-Pos, -exp) 

Rice_H3K27me3_qPCR <- Rice_H3K27me3_qPCR_exp2_pl1 %>%
  bind_rows(Rice_H3K27me3_qPCR_exp2_pl2)%>%
  bind_rows(Rice_H3K27me3_qPCR_pl3)%>%
  mutate(Cp=as.numeric(Cp))


#Data processing
ChIP_qPCR_input <- Rice_H3K27me3_qPCR %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

ChIP_qPCR_IP_input <- Rice_H3K27me3_qPCR %>% 
  filter(condition != "input") %>%
  left_join(ChIP_qPCR_input)


IP_input_H3 <- ChIP_qPCR_IP_input %>% 
  filter(condition == "H3") %>%
  mutate(Percent_input_H3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_input_H3K27me3 <- ChIP_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(Percent_input_H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

Percent_input_all <- IP_input_H3 %>%
  left_join(IP_input_H3K27me3) 



Percent_input_all %>%
  mutate(class = fct_reorder(position, Percent_input_H3K27me3, .fun='median')) %>%
  ggplot(aes(x=reorder(position, Percent_input_H3K27me3), y=Percent_input_H3K27me3)) + 
  geom_boxplot()+
  ylim(0,2)+
  geom_abline(slope=0, intercept=0.35,  col = "red",lty=2)

#LAB MEETING ----
Percent_input_all %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1", "03g01270-Pwz","MADS6","Actin1","Pos1","Pos2"
  )) %>%
  ggbarplot(., x = "position", y = "Percent_input_H3K27me3",add = "mean_se", fill="position") +
  theme_bw() +
  theme(legend.position="none") +
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
Percent_input_all %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1", "03g01270-Pwz","MADS6","Actin1","Pos1","Pos2"
  )) %>%
  ggbarplot(., x = "position", y = "Percent_input_H3",add = "mean_se", fill="position") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("IGG", "GFP")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
#





IP_H3K27me3 <- Percent_input_all %>%
  mutate(Enrichment=Percent_input_H3K27me3/Percent_input_H3) %>%
  dplyr::select(replicate,position,Enrichment)


IP_H3K27me3 %>%
  ggplot(., aes(position, Enrichment))+
  geom_boxplot()


IP_H3K27me3 %>%
  mutate(class = fct_reorder(position, Enrichment, .fun='median')) %>%
  ggplot(aes(x=reorder(position, Enrichment), y=Enrichment)) + 
  geom_boxplot()+
  geom_abline(slope=0, intercept=0.0156,  col = "red",lty=2)



IP_H3K27me3 %>%
  filter(position %in% c("Actin", "MADS87", 
                         "MADS22", "MADS34", "MADS1", "03g01270-Pwz","MADS6","Actin1","Pos1","Pos2")) %>%
  ggbarplot(., x = "position", y = "Enrichment",add = "mean_se", fill="position")  +
  theme_bw() +
  theme(legend.position="none")+
  xlab("")+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))




# ChIP test AB GFP Arabido ------


#Data import

GFP_DIRINDIR_pl1 <- read_excel("ChIP/in/GFP_DIRINDIR_pl1.xlsx",sheet=2)%>%
  select(-Pos)


#Data processing
GFP_testAB_qPCR_input <- GFP_DIRINDIR_pl1 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- GFP_DIRINDIR_pl1 %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_IGG_direct <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG_direct") %>%
  mutate(IGG_direct = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(IGG_direct, names_to = "Antibody", values_to = "% input")
IP_input_IGG_indirect <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG_indirect") %>%
  mutate(IGG_indirect = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(IGG_indirect, names_to = "Antibody", values_to = "% input")
IP_input_GFP_direct <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "GFP_direct") %>%
  mutate(GFP_direct = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(GFP_direct, names_to = "Antibody", values_to = "% input")
IP_input_GFP_indirect <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "GFP_indirect") %>%
  mutate(GFP_indirect = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(GFP_indirect, names_to = "Antibody", values_to = "% input")

IP_all <- IP_input_IGG_direct %>%
  bind_rows(IP_input_IGG_indirect) %>%
  bind_rows(IP_input_GFP_direct) %>%
  bind_rows(IP_input_GFP_indirect)


my_graph <- 
  IP_all %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge()) +
  theme_bw() 
my_graph

# Save 
ggsave(filename="ChIP/out/LHP1 direct indirect.pdf", plot=my_graph, width = 7, height = 4) 



IP_all$gene <- factor(IP_all$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all %>%
  filter(gene %in% c("MRN1_promoter", "FLC-2", "housekeeping", "actin2", "FT", "YUCCA7"),
         Antibody %in% c("IGG_direct","GFP_direct")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Direct ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Direct_LHP1.pdf", plot=my_graph, width = 4, height = 4) 


my_graph <- 
  IP_all %>%
  filter(gene %in% c("MRN1_promoter", "FLC-2", "housekeeping", "actin2", "FT", "YUCCA7"),
         Antibody %in% c("IGG_indirect","GFP_indirect")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Indirect ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Indirect_LHP1.pdf", plot=my_graph, width = 4, height = 4) 


# Fold enrichment ------



IP_input_IGG_direct <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG_direct") %>%
  mutate(IGG_direct = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_IGG_indirect <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG_indirect") %>%
  mutate(IGG_indirect = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_GFP_direct <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "GFP_direct") %>%
  mutate(GFP_direct = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
IP_input_GFP_indirect <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "GFP_indirect") %>%
  mutate(GFP_indirect = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

IP_all_Enrichment_direct <- 
  IP_input_IGG_direct %>%
  left_join(IP_input_IGG_direct) %>%
  left_join(IP_input_GFP_direct) %>%
  mutate(Enrichment=GFP_direct/IGG_direct)


IP_all_Enrichment_direct$gene <- factor(IP_all_Enrichment_direct$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_Enrichment_direct %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="grey", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Direct ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Direct_LHP1_Enrichment.pdf", plot=my_graph, width = 4, height = 4) 



IP_all_Enrichment_indirect <- 
  IP_input_IGG_direct %>%
  left_join(IP_input_IGG_indirect) %>%
  left_join(IP_input_GFP_indirect) %>%
  mutate(Enrichment=GFP_indirect/IGG_indirect)


IP_all_Enrichment_indirect$gene <- factor(IP_all_Enrichment_indirect$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_Enrichment_indirect %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="grey", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Indirect ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Indirect_LHP1_Enrichment.pdf", plot=my_graph, width = 4, height = 4) 






# Normalization to Actin2 ------

## DIRECT

#Cp mean ref
ref_gene <- "actin2"
ref_data <- IP_all %>%
  filter(gene == ref_gene,
         Antibody %in% c("IGG_direct","GFP_direct")) %>%
  rename(`% input ref`=`% input`) %>%
  dplyr::select(-gene)

gene_data <- IP_all %>%
  filter(gene != ref_gene,
         Antibody %in% c("IGG_direct","GFP_direct")) 

# 
IP_all_normalized <- gene_data %>%
  left_join(ref_data) %>%
  mutate(`% input normalized to Actin2`=`% input` - `% input ref`)


IP_all_normalized$gene <- factor(IP_all_normalized$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_normalized %>%
  ggbarplot(., x = "gene", y = "% input normalise to Actin2",add = "mean_se", fill="Antibody", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Direct ChIP normalized Actin 2") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Direct_LHP1_norm Actin2.pdf", plot=my_graph, width = 4, height = 4) 





## INDIRECT

#Cp mean ref
ref_gene <- "actin2"
ref_data <- IP_all %>%
  filter(gene == ref_gene,
         Antibody %in% c("IGG_indirect","GFP_indirect")) %>%
  rename(`% input ref`=`% input`) %>%
  dplyr::select(-gene)

gene_data <- IP_all %>%
  filter(gene != ref_gene,
         Antibody %in% c("IGG_indirect","GFP_indirect")) 

# 
IP_all_normalized <- gene_data %>%
  left_join(ref_data) %>%
  mutate(`% input normalized to Actin2`=`% input` - `% input ref`)


IP_all_normalized$gene <- factor(IP_all_normalized$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_normalized %>%
  ggbarplot(., x = "gene", y = "% input normalized to Actin2",add = "mean_se", fill="Antibody", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Indirect ChIP normalized Actin 2") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Indirect_LHP1_norm Actin2.pdf", plot=my_graph, width = 4, height = 4) 


# Fold enrichment normalized ----
## DIRECT

IP_all_Enrichment_direct <-
  IP_all_Enrichment_direct %>%
  dplyr::select(replicate, gene,Enrichment)

ref_data <- IP_all_Enrichment_direct %>%
  filter(gene == ref_gene) %>%
  rename(Enrichment_ref=Enrichment) %>%
  dplyr::select(-gene)

gene_data <- IP_all_Enrichment_direct %>%
  filter(gene != ref_gene) 

# 
IP_all_Enrichment_direct_normalized <- gene_data %>%
  left_join(ref_data) %>%
  mutate(`Enrichment normalized to Actin2`=Enrichment - Enrichment_ref)




IP_all_Enrichment_direct_normalized$gene <- factor(IP_all_Enrichment_direct_normalized$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_Enrichment_direct_normalized %>%
  ggbarplot(., x = "gene", y = "Enrichment normalized to Actin2",add = "mean_se", fill="grey", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Direct ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Direct_LHP1_Enrichment normalized.pdf", plot=my_graph, width = 4, height = 4) 


## INDIRECT

IP_all_Enrichment_indirect <-
  IP_all_Enrichment_indirect %>%
  dplyr::select(replicate, gene,Enrichment)

ref_data <- IP_all_Enrichment_indirect %>%
  filter(gene == ref_gene) %>%
  rename(Enrichment_ref=Enrichment) %>%
  dplyr::select(-gene)

gene_data <- IP_all_Enrichment_indirect %>%
  filter(gene != ref_gene) 

# 
IP_all_Enrichment_indirect_normalized <- gene_data %>%
  left_join(ref_data) %>%
  mutate(`Enrichment normalized to Actin2`=Enrichment - Enrichment_ref)




IP_all_Enrichment_indirect_normalized$gene <- factor(IP_all_Enrichment_indirect_normalized$gene, levels=c("FLC-2", "MRN1_promoter", "FT", "YUCCA7", "actin2", "housekeeping"))


my_graph <- 
  IP_all_Enrichment_indirect_normalized %>%
  ggbarplot(., x = "gene", y = "Enrichment normalized to Actin2",add = "mean_se", fill="grey", position=position_dodge(.8),
            label = TRUE, label.pos = "out", lab.nb.digits=2, lab.vjust=-1, lab.col = "blue") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.77),
        legend.background = element_rect(fill = "white", color = "black"))+
  ggtitle("Indirect ChIP") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1))
my_graph

# Save 
ggsave(filename="ChIP/out/Indirect_LHP1_Enrichment normalized.pdf", plot=my_graph, width = 4, height = 4) 



# Test H3K27me3 antibody in Arabidopsis -------

H3K27me3_testAB_pl1 <- read_excel("ChIP/in/H3K27me3_testAB_pl1.xlsx") %>%
  select(-Pos)



GFP_testAB_qPCR_input <- H3K27me3_testAB_pl1 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- H3K27me3_testAB_pl1 %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_H3 <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "H3") %>%
  mutate(H3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3, names_to = "Antibody", values_to = "% input")

IP_input_H3K27me3 <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(H3K27me3_invitrogen = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3K27me3_invitrogen, names_to = "Antibody", values_to = "% input")

IP_input_H3K27me3_good_tidy <- IP_input_H3K27me3_good %>%
  rename("% input"=Percent_input_H3K27me3,
         gene=region) %>%
  add_column(Antibody="H3K27me3_millipore")

IP_all <- IP_input_H3 %>%
  bind_rows(IP_input_H3K27me3) %>%
  bind_rows(IP_input_H3K27me3_good_tidy)

my_graph <- 
  IP_all %>% 
  filter(gene%in% c("MRN1", "housekeeping", "FLC-2"))%>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge()) +
  facet_wrap(~Antibody, scale="free")+
  theme_bw() 
my_graph


# Save 
ggsave(filename="ChIP/out/Invitrogen vs millipore H3K27me3.pdf", plot=my_graph, width = 7, height = 4) 



# ChIP SWI3B exp2 RNAse A -------
## Only Actin2, APOLO_G, FIL_a are ok for exp2 and Actin2 only for exp1


exp2SWI3BGFP_rChIP_qPCR_pl1 <- read_excel("ChIP/in/exp2SWI3BGFP_rChIP_qPCR_pl1.xlsx") %>%
  filter(gene =="actin2")%>%
  select(-Pos)
exp2SWI3BGFP_rChIP_qPCR_pl2 <- read_excel("ChIP/in/exp2SWI3BGFP_rChIP_qPCR_pl2.xlsx") %>%
  filter(gene %in% c("actin2", "APOLO_G"))%>%
  select(-Pos)





GFP_testAB_qPCR_input <- exp2SWI3BGFP_rChIP_qPCR_pl2 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- exp2SWI3BGFP_rChIP_qPCR_pl2 %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_IGG <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IGG") %>%
  mutate(IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(IGG, names_to = "Antibody", values_to = "% input")

IP_input_GFP <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("GFP")) %>%
  mutate(GFP = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(GFP, names_to = "Antibody", values_to = "% input")


IP_all <- IP_input_IGG %>%
  bind_rows(IP_input_GFP) 


stat_IP_all <- IP_all %>%
  group_by(RNAse_A, gene,Antibody) %>%
  summarise(mean=mean(`% input`), 
            median= median(`% input`),
            SD=sd(`% input`), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 





my_graph <- 
ggplot(stat_IP_all, aes(gene, mean, fill=Antibody)) +
  geom_col(position = 'dodge', colour = 'black') +
  geom_errorbar(aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std)), 
                size=.5, width=.2, colour = 'black',
                position = position_dodge(0.9))+
  facet_wrap(~RNAse_A, nrow=1)+
  xlab(label = "")+
  ylab(label = "% input")
my_graph


# ChIP WT as control -------

WTcontrol_pl1 <- read_excel("ChIP/in/WTcontrol_pl1.xlsx") %>%
  select(-Pos)







GFP_WT_input <- WTcontrol_pl1 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_WT_IP_input <- WTcontrol_pl1 %>% 
  filter(condition != "input") %>%
  left_join(GFP_WT_input)


IP_input_WT <- GFP_WT_IP_input %>% 
  filter(genotype == "WT") %>%
  mutate("% input" = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)

IP_input_GFP <- GFP_WT_IP_input %>% 
  filter(genotype == "GFP") %>%
  mutate("% input" = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


IP_all <- IP_input_WT %>%
  bind_rows(IP_input_GFP) 


stat_IP_all <- IP_all %>%
  group_by(genotype, gene) %>%
  summarise(mean=mean(`% input`), 
            median= median(`% input`),
            SD=sd(`% input`), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


stat_IP_all$gene <- factor(stat_IP_all$gene, levels=c("FLC-2", "MRN1", "FT", "actin2", "ref4", "YUCCA7"))



my_graph <- 
  stat_IP_all %>% filter(gene %in% c("actin2", "FT", "MRN1", "ref4", "YUCCA7")) %>%
  ggplot(., aes(gene, mean, fill=genotype)) +
  geom_col(position = 'dodge', colour = 'black') +
  geom_errorbar(aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std)), 
                size=.5, width=.2, colour = 'black',
                position = position_dodge(0.9))+
  xlab(label = "")+
  ylab(label = "% input")
my_graph

# Test Millipore H3K27me3 AB Indirect metehod Arabidopsps badly cross link ----

Arabidopsis_Millipore_H3K27me3_pl1 <- read_excel("ChIP/in/Arabidopsis_Millipore_H3K27me3_pl1.xls") %>%
  select(-Pos)



GFP_testAB_qPCR_input <- Arabidopsis_Millipore_H3K27me3_pl1 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- Arabidopsis_Millipore_H3K27me3_pl1 %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_H3 <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "H3") %>%
  mutate(H3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3, names_to = "Antibody", values_to = "% input")

IP_input_H3K27me3 <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3")) %>%
  mutate(H3K27me3 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3K27me3, names_to = "Antibody", values_to = "% input")


IP_all <- IP_input_H3 %>%
  bind_rows(IP_input_H3K27me3)



IP_all$gene <- factor(IP_all$gene, levels=c("MRN1", "FT", "Actin2", "ActinCourtney", "PP2A", "YUCCA7"))


my_graph <- 
  IP_all %>%
  filter(gene%in% c("MRN1", "FT", "Actin2", "ActinCourtney")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge()) +
  facet_wrap(~Antibody, scale="free")+
  theme_bw() 
my_graph


# Test H3K27me3 AB direct method Arabidopsps well cross link (Glass Bell) ----

Arabidopsis_H3K27me3_pl1 <- read_excel("2022_PostDoc_PreliminaryWorks/R/ChIP/in/Arabidopsis_Millipore_Thermo_GlassBell_H3K27me3_pl1.xls") %>%
  select(-Pos)



GFP_testAB_qPCR_input <- Arabidopsis_H3K27me3_pl1 %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

GFP_testAB_qPCR_IP_input <- Arabidopsis_H3K27me3_pl1 %>% 
  filter(condition != "input") %>%
  left_join(GFP_testAB_qPCR_input)


IP_input_IGG <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition == "IgG") %>%
  mutate(IgG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(IgG, names_to = "Antibody", values_to = "% input")

IP_input_H3K27me3MP <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3_Millipore")) %>%
  mutate(H3K27me3MP = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3K27me3MP, names_to = "Antibody", values_to = "% input")
IP_input_H3K27me3Thermo <- GFP_testAB_qPCR_IP_input %>% 
  filter(condition %in% c("H3K27me3_Thermo")) %>%
  mutate(H3K27me3Thermo = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)%>% 
  pivot_longer(H3K27me3Thermo, names_to = "Antibody", values_to = "% input")


IP_all <- IP_input_IGG %>%
  bind_rows(IP_input_H3K27me3MP) %>%
  bind_rows(IP_input_H3K27me3Thermo)



IP_all$gene <- factor(IP_all$gene, levels=c("Actin2", "ActinCourtney", "MRN1", "FT"))


my_graph <- 
  IP_all %>%
  filter(gene%in% c("Actin2", "ActinCourtney", "MRN1", "FT")) %>%
  ggbarplot(., x = "gene", y = "% input",add = "mean_se", fill="Antibody", position=position_dodge(.7)) +
  theme_bw() 
my_graph

