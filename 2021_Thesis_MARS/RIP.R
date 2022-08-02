# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/RIP/", showWarnings = FALSE, recursive = TRUE)

#Data import
RIP_qPCR_output <- read_excel("data/RIP/RIP_qPCR_output.xlsx",sheet=1)

#Data processing
RIP_percent <- RIP_qPCR_output %>%
  mutate(Cp_input=Cp_input-log2(10)) %>%
  mutate(Percent=2^-(Cp_sample-Cp_input)*100) 



RIP_percent_LHP1 <- RIP_percent %>%
  filter(condition =="GFP-LHP1") %>%
  mutate(Percent_LHP1=Percent) %>%
  dplyr::select(replicate, gene, Percent_LHP1)
RIP_percent_IGG <- RIP_percent %>%
  filter(condition =="Igg") %>%
  mutate(Percent_Igg=Percent)%>%
  dplyr::select(replicate, gene, Percent_Igg)

RIP_Percent_Ratio <- RIP_percent_LHP1 %>%
  left_join(RIP_percent_IGG) %>%
  mutate(ratio=Percent_LHP1/Percent_Igg)



# t.test


RIP_Percent_Ratio$gene <- factor(RIP_Percent_Ratio$gene, levels=c("PP2A", "MRN1", "APOLO","MARS"))

my_graph <- 
  RIP_Percent_Ratio %>%
  ggbarplot(., x = "gene", y = "ratio",add = "mean_se", fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("PP2A","MRN1"), c("PP2A", "APOLO"), c("PP2A", "MARS") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  scale_fill_manual(values=selected_gene_RIP$gene_color,) + 
  theme(legend.position="none")

my_graph



# Save 
ggsave(filename="out/RIP/RIP.pdf", plot=my_graph, width = 6, height = 4) 



# RIP repet isof-------


rip_repet_1 <- read_excel("data/RIP/rip_repet_1.xlsx") %>%
  mutate(Cp=as.numeric(Cp))
rip_repet_2 <- read_excel("data/RIP/rip_repet_2.xlsx")
rip_repet_3 <- read_excel("data/RIP/rip_repet_3.xlsx")
rip_repet_4 <- read_excel("data/RIP/rip_repet_4.xlsx")




rip_repet <- rip_repet_1 %>%
  bind_rows(rip_repet_2) %>%
  bind_rows(rip_repet_3) %>%
  bind_rows(rip_repet_4) %>%
  select(-Pos)%>%
  mutate(replicate=as.character(replicate))%>%
  drop_na()


rip_repet <- read_excel("data/RIP/rip_repet_all.xlsx")




rip <- rip_repet %>%
  group_by(condition,replicate,gene)%>%
  summarise(Cp=mean(Cp)) %>%
  ungroup()





rip_input <- rip %>% 
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10)) %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)


# Value input in vivo RIP for Referee MP -----------
rip_input %>%
  filter(replicate %in% c(1,2,4)) %>%
  group_by(gene) %>%
  summarise(mean=mean(Cp_input))
  
  
#

rip_input_IGG <- rip %>% 
  left_join(rip_input)%>%
  filter(condition == "IGG") %>%
  mutate(Percent_input_IGG = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
rip_input_GFP <-  rip %>% 
  left_join(rip_input)%>%
  filter(condition == "GFP") %>%
  mutate(Percent_input_GFP = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


rip_IGG_GFP <- rip_input_IGG %>%
  left_join(rip_input_GFP) %>%
  mutate(ratio=Percent_input_GFP/Percent_input_IGG)




rip_IGG_GFP$gene <- factor(rip_IGG_GFP$gene, levels=c("PP2A","ref2", "ref3","ref4", "MRN1", "APOLO", "MARS.1","MARS.2","MARS.3","MARS"))

rip_IGG_GFP %>%
  filter(replicate %in% c(1,2,4)) %>%
  ggplot(., aes(gene,ratio)) +
  geom_boxplot()+
  geom_jitter(aes(color=replicate))
  

rip_IGG_GFP  %>%
  filter(gene %in% c("ref1","APOLO","MRN1","MARS.1","MARS.2","MARS.3","all2")) %>%
  ggplot(., aes(gene,ratio)) +
  geom_boxplot()+
  geom_jitter(aes(color=replicate))



rip_IGG_GFP %>%
  filter(gene %in% c("ref1","APOLO","MRN1","MARS.1","MARS.2","MARS.3","all2"),
    replicate %in% c(1,2,4)) %>%
  ggplot(., aes(gene,ratio)) +
  geom_boxplot()+
  geom_jitter(aes(color=replicate))





rip_IGG_GFP %>%
  filter(gene %in% c("PP2A","APOLO","MRN1","MARS.1","MARS.2","MARS.3","MARS")) %>%
  ggbarplot(., x = "gene", y = "ratio",add = c("mean_se","point"), fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("PP2A","MRN1"), c("PP2A", "APOLO"), c("PP2A", "MARS.1"), c("PP2A","MARS.2"),
                                        c("PP2A","MARS.3"), c("PP2A","MARS")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")



my_graph <- 
rip_IGG_GFP %>%
  filter(gene %in% c("ref2","APOLO","MRN1","MARS.1","MARS.2","MARS.3","MARS"),
         replicate %in% c(1,2,4)) %>%
  ggbarplot(., x = "gene", y = "ratio",add = c("mean_se"), fill="gene", width=0.95) +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref2","MRN1"), c("ref2", "APOLO"), c("ref2", "MARS.1"), c("ref2","MARS.2"),
                                        c("ref2","MARS.3"), c("ref2","MARS")), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("Enrichment") +
  xlab("") +
  theme(legend.position="none")
my_graph



##ref2 is the PP2A !!

# Save 
ggsave(filename="out/RIP/RIP_isoform.pdf", plot=my_graph, width = 2, height = 3) 






# in-vitro ------


read_excel("data/RIP/qPCR_output_invitroRIP.xlsx",sheet=2)%>%
  mutate(Cp=as.numeric(Cp))%>%
  select(-Position)%>%
  group_by(gene,replicate,condition)%>%
  summarise(Cp=mean(Cp))%>%
  drop_na() %>%
  filter(condition!="input",gene!="MARS")%>%
  ggplot(., aes(gene,Cp))+
  geom_point(aes(color=replicate),position="dodge")+
  facet_grid(~condition)





#R3 LHP1 and R1 NSR is shit
qPCR_output_invitroRIP <- read_excel("data/RIP/qPCR_output_invitroRIP.xlsx",sheet=1)%>%
  mutate(Cp=as.numeric(Cp))%>%
  select(-Position)%>%
  group_by(gene,replicate,condition)%>%
  summarise(Cp=mean(Cp))%>%
  drop_na()

write.csv(qPCR_output_invitroRIP, file = "qPCR_output_invitroRIP.csv")

rip_input <- qPCR_output_invitroRIP %>% 
  filter(condition == "input") %>%
  group_by(gene)%>%
  summarise(Cp=mean(Cp))%>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp)

rip_input_LHP1 <- qPCR_output_invitroRIP %>% 
  left_join(rip_input)%>%
  filter(condition == "LHP1") %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
rip_input_NSR <-  qPCR_output_invitroRIP %>% 
  left_join(rip_input)%>%
  filter(condition == "NSR") %>%
  mutate(Percent_input_NSR = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


rip_LHP1 <- rip_input_LHP1 %>%
  pivot_longer(Percent_input_LHP1, names_to = "condition")
rip_NSR <- rip_input_NSR %>%
  pivot_longer(Percent_input_NSR, names_to = "condition")

rip_LHP1_NSR <- rip_LHP1 %>%
  bind_rows(rip_NSR)
  
rip_LHP1_NSR$gene <- factor(rip_LHP1_NSR$gene, levels=c("ASCO","APOLO","MARS_iso1","MARS_iso2"))

 
 
 rip_LHP1_NSR %>%
   filter(gene !="MARS",condition=="Percent_input_LHP1") %>%
   ggbarplot(., x = "gene", y = "value",add = c("mean_se"), fill="gene")+
   stat_compare_means(method="t.test", 
                      comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                      label = "p.label", bracket.size= 0.5) +
   theme_bw() +
   ylab("% input") +
   xlab("") +
   theme(legend.position="none")+
   ggtitle("RIP LHP1")
 
 


rip_LHP1_NSR %>%
  filter(gene !="MARS",condition=="Percent_input_NSR") %>%
  ggbarplot(., x = "gene", y = "value",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP NSR")

# as ratio

rip_LHP1_NSR_LHP1_input <- rip_LHP1_NSR %>%
  filter(condition == "Percent_input_LHP1")%>%
  ungroup()%>%
  select(-condition,-replicate)%>%
  rename(LHP1_input=value)

rip_LHP1_NSR_NSR_input <- rip_LHP1_NSR %>%
  filter(condition == "Percent_input_NSR")%>%
  ungroup()%>%
  select(-condition,-replicate)%>%
  rename(NSR_input=value)


rip_LHP1_NSR_ratio_input <- rip_LHP1_NSR_LHP1_input %>%
  left_join(rip_LHP1_NSR_NSR_input)%>%
  mutate(ratio=LHP1_input/NSR_input)


rip_LHP1_NSR_ratio_input %>%
  filter(gene !="MARS") %>%
  ggbarplot(., x = "gene", y = "ratio",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("ratio") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1 / RIP NSR")



#in vitro repet with dilution and noise control -------------

qPCR_output_invitro3 <- read_excel("data/RIP/qPCR_output_invitro3.xlsx",sheet=2) %>%
  select(-Position)%>%
  mutate(Cp=as.numeric(Cp))%>%
  drop_na()

# percent input LHP and NSR

qPCR_output_invitro3_input_sample <- qPCR_output_invitro3 %>%
  filter(condition == "input_sample")%>%
  select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))



rip_input_LHP1 <- qPCR_output_invitro3 %>% 
  left_join(qPCR_output_invitro3_input_sample)%>%
  filter(condition == "LHP1") %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


rip_input_NSR <-  qPCR_output_invitro3 %>% 
  left_join(qPCR_output_invitro3_input_sample)%>%
  filter(condition == "NSR") %>%
  mutate(Percent_input_NSR = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)


rip_LHP1 <- rip_input_LHP1 %>%
  pivot_longer(Percent_input_LHP1, names_to = "condition")
rip_NSR <- rip_input_NSR %>%
  pivot_longer(Percent_input_NSR, names_to = "condition")

rip_LHP1_NSR <- rip_LHP1 %>%
  bind_rows(rip_NSR)





rip_LHP1_NSR$gene <- factor(rip_LHP1_NSR$gene, levels=c("ASCO","APOLO","MARS_iso1","MARS_iso2"))



rip_LHP1_NSR %>%
  filter(gene !="MARS",condition=="Percent_input_LHP1") %>%
  ggbarplot(., x = "gene", y = "value",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1")




rip_LHP1_GFP %>%
  filter(gene !="MARS",condition=="Percent_input_NSR") %>%
  ggbarplot(., x = "gene", y = "value",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP NSR")

# percent input noise


qPCR_output_invitro3_input_gfp <- qPCR_output_invitro3 %>%
  filter(condition == "input_gfp")%>%
  select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))

qPCR_output_invitro3_GFP <- qPCR_output_invitro3%>%
  filter(condition =="CTRL")%>%
  select(-replicate)%>%
  group_by(gene)%>%
  mutate(Cp=mean(Cp)) %>%
  unique()

rip_input_GFP <- qPCR_output_invitro3_GFP %>% 
  left_join(qPCR_output_invitro3_input_gfp)%>%
  filter(condition == "CTRL") %>%
  mutate(Percent_input_GFP = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)



rip_all <- rip_input_GFP %>%
  left_join(rip_input_LHP1) %>%
  left_join(rip_input_NSR)%>%
  mutate(Enrichment_LHP1=Percent_input_LHP1/Percent_input_GFP,
         Enrichment_NSR=Percent_input_NSR/Percent_input_GFP)




rip_all %>%
  ggbarplot(., x = "gene", y = "Enrichment_",add = c("mean_se"))+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1")



# Use the first qPCR ------------

## from qPCR 1:

rip_input_LHP1
rip_input_NSR


# add control from qPCR3:


rip_input_LHP1 %>%
  left_join(rip_input_GFP)%>%
  mutate(ratio_LHP1=Percent_input_LHP1/Percent_input_GFP)%>%
  ggbarplot(., x = "gene", y = "ratio_LHP1",add = c("mean_se"))

rip_input_NSR %>%
  left_join(rip_input_GFP)%>%
  mutate(ratio_LHP1=Percent_input_NSR/Percent_input_GFP)%>%
  ggbarplot(., x = "gene", y = "ratio_LHP1",add = c("mean_se"))




# plate 4 ----------


qPCR_output_invitro4 <- read_excel("data/RIP/qPCR_output_invitro4.xlsx",sheet=1)%>%
  mutate(Cp=as.numeric(Cp))%>%
  select(-Position)

qPCR_output_invitro4_ctrl <- qPCR_output_invitro4 %>%
  filter(condition == "input_sample")%>%
  select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))



rip_input_ctrl <- qPCR_output_invitro4 %>% 
  left_join(qPCR_output_invitro4_ctrl)%>%
  filter(condition == "GFP") %>%
  mutate(Percent_input_CTRL = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)



qPCR_output_invitro4_LHP1 <- qPCR_output_invitro4 %>%
  filter(condition == "input")%>%
  select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))




rip_input_LHP1 <- qPCR_output_invitro4 %>% 
  left_join(qPCR_output_invitro4_LHP1)%>%
  filter(condition == "LHP1") %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)













rip_LHP1 <- rip_input_LHP1 %>%
  pivot_longer(Percent_input_LHP1, names_to = "condition")
rip_NSR <- rip_input_NSR %>%
  pivot_longer(Percent_input_NSR, names_to = "condition")

rip_LHP1_NSR <- rip_LHP1 %>%
  bind_rows(rip_NSR)





rip_LHP1_NSR$gene <- factor(rip_LHP1_NSR$gene, levels=c("ASCO","APOLO","MARS_iso1","MARS_iso2"))



rip_LHP1_NSR %>%
  filter(gene !="MARS",condition=="Percent_input_LHP1") %>%
  ggbarplot(., x = "gene", y = "value",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1")


### RIP invitro all ------------
RIP_invitro_all <- read_excel("data/RIP/RIP_invitro_all.xlsx") %>%
  mutate(Cp=as.numeric(Cp))


RIP_invitro_all_ctrl <- RIP_invitro_all %>%
  filter(condition == "input_ctrl")%>%
  dplyr::select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))



RIP_invitro_all_input_ctrl <- RIP_invitro_all %>% 
  left_join(RIP_invitro_all_ctrl)%>%
  filter(condition == "CTRL") %>%
  mutate(Percent_input_CTRL = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)



RIP_invitro_all_LHP1 <- RIP_invitro_all %>%
  filter(condition == "input")%>%
  dplyr::select(-replicate,-condition) %>%
  group_by(gene) %>%
  summarise(Cp_input=mean(Cp))


RIP_invitro_all_input_LHP1 <- RIP_invitro_all %>% 
  left_join(RIP_invitro_all_LHP1)%>%
  filter(condition == "LHP1") %>%
  mutate(Percent_input_LHP1 = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)
RIP_invitro_all_input_NSR <- RIP_invitro_all %>% 
  left_join(RIP_invitro_all_LHP1)%>%
  filter(condition == "NSR") %>%
  mutate(Percent_input_NSR = 2^-(Cp-Cp_input)*100) %>%
  dplyr::select(-Cp, -Cp_input,-condition)



RIP_input_GFP_LHP1_NSR <- RIP_invitro_all_input_ctrl %>%
  left_join(RIP_invitro_all_input_LHP1) %>%
  left_join(RIP_invitro_all_input_NSR) %>%
  mutate(Enrichment_LHP1=Percent_input_LHP1/Percent_input_CTRL,
         Enrichment_NSR=Percent_input_NSR/Percent_input_CTRL)

RIP_input_GFP_LHP1_NSR$gene <- factor(RIP_input_GFP_LHP1_NSR$gene, levels=c("APOLO","ASCO","MARS_iso1","MARS_iso2","GFP"))



RIP_input_GFP_LHP1_NSR %>%
  ggbarplot(., x = "gene", y = "Percent_input_LHP1",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1")

RIP_input_GFP_LHP1_NSR %>%
  ggbarplot(., x = "gene", y = "Percent_input_NSR",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("% input") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP NSR")

RIP_input_GFP_LHP1_NSR %>%
  ggbarplot(., x = "gene", y = "Enrichment_LHP1",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("Enrichment") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP LHP1")

RIP_input_GFP_LHP1_NSR %>%
  ggbarplot(., x = "gene", y = "Enrichment_NSR",add = c("mean_se"), fill="gene")+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ASCO","APOLO"), c("ASCO", "MARS_iso1"), c("ASCO", "MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("Enrichment") +
  xlab("") +
  theme(legend.position="none")+
  ggtitle("RIP NSR")



RIP_input_GFP_tidy <- RIP_input_GFP_LHP1_NSR %>%
  pivot_longer(cols = "Percent_input_CTRL", names_to = "condition", values_to = "Percent_input")%>%
  select(gene,replicate,condition,Percent_input)
RIP_input_LHP1_tidy <- RIP_input_GFP_LHP1_NSR %>%
  pivot_longer(cols = "Percent_input_LHP1", names_to = "condition", values_to = "Percent_input")%>%
  select(gene,replicate,condition,Percent_input)
RIP_input_NSR_tidy <- RIP_input_GFP_LHP1_NSR %>%
  pivot_longer(cols = "Percent_input_NSR", names_to = "condition", values_to = "Percent_input")%>%
  select(gene,replicate,condition,Percent_input)

RIP_tidy <- RIP_input_GFP_tidy %>%
  bind_rows(RIP_input_LHP1_tidy)%>%
  bind_rows(RIP_input_NSR_tidy)

RIP_tidy_stat <- RIP_tidy%>%
  group_by(gene,condition)%>%
  summarise(mean=mean(Percent_input), 
            median= median(Percent_input),
            SD=sd(Percent_input), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

RIP_tidy_stat%>%
  filter(condition!="Percent_input_LHP1")%>%
  ggplot(., aes(gene,mean,fill=condition))+
  geom_col(position = 'dodge') +
  geom_errorbar(aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std)), position = position_dodge(0.9),width=.4)+
  ylab("% input")



# paper:

RIP_input_GFP_LHP1_NSR$gene <- factor(RIP_input_GFP_LHP1_NSR$gene, levels=c("GFP", "ASCO","APOLO", "MARS_iso1","MARS_iso2"))


my_graph <- 
RIP_input_GFP_LHP1_NSR %>%
  filter(gene%in%c("GFP", "ASCO","APOLO", "MARS_iso1", "MARS_iso2"))%>%
  ggbarplot(., x = "gene", y = "Enrichment_LHP1",add = c("mean_se"), fill="gene",width=0.95)+
  stat_compare_means(method="t.test", 
                     comparisons = list(c("GFP","ASCO"),c("GFP","APOLO"), c("GFP","MARS_iso1"), c("GFP", "MARS_iso2"), c("MARS_iso1","MARS_iso2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  ylab("Enrichment") +
  xlab("") +
  theme(legend.position="none")
my_graph


# Save 
ggsave(filename="out/RIP/RIP_invitro.pdf", plot=my_graph, width = 2, height = 3) 



