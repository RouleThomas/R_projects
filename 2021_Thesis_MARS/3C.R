# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/3C/", showWarnings = FALSE, recursive = TRUE)


#Data import
Col_RNAi <- read_excel("data/3C/Col_RNAi.xlsx") %>% clean_genotype_code
Col_ShortTime <- read_excel("data/3C/Col_ShortTime.xlsx") %>% clean_genotype_code

#Data processing_Col_RNAi
ref_gene <- "ctrl"
ref_data <- Col_RNAi %>%
  filter(position == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position)
gene_data <- Col_RNAi %>%
  filter(position != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>% #time 0 as control condition
  filter(genotype =="col",
         time == "0") %>%
  group_by(position) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_Col_RNAi <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt),
         time=as.character(time)) %>%
  dplyr::select(genotype,time,position,ddCt,RealTime)

#Data processing_Col_ShortTime
ref_gene <- "ctrl"
ref_data <- Col_ShortTime %>%
  filter(position == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position)
gene_data <- Col_ShortTime %>%
  filter(position != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>% #col and time 0 as control condition
  filter(genotype =="col",
         time == "0") %>%
  group_by(position) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_Col_ShortTime <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt),
         time=as.character(time)) %>%
  dplyr::select(genotype,time,position,ddCt,RealTime) 


#Combine Col_RNAi and Col_ShortTime
ddCt_data_3C <- ddCt_data_Col_RNAi %>%
  bind_rows(ddCt_data_Col_ShortTime)

# ANOVA-Col only ----
ddCt_data_3C_col <- ddCt_data_3C %>%
  filter(genotype =="col")
  
stat_ddCt_data_3C <- ddCt_data_3C_col %>%
  group_by(position, time, genotype, RealTime) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  
  
  

model_gene <- aov(ddCt~time,data=ddCt_data_3C_col)
anova(model_gene)
TukeyHSD(model_gene)

resume_lettre_gene <- multcompLetters(extract_p(TukeyHSD(model_gene)$`time`))
resume_lettre_df_gene <- data.frame(genotype=names(resume_lettre_gene$Letters),
                                    letter=resume_lettre_gene$Letters) %>%
  separate(col=genotype, into=c("time"), sep=":")
# Join letters to stat data
stat_all_gene <- stat_ddCt_data_3C %>% 
  left_join(resume_lettre_df_gene, by=c("time")) 
###

my_graph <- 
  stat_all_gene %>% left_join(genotype_metadata) %>%
  ggplot(.,aes(RealTime,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  geom_text(aes(label=letter), nudge_y = 0.2, vjust=1.5, hjust=-0.5, size=4)+
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=genotype_metadata$genotype,
                     values=genotype_metadata$genotype_color,
                     labels=genotype_metadata$genotype_name,)
my_graph

# Save 
ggsave(filename="out/3C/3C_col_kinetic.pdf", plot=my_graph, width = 5, height = 3)


# ANOVA-Col vs RNAi time0 ----
my_graph <- 
  ddCt_data_3C %>% filter(time =="0") %>%
  mutate(ddCt=-ddCt, time=as.numeric(time)) %>%
  left_join(selected_genotype_3C) %>%
  ggbarplot(., x = "genotype", y = "ddCt",add = "mean_se", fill="genotype") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("col","rnai38"),c("col", "rnai92")), 
                     label = "p.format", bracket.size= 0.5) +
  theme_bw() +
  xlab("") +
  ylab("log2(FC)") +
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,)
my_graph


# Save 
ggsave(filename="out/3C/3C_col vs RNAi t0.pdf", plot=my_graph, width = 5, height = 5)

# ANOVA-Col vs RNAi kinetic ----
ddCt_data_3C_col_RNAi_kin <- ddCt_data_3C %>%
  filter(RealTime %in% c(0,4))

stat_ddCt_data_3C_col_RNAi_kin <- ddCt_data_3C_col_RNAi_kin %>%
  group_by(position, time, genotype, RealTime) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  




model_gene <- aov(ddCt~time*genotype,data=ddCt_data_3C_col_RNAi_kin)
anova(model_gene)
TukeyHSD(model_gene)

resume_lettre_gene <- multcompLetters(extract_p(TukeyHSD(model_gene)$`time:genotype`))
resume_lettre_df_gene <- data.frame(genotype=names(resume_lettre_gene$Letters),
                                    letter=resume_lettre_gene$Letters) %>%
  separate(col=genotype, into=c("time","genotype"), sep=":")
# Join letters to stat data
stat_all_gene <- stat_ddCt_data_3C_col_RNAi_kin %>% 
  left_join(resume_lettre_df_gene, by=c("time","genotype")) %>%
  add_column(letter_corr= c("b","a","ab","","a","")) 


my_graph <- 
  stat_all_gene %>% left_join(selected_genotype_3C) %>%
  ggplot(.,aes(RealTime,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  geom_text(aes(label=letter_corr), nudge_y = 0.3, nudge_x = 0.1, size=4)+
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=genotype_metadata$genotype,
                     values=genotype_metadata$genotype_color,
                     labels=genotype_metadata$genotype_name,)
my_graph

# Save 
ggsave(filename="out/3C/3C_col vs RNAi kinetic.pdf", plot=my_graph, width = 5, height = 3)




# ANOVA-Col vs RNAi kinetic all points----
ddCt_data_3C_col_RNAi_kin <- ddCt_data_3C %>%
  filter(time != "8 hours", genotype != "rnai92") 

ddCt_data_3C_col_RNAi_kin_rnai92 <- ddCt_data_3C %>% filter(genotype =="rnai92", time != "2 hours", time != "8 hours")

ddCt_data_3C_col_RNAi_kin <- ddCt_data_3C_col_RNAi_kin_rnai92 %>% bind_rows(ddCt_data_3C_col_RNAi_kin)

stat_ddCt_data_3C_col_RNAi_kin <- ddCt_data_3C_col_RNAi_kin %>%
  group_by(position, time, genotype, RealTime) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  




model_gene <- aov(ddCt~time*genotype,data=ddCt_data_3C_col_RNAi_kin)
anova(model_gene)
TukeyHSD(model_gene)

resume_lettre_gene <- multcompLetters(extract_p(TukeyHSD(model_gene)$`time:genotype`))
resume_lettre_df_gene <- data.frame(genotype=names(resume_lettre_gene$Letters),
                                    letter=resume_lettre_gene$Letters) %>%
  separate(col=genotype, into=c("time","genotype"), sep=":")
# Join letters to stat data
stat_all_gene <- stat_ddCt_data_3C_col_RNAi_kin %>% 
  left_join(resume_lettre_df_gene, by=c("time","genotype")) %>%
  add_column(letter_corr= c("b","a","ab","","a","")) 



my_graph <- 
  stat_ddCt_data_3C_col_RNAi_kin %>% left_join(selected_genotype_3C) %>%
  ggplot(.,aes(RealTime,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=selected_genotype_3C$genotype,
                     values=selected_genotype_3C$genotype_color,
                     labels=selected_genotype_3C$genotype_name,)+ theme(legend.position = "none")
my_graph

# Save 
ggsave(filename="out/3C/3C_col vs RNAi kinetic_all.pdf", plot=my_graph, width = 3, height = 3)


# Col and lhp1 ABA -----

Col_lhp1 <- read_excel("data/3C/Col_lhp1.xlsx")

#Data processing_Col_RNAi
ref_gene <- "ctrl"
ref_data <- Col_lhp1 %>%
  filter(position == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position)
gene_data <- Col_lhp1 %>%
  filter(position != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>% #time 0 as control condition
  filter(genotype =="col",
         time == "0") %>%
  group_by(position) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_Col_lhp1 <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt),
         time=as.character(time)) %>%
  dplyr::select(genotype,time,position,ddCt)


stat_ddCt_data_3C_col_lhp1 <- ddCt_data_Col_lhp1 %>%
  mutate(time=as.character(time)) %>%
  group_by(position, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))
  


model_gene <- aov(ddCt~time*genotype,data=ddCt_data_Col_lhp1)
anova(model_gene)
TukeyHSD(model_gene)

resume_lettre_gene <- multcompLetters(extract_p(TukeyHSD(model_gene)$`time:genotype`))
resume_lettre_df_gene <- data.frame(genotype=names(resume_lettre_gene$Letters),
                                    letter=resume_lettre_gene$Letters) %>%
  separate(col=genotype, into=c("time","genotype"), sep=":")
# Join letters to stat data
stat_all_gene <- stat_ddCt_data_3C_col_lhp1 %>%
  left_join(resume_lettre_df_gene, by=c("time","genotype")) 



my_graph <- 
  stat_all_gene %>%
  ggplot(.,aes(time,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  geom_text(aes(label=letter), nudge_y = 0.3, nudge_x = 0.1, size=4)+
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=genotype_metadata$genotype,
                     values=genotype_metadata$genotype_color,
                     labels=genotype_metadata$genotype_name,)
my_graph

my_graph <- 
  stat_ddCt_data_3C_col_lhp1 %>%
ggplot(.,aes(time,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=genotype_metadata$genotype,
                     values=genotype_metadata$genotype_color,
                     labels=genotype_metadata$genotype_name,)+ theme(legend.position = "none")
my_graph

ggsave(filename="out/3C/3C_lhp1_kin.pdf", plot=my_graph, width = 3, height = 3)



my_graph <- 
  ggbarplot(ddCt_data_Col_lhp1, x = "time", y = "ddCt", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 2.5) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/3C/3C_lhp1.pdf", plot=my_graph, width = 3, height = 3)



# Combine RNAi and lhp1 in one plot ------

my_graph <- 
stat_ddCt_data_3C_col_RNAi_kin %>%
  ungroup()%>%
  select(position, genotype, RealTime, mean, median, SD, n ,erreur_std)%>%
  rename(time=RealTime) %>%
  bind_rows(stat_ddCt_data_3C_col_lhp1%>%
              mutate(time=as.numeric(time))) %>% left_join(selected_genotype_3C) %>%
  filter(genotype%in%c("col","lhp1","rnai38"))%>%
  ggplot(.,aes(time,mean, group=genotype))+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "time")+
  ylab(label = "log2(FC)") +
  scale_color_manual(breaks=selected_genotype_3C$genotype,
                     values=selected_genotype_3C$genotype_color,
                     labels=selected_genotype_3C$genotype_name,)+ theme(legend.position = "none")
my_graph

# Save 
ggsave(filename="out/3C/3C_rnai 1 and lhp1.pdf", plot=my_graph, width = 3, height = 3)




# loop detail ------

## exp A----
X3C_loop_detail <- read_excel("data/3C/loop_detail/3C_loop_detail.xlsx") %>%
  select(-Well,-orientation) %>%
  mutate(Cp=as.numeric(Cp)) %>%
  filter(exp=="C") 

ref_gene <- "control"
ref_data <- X3C_loop_detail %>%
  filter(position == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position,-position1)
gene_data <- X3C_loop_detail %>%
  filter(position != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(position1=="11")%>%
  group_by(exp)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_detail <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_3C_col_lhp1 <- ddCt_data_detail %>%
  group_by(position) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))


stat_ddCt_data_3C_col_lhp1 %>%
  left_join(read_excel("data/3C/loop_detail/3C_loop_detail.xlsx") %>%
              select(position,orientation,position1))%>%
  filter(position1 %in% c(2,3,4,5,6,7,8,9))%>%
  ggplot(.,aes(position1,mean,color=orientation), group=orientation)+
  geom_line(aes(color=orientation),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") 

## pull exp -----


X3C_loop_detail <- read_excel("data/3C/loop_detail/3C_loop_detail.xlsx")%>%
  drop_na(Cp)%>%
  select(-Well,-orientation) %>%
  mutate(Cp=as.numeric(Cp)) %>%
  group_by(exp,position,position1)%>%
  summarise(Cp=mean(Cp)) %>%
  ungroup()


ref_gene <- "control"
ref_data <- X3C_loop_detail %>%
  filter(position == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position,-position1)
gene_data <- X3C_loop_detail %>%
  filter(position != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  group_by(exp)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_detail <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_3C_col_lhp1 <- ddCt_data_detail %>%
  group_by(position) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
stat_ddCt_data_3C_col_lhp1 %>%
  left_join(read_excel("data/3C/loop_detail/3C_loop_detail.xlsx") %>%
              select(orientation,position1,position))%>%
  filter(position1 %in% c(1,2,3,4,5,6,7,8,9))%>%
  ggplot(.,aes(position1,mean,color=orientation), group=orientation)+
  geom_line(aes(color=orientation),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=2) +
  xlab(label = "")+
  ylab(label = "log2(FC)") 
my_graph


ggsave(filename="out/3C/3C_loopdetail_v2.pdf", plot=my_graph, width = 5, height = 3)











X3C_loop_detail <- read_excel("data/3C/loop_detail/3C_loop_detail.xlsx",sheet=2)%>%
  select(Cp,exp,replicate,orientation,position1)%>%
  drop_na(Cp)%>%
  mutate(Cp=as.numeric(Cp)) %>%
  group_by(exp,position1)%>%
  summarise(Cp=mean(Cp)) %>%
  ungroup()


ref_gene <- "control"
ref_data <- X3C_loop_detail %>%
  filter(position1 == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-position1)
gene_data <- X3C_loop_detail %>%
  filter(position1 != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(position1=="11")%>%
  group_by(exp)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_detail <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_3C_col_lhp1 <- ddCt_data_detail %>%
  group_by(position1) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
  stat_ddCt_data_3C_col_lhp1 %>%
  left_join(read_excel("data/3C/loop_detail/3C_loop_detail.xlsx") %>%
              select(orientation,position1,position))%>%
  filter(position1 %in% c(1,2,3,4,5,6,7,8,9,10))%>%
  ggplot(.,aes(position1,mean,color=orientation), group=orientation)+
  geom_line(aes(color=orientation),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") 
my_graph

## CL3 -----

#short------
X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx")%>%
  filter(exp=="short")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))
  

ref_gene <- "control"
ref_data <- X3C_loop1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop1 %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3_short <- ddCt_data_short %>%
  group_by(gene,time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
  stat_ddCt_data_CL3_short %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") 
my_graph



#lhp1------
X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx")%>%
  filter(exp=="lhp1")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


ref_gene <- "control"
ref_data <- X3C_loop1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop1 %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
  stat_ddCt_data_CL3 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph



# RNAI92 -----

X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx",sheet=2)%>%
  filter(exp=="92")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


ref_gene <- "control"
ref_data <- X3C_loop1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop1 %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3_92 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
  stat_ddCt_data_CL3_92 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph


# RNAI38 -----

X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx",sheet=2)%>%
  filter(exp=="38")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


ref_gene <- "control"
ref_data <- X3C_loop1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop1 %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3_38 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

my_graph <- 
  stat_ddCt_data_CL3_38 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph





# pull all exp ----------

stat_ddCt_data_CL3_38_pull <- stat_ddCt_data_CL3_38 %>%
  filter(genotype=="rnai38",
         gene=="CL3.5")
stat_ddCt_data_CL3_92_pull<- stat_ddCt_data_CL3_92 %>%
  filter(genotype=="rnai92",
         gene=="CL3.5")
stat_ddCt_data_CL3_short_col_pull <- stat_ddCt_data_CL3_short %>%
  filter(gene=="CL3.5")%>%
  add_column(genotype="col")

all_3C <- stat_ddCt_data_CL3_38_pull %>%
  bind_rows(stat_ddCt_data_CL3_92_pull)%>%
  bind_rows(stat_ddCt_data_CL3_short_col_pull)

my_graph <- 
  all_3C %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2)+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  scale_color_manual(breaks=selected_genotype_3C$genotype,
                     values=selected_genotype_3C$genotype_color,
                     labels=selected_genotype_3C$genotype_name,)+ theme(legend.position = "none")+
  ylim(-1.5,2.25)
my_graph



ggsave(filename="out/3C/CL_3.pdf", plot=my_graph, width = 5, height = 3)


# other CL1 valdiation ---------
##short -----

X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx")%>%
  filter(exp=="short")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


X3C_loop3 <- read_excel("data/3C/CL3/3C loop3.xlsx")%>%
  filter(exp=="short") %>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))

X3C_loop <- X3C_loop1 %>%
  bind_rows(X3C_loop3)


ref_gene <- "control"
ref_data <- X3C_loop %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3_short <- ddCt_data_short %>%
  group_by(gene,time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

stat_ddCt_data_short_1F8R_col <- stat_ddCt_data_CL3_short %>%
  filter(gene =="CL1_1F8R")

my_graph <- 
  stat_ddCt_data_CL3_short %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") 
my_graph


#lhp1 ----
# only smthing for 1F8R


X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx")%>%
  filter(exp=="lhp1")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


X3C_loop3 <- read_excel("data/3C/CL3/3C loop3.xlsx")%>%
  filter(exp=="lhp1") %>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))

X3C_loop <- X3C_loop1 %>%
  bind_rows(X3C_loop3)


ref_gene <- "control"
ref_data <- X3C_loop %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))



stat_ddCt_data_CL3 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

stat_ddCt_data_lhp1 <- stat_ddCt_data_CL3 %>%
  filter(genotype=="lhp1",gene=="CL1_1F8R")

my_graph <- 
  stat_ddCt_data_CL3 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph



# RNAi92------



X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx",sheet=2)%>%
  filter(exp=="92")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


X3C_loop3 <- read_excel("data/3C/CL3/3C loop3.xlsx",sheet=2)%>%
  filter(exp=="92")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))



X3C_loop <- X3C_loop1 %>%
  bind_rows(X3C_loop3)


ref_gene <- "control"
ref_data <- X3C_loop %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))

stat_ddCt_data_1F8R_92 <- stat_ddCt_data_CL3 %>%
  filter(genotype=="rnai92", gene =="CL1_1F8R")

my_graph <- 
  stat_ddCt_data_CL3 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph




# RNAi38------



X3C_loop1 <- read_excel("data/3C/CL3/3C loop1.xlsx",sheet=2)%>%
  filter(exp=="38")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))


X3C_loop3 <- read_excel("data/3C/CL3/3C loop3.xlsx",sheet=2)%>%
  filter(exp=="38")%>%
  select(-Pos,-exp)%>%
  mutate(Cp=as.numeric(Cp))



X3C_loop <- X3C_loop1 %>%
  bind_rows(X3C_loop3)


ref_gene <- "control"
ref_data <- X3C_loop %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp)%>%
  dplyr::select(-gene)
gene_data <- X3C_loop %>%
  filter(gene != ref_gene)

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         time=="0")%>%
  group_by(gene)%>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_short <- dCt_data %>%
  left_join(dCt_data_control)%>%
  mutate(ddCt = (dCp - dCt_control),
         FC = 2^(-ddCt))

stat_ddCt_data_CL3 <- ddCt_data_short %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))


stat_ddCt_data_1F8R_38 <- stat_ddCt_data_CL3 %>%
  filter(genotype=="rnai38", gene =="CL1_1F8R")

my_graph <- 
  stat_ddCt_data_CL3 %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=gene),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype)
my_graph

# pull 1F8R together -------



stat_ddCt_data_lhp1




all_3C <- stat_ddCt_data_short_1F8R_col %>%
  add_column(genotype ="col")%>%
  bind_rows(stat_ddCt_data_1F8R_92)%>%
  bind_rows(stat_ddCt_data_1F8R_38)

my_graph <- 
  all_3C %>%
  ggplot(.,aes(time,mean), group=gene)+
  geom_line(aes(color=genotype),size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std),width=.2)+
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2)+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  scale_color_manual(breaks=selected_genotype_3C$genotype,
                     values=selected_genotype_3C$genotype_color,
                     labels=selected_genotype_3C$genotype_name,)+ theme(legend.position = "none")
my_graph



ggsave(filename="out/3C/CL_1_1F8R.pdf", plot=my_graph, width = 5, height = 3)











