# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/FAIRE/", showWarnings = FALSE, recursive = TRUE)

# Data import
FAIRE_qPCR_output <- read_excel("data/FAIRE/FAIRE_qPCR_output.xlsx",sheet=2) %>% clean_genotype_code

# Data processing
qPCR_FAIRE_input <- FAIRE_qPCR_output %>% filter(condition == "input") %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(Cp_input, time, replicate, position, genotype)

qPCR_FAIRE_sample <- FAIRE_qPCR_output %>% filter(condition == "sample") %>%
  mutate(Cp_sample=Cp) %>%
  dplyr::select(Cp_sample, time, replicate, position, genotype)

qPCR_FAIRE_Delta <- qPCR_FAIRE_input %>% left_join(qPCR_FAIRE_sample) %>%
  mutate(Delta=2^-(Cp_sample-Cp_input)*100,
         time=as.character(time)) %>% 
  unique()



#Col time 0 vs treated----
qPCR_FAIRE_Delta$time <- factor(qPCR_FAIRE_Delta$time, c("0","30 min","1 hour","4 hours"))
qPCR_FAIRE_Delta$position <- factor(qPCR_FAIRE_Delta$position, c("CYP705A12.3'","CYP705A12.5'", "CYP71A16.3'","CYP71A16.5'","intergenic1","intergenic2","intergenic3",
                                                             "intergenic4","MHAL.5'","MHAL.3'", "MRN1.promotor", "MRN1.5'","MRN1.3'"))

#30min
qPCR_FAIRE_Delta_col <- qPCR_FAIRE_Delta %>%
  filter(genotype=="col", time %in% c("0","30 min"))

stat_dCt_data_col <- qPCR_FAIRE_Delta_col %>%
  group_by(time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_col, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 15) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_col_0vs30min.pdf", plot=my_graph, width = 5, height = 3)

#1hour
qPCR_FAIRE_Delta_col <- qPCR_FAIRE_Delta %>%
  filter(genotype=="col", time %in% c("0","1 hour"))

stat_dCt_data_col <- qPCR_FAIRE_Delta_col %>%
  group_by(coordinate, time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_col, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.signif", label.y = 24) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_col_0vs1hour.pdf", plot=my_graph, width = 5, height = 3)

#30min & 1 hour
qPCR_FAIRE_Delta_col <- qPCR_FAIRE_Delta %>%
  filter(genotype=="col", time %in% c("0","30 min", "1 hour"))

stat_dCt_data_col <- qPCR_FAIRE_Delta_col %>%
  group_by(time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_col, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey") +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_col_0vs30min&1hour&4hours.pdf", plot=my_graph, width = 5, height = 3)

#4hour
qPCR_FAIRE_Delta_col <- qPCR_FAIRE_Delta %>%
  filter(genotype=="col", time %in% c("0","4 hours"))

stat_dCt_data_col <- qPCR_FAIRE_Delta_col %>%
  group_by(coordinate, time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_col, x = "position", y = "Delta", add = "mean_se",
          fill = "time", 
          position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 14.5) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_col_0vs4h.pdf", plot=my_graph, width = 5, height = 3)

#Col vs RNAi time 0----
qPCR_FAIRE_Delta_RNAi_t0 <- qPCR_FAIRE_Delta %>%
  filter(time %in% c("0"))

stat_dCt_data_col <-
    qPCR_FAIRE_Delta_RNAi_t0 %>%
    group_by(position) %>%
    nest %>%
    mutate(stat = map(data, one_way_anova, "Delta", "genotype")) %>%
    dplyr::select(-data) %>%
    unnest(stat)

# Statistics and Graph
# Anova

# Calculate the shit of letter based on the maximum value
label_y_shift <-
  max(stat_dCt_data_col$mean + stat_dCt_data_col$sem) * 0.05
my_graph <- 
  ggplot(stat_dCt_data_col, aes(x = position, y = mean, fill = genotype)) +
  geom_col(position = position_dodge(0.9), color="black")+
  geom_errorbar(aes(ymin = mean - sem, ymax=mean + sem),
                position = position_dodge(0.9), width=.2)+
  geom_text(aes(label = letter, y = mean + sem + label_y_shift),
            position = position_dodge(0.9), size=4) +
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")+theme(legend.position="none")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_colvsRNAi.pdf", plot=my_graph, width = 5, height = 3)

# RNAi 3.8 time 0 vs treated ----

qPCR_FAIRE_Delta_rnai38 <- qPCR_FAIRE_Delta %>%
  filter(genotype=="rnai38", time %in% c("0","4 hours"))

stat_dCt_data_rnai38 <- qPCR_FAIRE_Delta_rnai38 %>%
  group_by(coordinate, time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_rnai38, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 25) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_rnai38_0vs4h.pdf", plot=my_graph, width = 5, height = 3)

# RNAi 9.2 time 0 vs treated ----

qPCR_FAIRE_Delta_rnai92 <- qPCR_FAIRE_Delta %>%
  filter(genotype=="rnai92", time %in% c("0","4 hours"))

stat_dCt_data_rnai92 <- qPCR_FAIRE_Delta_rnai92 %>%
  group_by(time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova
my_graph <- 
  ggbarplot(qPCR_FAIRE_Delta_rnai92, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 35) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")
my_graph

# Save 
ggsave(filename="out/FAIRE/FAIRE_rnai92_0vs4h.pdf", plot=my_graph, width = 5, height = 3)


#Graph/position genotype per genotype FAIRE -----
#
qPCR_FAIRE_intergenic2 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic2", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic2, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic2.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic4 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic4", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic4, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 36) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic4.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MRN1p5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1p5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_MRN1p5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MRN1p3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1p3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_MRN1p3.pdf", plot=my_graph, width = 2.5, height = 2)
#
#
qPCR_FAIRE_MRN1pro <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.promotor", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1pro, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_MRN1pP.pdf", plot=my_graph, width = 2.5, height = 2)
#
#
qPCR_FAIRE_intergenic1 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic1", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic1, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic1.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic3 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic3", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 34) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic4 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic4", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic4, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 35) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic4.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MHAL5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MHAL5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 25) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_MHAL5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MHAL3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MHAL3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 22) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_MHAL3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP71_3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP71_3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 22) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP71_3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP71_5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP71_5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 22) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP71_5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP705_3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP705_3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 30) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP705_3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP705_5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP705_5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 25) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP705_5.pdf", plot=my_graph, width = 2.5, height = 2)


#Graph/position rnai 38 and 92 FAIRE -----
#
qPCR_FAIRE_intergenic2 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic2", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_intergenic2 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_intergenic2)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)



my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph


ggsave(filename="out/FAIRE/FAIRE_intergenic2_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_intergenic4 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic4", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_intergenic4 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_intergenic4)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph


ggsave(filename="out/FAIRE/FAIRE_intergenic4_all.pdf", plot=my_graph, width = 2, height = 3)

#
qPCR_FAIRE_MRN1p5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_MRN1p5 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_MRN1p5)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph


ggsave(filename="out/FAIRE/FAIRE_MRN1p5_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_MRN1p3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_MRN1p3 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_MRN1p3)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph


ggsave(filename="out/FAIRE/FAIRE_MRN1p3_all.pdf", plot=my_graph, width = 2, height = 3)
#
#
qPCR_FAIRE_MRN1pro <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.promotor", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_MRN1pro %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_MRN1pro)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph


ggsave(filename="out/FAIRE/FAIRE_MRN1pP_all.pdf", plot=my_graph, width = 2, height = 3)
#
#
qPCR_FAIRE_intergenic1 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic1", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_intergenic1 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_intergenic1)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic1_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_intergenic3 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic3", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_intergenic3 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_intergenic3)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic3_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_intergenic4 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic4", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_intergenic4 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_intergenic4)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_intergenic4_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_MHAL5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_MHAL5 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_MHAL5)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_MHAL5_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_MHAL3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_MHAL3 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_MHAL3)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_MHAL3_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_CYP71_3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_CYP71_3 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_CYP71_3)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP71_3_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_CYP71_5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_CYP71_5 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_CYP71_5)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP71_5_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_CYP705_3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.3'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_CYP705_3 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_CYP705_3)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP705_3_all.pdf", plot=my_graph, width = 2, height = 3)
#
qPCR_FAIRE_CYP705_5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.5'", time %in% c("0","4 hours"), genotype %in% c("col","rnai92","rnai38")) %>%
  dplyr::select(time,genotype,Delta) %>%
  mutate(time =as.character(time))


stat_total_root_length <- qPCR_FAIRE_CYP705_5 %>%
  group_by(time, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n))  



model_root_length <- aov(Delta ~ time*genotype, data=qPCR_FAIRE_CYP705_5)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`time:genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(col=genotype, into=c("time", "genotype"), sep=":")




#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("time","genotype")) 

#graph 
position=position_dodge(.9)


my_graph <- 
  stat_all_total_root_length %>%
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position, colour="black") +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")

my_graph

ggsave(filename="out/FAIRE/FAIRE_CYP705_5_all.pdf", plot=my_graph, width = 2, height = 3)

















# lhp1 ABA ----
FAIRE_qPCR_output_lhp1 <- read_excel("data/FAIRE/output_FAIRE_lhp1.xlsx") %>% clean_genotype_code


# Data processing
qPCR_FAIRE_input <- FAIRE_qPCR_output_lhp1 %>% filter(condition == "input") %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(Cp_input, time, replicate, position, genotype)

qPCR_FAIRE_sample <- FAIRE_qPCR_output_lhp1 %>% filter(condition == "sample") %>%
  mutate(Cp_sample=Cp) %>%
  dplyr::select(Cp_sample, time, replicate, position, genotype)

qPCR_FAIRE_Delta <- qPCR_FAIRE_input %>% left_join(qPCR_FAIRE_sample) %>%
  mutate(Delta=2^-(Cp_sample-Cp_input)*100,
         time=as.character(time)) %>% 
  unique()



#Col time 0 vs treated----

qPCR_FAIRE_Delta_genotype <- qPCR_FAIRE_Delta %>%
  filter(genotype=="col")

stat_dCt_data_col <- qPCR_FAIRE_Delta_col %>%
  group_by(time, position, genotype) %>%
  summarise(mean=mean(Delta), 
            median= median(Delta),
            SD=sd(Delta), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

# Statistics and Graph
# Anova

ggbarplot(qPCR_FAIRE_Delta_genotype, x = "position", y = "Delta", add = "mean_se",
            fill = "time", 
            position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = time), method="anova", label = "p.format", label.y = 15) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")



# col vs lhp1 ----

qPCR_FAIRE_Delta_time <- qPCR_FAIRE_Delta %>%
  filter(time =="0")




# Statistics and Graph
# Anova

ggbarplot(qPCR_FAIRE_Delta_time, x = "position", y = "Delta", add = "mean_se",
          fill = "genotype", 
          position = position_dodge(0.8), palette="grey")+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 50) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))+
  xlab("")+
  ylab("% input")

#
qPCR_FAIRE_intergenic2 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic2", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic2, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 17) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_intergenic2.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MRN1pP <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.promotor", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1pP, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_MRN1pP.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MRN1p5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.5'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1p5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 18) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_MRN1p5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MRN1p3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MRN1.3'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MRN1p3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 15) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_MRN1p3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP71A16p3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.3'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP71A16p3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 21) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_CYP71A16p3.pdf", plot=my_graph, width = 2.5, height = 2)

#
qPCR_FAIRE_CYP71A16p5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP71A16.5'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP71A16p5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 17) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_CYP71A16p5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP705p3 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.3'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP705p3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_CYP705p3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_CYP705p5 <- qPCR_FAIRE_Delta %>%
  filter(position =="CYP705A12.5'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_CYP705p5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 18) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_CYP705p5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic1 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic1", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic1, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 23) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_intergenic1.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic3 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic3", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 25) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_intergenic3.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_intergenic4 <- qPCR_FAIRE_Delta %>%
  filter(position =="intergenic4", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_intergenic4, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 22.5) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_intergenic4.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MHALp5 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.5'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MHALp5, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 18) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_MHALp5.pdf", plot=my_graph, width = 2.5, height = 2)
#
qPCR_FAIRE_MHALp3 <- qPCR_FAIRE_Delta %>%
  filter(position =="MHAL.3'", time %in% c("0","4 hours"))

my_graph <- 
  ggbarplot(qPCR_FAIRE_MHALp3, x = "time", y = "Delta", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 20) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("% input")+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) 
my_graph

ggsave(filename="out/FAIRE/FAIRE_lhp1_MHALp3.pdf", plot=my_graph, width = 2.5, height = 2)
