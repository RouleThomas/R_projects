# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/Root/", showWarnings = FALSE, recursive = TRUE)


# data processing

#exp1 ------

plate_description_ABA <- 
  read_excel("data/Root/plate_description_ABA.xlsx")

D9 <- read_csv("data/Root/D9.csv") %>%
  dplyr::select(-X6) %>%
  separate(Tag, into = c("plate_nb","plant_nb"), sep = ":") %>%
  mutate(plate_nb=as.numeric(gsub("plate_", "", plate_nb)),
         plant_nb=as.numeric(plant_nb)) %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  left_join(plate_description_ABA) %>%
  dplyr::select(genotype, condition,total_root_cm, primary_root_cm,lateral_root_cm,lateral_root_count,lateral_root_sum,density)



# exp 2 ------

plate_description_ABA <- 
  read_excel("data/Root/root_info_ABA.xlsx")

D9_2 <- read_csv("data/Root/D9_exp2.csv") %>%
  dplyr::select(-X6) %>%
  separate(Tag, into = c("plate_nb","plant_nb"), sep = ":") %>%
  mutate(plate_nb=as.numeric(gsub("plate_", "", plate_nb)),
         plant_nb=as.numeric(plant_nb)) %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  left_join(plate_description_ABA) %>%
  dplyr::select(genotype, condition,total_root_cm, primary_root_cm,lateral_root_cm,lateral_root_count,lateral_root_sum,density)




# exp1 and exp2 pulled ------

D9_all <- D9 %>%
  bind_rows(D9_2)


# exp 3 drought salt ------

plate_description_salt_drought <- read_excel("data/Root/plate_description_salt_drought.xlsx")



D9_drought_salt <- read_csv("data/Root/salt_drought_measure.csv") %>%
  dplyr::select(-X6) %>%
  separate(Tag, into = c("plate_nb","plant_nb"), sep = ":") %>%
  mutate(plate_nb=as.numeric(gsub("plate_", "", plate_nb)),
         plant_nb=as.numeric(plant_nb)) %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  left_join(plate_description_salt_drought) %>%
  dplyr::select(genotype, condition,total_root_cm, primary_root_cm,lateral_root_cm,lateral_root_count,lateral_root_sum,density)




# pull all exp --------


D9_all_all <- D9_all %>%
  bind_rows(D9_drought_salt)


#Pr ABA -----

stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "primary_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Primary root")

my_graph

# Save 
ggsave(filename="out/Root/Primary_root_ABA.pdf", plot=my_graph, width = 7, height = 5)


#LR count -----


stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_count", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Lateral root count")

my_graph

# Save 
ggsave(filename="out/Root/LR_count_ABA.pdf", plot=my_graph, width = 7, height = 5)





#LR lenght -----


stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Mean lenght LR")

my_graph


# Save 
ggsave(filename="out/Root/LR_mean_lenght_ABA.pdf", plot=my_graph, width = 7, height = 5)


# density ------


stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "density", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Density")

my_graph


# Save 
ggsave(filename="out/Root/density_ABA.pdf", plot=my_graph, width = 7, height = 5)


# total root -----


stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "total_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Total root lenght")

my_graph

# Save 
ggsave(filename="out/Root/total_ABA.pdf", plot=my_graph, width = 7, height = 5)


# lateral root lenght total -----


stat_total_root_length <- 
  D9_all_all %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_sum", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_total_root_length$condition <- factor(stat_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))
position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  left_join(genotype_metadata_root) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  scale_fill_manual(values=genotype_metadata_root$genotype_color)+ 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Lateral root length")

my_graph


# Save 
ggsave(filename="out/Root/LR_lenght_total_ABA.pdf", plot=my_graph, width = 7, height = 5)












#ratio -----




# PR ----

D9_ABA <- D9_all_all %>%
  filter(condition =="ABA") %>%
  dplyr::select(-condition) %>%
  rename(PR_ABA = "primary_root_cm") %>%
  dplyr::select(genotype,PR_ABA)
D9_control <- D9_all_all %>%
  filter(condition =="control") %>%
  dplyr::select(-condition) %>%
  rename(PR_NoABA = "primary_root_cm") %>%
  dplyr::select(genotype,PR_NoABA)

D9_ratio <- D9_ABA %>% 
  left_join(D9_NoABA) %>%
  mutate(ratio= PR_NoABA-PR_ABA) %>%
  filter(ratio>0.25) 


stat_total_root_length <- D9_ratio %>%
  group_by(genotype) %>%
  summarise(mean_ratio=mean(ratio), 
            median_ratio =median(ratio),
            ecart_type=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart



#COMPARAISON PARAMETRIQUE: TESTS ANOVA (comparaison de plusieurs param?tres)#


model_root_length <- aov(ratio ~ genotype, data=D9_ratio)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 



#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)

stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))


my_graph <- 
  stat_all_total_root_length %>% 
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean_ratio, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean_ratio - erreur_std, ymax=mean_ratio + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean_ratio + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  scale_fill_manual(values=genotype_metadata$genotype_color)

my_graph


# Save 
ggsave(filename="out/Root/ratio_PR.pdf", plot=my_graph, width = 5, height = 3)


# Lr lenght total ------
D9_ABA <- D9_all %>%
  filter(condition =="ABA") %>%
  dplyr::select(-condition) %>%
  rename(PR_ABA = "lateral_root_sum") %>%
  dplyr::select(genotype,PR_ABA)
D9_NoABA <- D9_all %>%
  filter(condition =="control") %>%
  dplyr::select(-condition) %>%
  rename(PR_NoABA = "lateral_root_sum") %>%
  dplyr::select(genotype,PR_NoABA)

D9_ratio <- D9_ABA %>% 
  left_join(D9_NoABA) %>%
  mutate(ratio= PR_NoABA-PR_ABA) %>%
  filter(ratio>0.25) 



stat_total_root_length <- D9_ratio %>%
  group_by(genotype) %>%
  summarise(mean_ratio=mean(ratio), 
            median_ratio =median(ratio),
            ecart_type=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart



#COMPARAISON PARAMETRIQUE: TESTS ANOVA (comparaison de plusieurs param?tres)#


model_root_length <- aov(ratio ~ genotype, data=D9_ratio)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 



#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)

stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))


my_graph <- 
  stat_all_total_root_length %>% 
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean_ratio, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean_ratio - erreur_std, ymax=mean_ratio + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean_ratio + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  scale_fill_manual(values=genotype_metadata$genotype_color)

my_graph


# Save 
ggsave(filename="out/Root/ratio_LR_lenght_total.pdf", plot=my_graph, width = 5, height = 3)



# total root -----

D9_ABA <- D9_all %>%
  filter(condition =="ABA") %>%
  dplyr::select(-condition) %>%
  rename(PR_ABA = "total_root_cm") %>%
  dplyr::select(genotype,PR_ABA)
D9_NoABA <- D9_all %>%
  filter(condition =="control") %>%
  dplyr::select(-condition) %>%
  rename(PR_NoABA = "total_root_cm") %>%
  dplyr::select(genotype,PR_NoABA)

D9_ratio <- D9_ABA %>% 
  left_join(D9_NoABA) %>%
  mutate(ratio= PR_NoABA-PR_ABA) %>%
  filter(ratio>0.25) 



stat_total_root_length <- D9_ratio %>%
  group_by(genotype) %>%
  summarise(mean_ratio=mean(ratio), 
            median_ratio =median(ratio),
            ecart_type=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart



#COMPARAISON PARAMETRIQUE: TESTS ANOVA (comparaison de plusieurs param?tres)#


model_root_length <- aov(ratio ~ genotype, data=D9_ratio)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 



#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)

stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))


my_graph <- 
  stat_all_total_root_length %>% 
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean_ratio, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean_ratio - erreur_std, ymax=mean_ratio + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean_ratio + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  scale_fill_manual(values=genotype_metadata$genotype_color)

my_graph


# Save 
ggsave(filename="out/Root/ratio_lenght_total.pdf", plot=my_graph, width = 5, height = 3)


# density ------



D9_ABA <- D9_all %>%
  filter(condition =="ABA") %>%
  dplyr::select(-condition) %>%
  rename(PR_ABA = "density") %>%
  dplyr::select(genotype,PR_ABA)
D9_NoABA <- D9_all %>%
  filter(condition =="control") %>%
  dplyr::select(-condition) %>%
  rename(PR_NoABA = "density") %>%
  dplyr::select(genotype,PR_NoABA)

D9_ratio <- D9_ABA %>% 
  left_join(D9_NoABA) %>%
  mutate(ratio= PR_NoABA-PR_ABA) %>%
  filter(ratio>0.25) 



stat_total_root_length <- D9_ratio %>%
  group_by(genotype) %>%
  summarise(mean_ratio=mean(ratio), 
            median_ratio =median(ratio),
            ecart_type=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart



#COMPARAISON PARAMETRIQUE: TESTS ANOVA (comparaison de plusieurs param?tres)#


model_root_length <- aov(ratio ~ genotype, data=D9_ratio)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 



#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)

stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))


my_graph <- 
  stat_all_total_root_length %>% 
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean_ratio, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean_ratio - erreur_std, ymax=mean_ratio + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean_ratio + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  scale_fill_manual(values=genotype_metadata$genotype_color)

my_graph


# Save 
ggsave(filename="out/Root/ratio_density.pdf", plot=my_graph, width = 5, height = 3)

# mean leanght LR ------




D9_ABA <- D9_all %>%
  filter(condition =="ABA") %>%
  dplyr::select(-condition) %>%
  rename(PR_ABA = "lateral_root_cm") %>%
  dplyr::select(genotype,PR_ABA)
D9_NoABA <- D9_all %>%
  filter(condition =="control") %>%
  dplyr::select(-condition) %>%
  rename(PR_NoABA = "lateral_root_cm") %>%
  dplyr::select(genotype,PR_NoABA)

D9_ratio <- D9_ABA %>% 
  left_join(D9_NoABA) %>%
  mutate(ratio= PR_NoABA-PR_ABA) %>%
  filter(ratio>0.25) 



stat_total_root_length <- D9_ratio %>%
  group_by(genotype) %>%
  summarise(mean_ratio=mean(ratio), 
            median_ratio =median(ratio),
            ecart_type=sd(ratio), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart



#COMPARAISON PARAMETRIQUE: TESTS ANOVA (comparaison de plusieurs param?tres)#


model_root_length <- aov(ratio ~ genotype, data=D9_ratio)
anova(model_root_length)



# Comparaison 2 a 2 toutes d'un coup
TukeyHSD(model_root_length)

# Resume avec padj < .05 pour (genotype:condition)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))

resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 



#croiser resume_letter_df et stat par "genotype" et "condition"
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)

stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("control","ABA", "mannitol","salt"))


my_graph <- 
  stat_all_total_root_length %>% 
  left_join(genotype_metadata) %>% 
  ggplot(data=., mapping=aes(x=genotype_name, y=mean_ratio, fill=genotype_name)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean_ratio - erreur_std, ymax=mean_ratio + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean_ratio + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  scale_fill_manual(values=genotype_metadata$genotype_color)

my_graph


# Save 
ggsave(filename="out/Root/ratio_LR_mean.pdf", plot=my_graph, width = 5, height = 3)






