# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/Root phenotype/", showWarnings = FALSE, recursive = TRUE)

#exp1 SALK_140401 ----

#Data import
germination_IAA_exp1 <- read_excel("data/exp1_MP/germination_IAA_exp1.xlsx")
Plant_measurements_exp1 <- read_csv("data/exp1_MP/Plant measurements_exp1.csv") %>%
  separate(Tag, c("plate", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(gsub("iaa_", "", plate)),
         plant_nb=as.numeric(plant_nb)) %>%
  dplyr::select(-plate)

# Processing data

measures_J11 <- Plant_measurements_exp1 %>%
  left_join(germination_IAA_exp1)

full_data <- measures_J11 %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  filter(genotype %in% c("col","SALK_140401"))




# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(total_root_cm), 
            (mediane_root_length_cm =median(total_root_cm)),
            ecart_type=sd(total_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph


# Save 
ggsave(filename="out/Root phenotype/Total_root_SALK140401.pdf", plot=my_graph, width = 5, height = 3)




#Lateral root
stat_lateral_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_cm), 
            (mediane_root_length_cm =median(lateral_root_cm)),
            ecart_type=sd(lateral_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")+
  theme_bw()
my_graph


# Save 
ggsave(filename="out/Root phenotype/LR_lenght_SALK140401.pdf", plot=my_graph, width = 3, height = 3)


#Primary root
stat_primary_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_primary_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph


# Save 
ggsave(filename="out/Root phenotype/PR_SALK140401.pdf", plot=my_graph, width = 5, height = 3)




#Density
stat_density <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(density), 
            (mediane_root_length_cm =median(density)),
            ecart_type=sd(density), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(density ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_density %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Density")+
  xlab(label = "")+
  ylab(label = "count")+
  theme_bw()
my_graph


# Save 
ggsave(filename="out/Root phenotype/density_SALK140401.pdf", plot=my_graph, width = 3, height = 3)






#exp3_useless----

#Data import
germination_IAA_exp3 <- read_excel("data/exp1_MP/germination_IAA_exp3.xlsx")
Plant_measurements_exp3 <- read_csv("data/exp1_MP/Plant measurements_exp3.csv") %>%
  dplyr::select(-X6) %>% 
  separate(Tag, c("plate", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(gsub("iaa_", "", plate)),
         plant_nb=as.numeric(plant_nb)) %>%
  dplyr::select(-plate)

# Processing data

measures_J11 <- germination_IAA_exp3 %>%
  left_join(Plant_measurements_exp3)

full_data <- measures_J11 %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1)



# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(total_root_cm), 
            (mediane_root_length_cm =median(total_root_cm)),
            ecart_type=sd(total_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")

#Lateral root
stat_lateral_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_cm), 
            (mediane_root_length_cm =median(lateral_root_cm)),
            ecart_type=sd(lateral_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")


#Primary root
stat_primary_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_primary_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")


#Density
stat_density <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(density), 
            (mediane_root_length_cm =median(density)),
            ecart_type=sd(density), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(density ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_density %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Density")+
  xlab(label = "")+
  ylab(label = "count")


#exp4 ----

#Data import
germination_IAA_exp4 <- read_excel("data/exp1_MP/germination_IAA_exp4.xlsx")
Plant_measurements_exp4 <- read_csv("data/exp1_MP/Plant measurements_exp4.csv") %>%
  separate(Tag, c("plate", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(gsub("iaa_", "", plate)),
         plant_nb=as.numeric(plant_nb)) %>%
  dplyr::select(-plate)

# Processing data

measures_J11 <- germination_IAA_exp4 %>%
  left_join(Plant_measurements_exp4)

full_data <- measures_J11 %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1)



# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(total_root_cm), 
            (mediane_root_length_cm =median(total_root_cm)),
            ecart_type=sd(total_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))

#graph 
position=position_dodge(.9)
my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/total_root_SALK113294.pdf", plot=my_graph, width = 5, height = 3)


#Lateral root
stat_lateral_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_cm), 
            (mediane_root_length_cm =median(lateral_root_cm)),
            ecart_type=sd(lateral_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")+
  theme_bw()
my_graph

# Save 
ggsave(filename="out/Root phenotype/LR_SALK113294.pdf", plot=my_graph, width = 3, height = 3)



#Primary root
stat_primary_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_primary_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph


# Save 
ggsave(filename="out/Root phenotype/PR_SALK113294.pdf", plot=my_graph, width = 3, height = 3)



#Density
stat_density <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(density), 
            (mediane_root_length_cm =median(density)),
            ecart_type=sd(density), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(density ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_density %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Density")+
  xlab(label = "")+
  ylab(label = "count")+
  theme_bw()
my_graph


# Save 
ggsave(filename="out/Root phenotype/density_SALK113294.pdf", plot=my_graph, width = 3, height = 3)


#exp1 SALK_140401 with NAA day12 -----

#Data import
germination_IAA_exp1NAA <- read_excel("data/exp1_MP/germination_IAA_exp1NAA.xlsx")
germination_IAA_exp1NAA_nutt <- read_excel("data/exp1_MP/germination_IAA_exp1NAA_nutt.xlsx")


Plant_measurements_exp1NAA <- read_excel("data/exp1_MP/Plant measurements_exp1NAA.xlsx") %>% 
  separate(Tag, c("plate", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(gsub("naa", "", plate)),
         plant_nb=as.numeric(plant_nb)) %>%
  dplyr::select(-plate)

Plant_measurements_D6 <- read_csv("data/exp1_MP/D6.csv") %>%
  dplyr::select(-X6) %>% 
  separate(Tag, c("plate_nb", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(plate_nb), plant_nb=as.numeric((plant_nb)))

Plant_measurements_D9 <- read_csv("data/exp1_MP/D9.csv") %>%
  dplyr::select(-X3) %>% 
  separate(Tag, c("plate_nb", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(plate_nb), plant_nb=as.numeric((plant_nb)))



# Processing data

measures_J11 <- germination_IAA_exp1NAA %>%
  left_join(Plant_measurements_exp1NAA) %>%
  mutate(day="12")
measures_J6 <- germination_IAA_exp1NAA_nutt %>%
  left_join(Plant_measurements_D6) %>%
  mutate(day="6")
measures_J9 <- germination_IAA_exp1NAA_nutt %>%
  left_join(Plant_measurements_D9) %>%
  mutate(day="9", `Average Length - Primary roots`=`Average Length - Primary roots`/3)
  
measures <- measures_J11 %>%
  bind_rows(measures_J6,measures_J9)

measures <- measures_J11


full_data <- measures %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>0.5)

#full_data$day <- factor(full_data$day, levels=c("6","12")) # Choisir ordrer du facte_wrap



# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype, condition) %>%
  summarise(moyenne_root_length_cm=mean(total_root_cm), 
            (mediane_root_length_cm =median(total_root_cm)),
            ecart_type=sd(total_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_root_cm ~ genotype*condition, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype:condition`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(genotype, into= c("genotype","condition"), sep = ":")

stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype", "condition"))
#graph 
#stat_all_total_root_length$day <- factor(stat_all_total_root_length$day, levels=c("6","12")) # Choisir ordrer du facte_wrap
stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("NoNAA","NAA")) # Choisir ordrer du facte_wrap

position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total root length")+
  xlab(label = "")+
  ylab(label = "length (cm)") +
  facet_wrap(~condition)
my_graph


# Save 
ggsave(filename="out/Root phenotype/TR_NAA_SALK140401.pdf", plot=my_graph, width = 5, height = 3)




#Lateral root
stat_lateral_root_length <- full_data %>%
  group_by(genotype, condition) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_cm), 
            (mediane_root_length_cm =median(lateral_root_cm)),
            ecart_type=sd(lateral_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_cm ~ genotype*condition, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype:condition`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(genotype, into= c("genotype","condition"), sep = ":")
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype","condition"))
#graph 
stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("NoNAA","NAA")) # Choisir ordrer du facte_wrap

position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)") +
  facet_wrap(~condition)
my_graph

# Save 
ggsave(filename="out/Root phenotype/LR_NAA_SALK140401.pdf", plot=my_graph, width = 5, height = 3)


#Primary root
stat_primary_root_length <- full_data %>%
  group_by(genotype,condition) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype*condition, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype:condition:day`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters)  %>%
  separate(genotype, into= c("genotype","condition"), sep = ":")
stat_all_total_root_length <- stat_primary_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype","condition"))
#graph 
stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("NoNAA","NAA")) # Choisir ordrer du facte_wrap

position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)") +
  facet_wrap(~condition)
my_graph

# Save 
ggsave(filename="out/Root phenotype/PR_NAA_SALK140401.pdf", plot=my_graph, width = 5, height = 3)


#Density
stat_density <- full_data %>%
  group_by(genotype,condition) %>%
  summarise(moyenne_root_length_cm=mean(density), 
            (mediane_root_length_cm =median(density)),
            ecart_type=sd(density), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(density ~ genotype*condition, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype:condition`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(genotype, into= c("genotype","condition"), sep = ":")
stat_all_total_root_length <- stat_density %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype","condition"))
#graph 
stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("NoNAA","NAA")) # Choisir ordrer du facte_wrap

position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Density")+
  xlab(label = "")+
  ylab(label = "count") +
  facet_wrap(~condition)
my_graph

# Save 
ggsave(filename="out/Root phenotype/density_NAA_SALK140401.pdf", plot=my_graph, width = 5, height = 3)


# Exp4 Nutt RNAi control THE ONE USE IN PAPER for RNAi DAY 11 ------

#Data import

germination_EXP4 <- read_excel("data/exp1_MP/germination_EXP4.xlsx")
Plant_measurements_D11 <- read_csv("data/exp1_MP/D11.csv") %>%
  separate(Tag, c("plate", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(gsub("iaa_", "", plate)),
         plant_nb=as.numeric(plant_nb)) %>%
  dplyr::select(-plate)


# Processing data

measures_J11 <- germination_EXP4 %>%
  left_join(Plant_measurements_D11)

full_data <- measures_J11 %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1)



# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(total_root_cm), 
            (mediane_root_length_cm =median(total_root_cm)),
            ecart_type=sd(total_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/TR_RNAi.pdf", plot=my_graph, width = 5, height = 3)


#Lateral root
stat_lateral_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_cm), 
            (mediane_root_length_cm =median(lateral_root_cm)),
            ecart_type=sd(lateral_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/LR_RNAi.pdf", plot=my_graph, width = 3.5, height = 3)





#Total Lateral root (Total - Primary)
stat_lateral_root_length <- full_data %>%
  mutate(total_LR_root=total_root_cm-primary_root_cm) %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(total_LR_root), 
            (mediane_root_length_cm =median(total_LR_root)),
            ecart_type=sd(total_LR_root), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(total_LR_root ~ genotype, data=full_data %>%
                           mutate(total_LR_root=total_root_cm-primary_root_cm) )
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
  ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/Total_LR_RNAi.pdf", plot=my_graph, width = 3.5, height = 3)




#Total Lateral root (Total - Primary)
stat_lateral_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(lateral_root_sum), 
            (mediane_root_length_cm =median(lateral_root_sum)),
            ecart_type=sd(lateral_root_sum), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(lateral_root_sum ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_lateral_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
  ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Total lateral root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/Total_LR_RNAi.pdf", plot=my_graph, width = 3.5, height = 3)









#Primary root
stat_primary_root_length <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_primary_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)")
my_graph

# Save 
ggsave(filename="out/Root phenotype/PR_RNAi.pdf", plot=my_graph, width = 3.5, height = 3)



#Density
stat_density <- full_data %>%
  group_by(genotype) %>%
  summarise(moyenne_root_length_cm=mean(density), 
            (mediane_root_length_cm =median(density)),
            ecart_type=sd(density), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(density ~ genotype, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) 
stat_all_total_root_length <- stat_density %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype"))
#graph 
position=position_dodge(.9)

my_graph <- 
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Density")+
  xlab(label = "")+
  ylab(label = "count")
my_graph



# Save 
ggsave(filename="out/Root phenotype/density_RNAi.pdf", plot=my_graph, width = 3.5, height = 3)



# Clean tidy stat plots with pvalue -----



#Total Lateral root (Total - Primary) * PLOT stat other ------
lateral_root_length <- full_data %>%
  mutate(total_LR_root=total_root_cm-primary_root_cm) %>%
  select(genotype,total_LR_root)

lateral_root_length%>%
ggbarplot(., x = "genotype", y = "total_LR_root",add = "mean_se") +
  stat_compare_means(method="wilcox.test", 
                     comparisons = list(c("col", "RNAi 1.2"), c("col", "RNAi 3.2"), c("col", "RNAi 6.2") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() + 
  theme(legend.position="none")+
  ylab("length (cm)")+
  xlab("")


my_graph


# Save 
ggsave(filename="C:/Users/roule/Box/Ongoing projects/ARES/R nucleus/out/nucleus_ARES.pdf", plot=my_graph, width = 3, height = 3) 




p.adjust(c(0.11,0.32,0.071), method="fdr")
























# NAA RNAi 1.2 and SALK140401_useless walah or need to draw the other root parameters -----

germination_file_rnai_salk <- read_excel("data/exp1_MP/germination_file_rnai_salk.xlsx")
D11_rnai_salk <- read_csv("data/exp1_MP/D11_rnai_salk.csv") %>%
  dplyr::select(-X3) %>% 
  separate(Tag, c("plate_nb", "plant_nb"), sep=":") %>%
  mutate(plate_nb=as.numeric(plate_nb), plant_nb=as.numeric((plant_nb)))


# Processing data

measures_J11 <- germination_file_rnai_salk %>%
  left_join(D11_rnai_salk)


full_data <- measures_J11 %>%
  mutate(primary_root_cm = `Average Length - Primary roots` / 98) %>%
  filter(primary_root_cm>0.5)


# Stat
#Total root
stat_total_root_length <- full_data %>%
  group_by(genotype, condition) %>%
  summarise(moyenne_root_length_cm=mean(primary_root_cm), 
            (mediane_root_length_cm =median(primary_root_cm)),
            ecart_type=sd(primary_root_cm), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

model_root_length <- aov(primary_root_cm ~ genotype*condition, data=full_data)
anova(model_root_length)

TukeyHSD(model_root_length)
resume_lettre_total_root_length <- multcompLetters(extract_p(TukeyHSD(model_root_length)$`genotype:condition`))
resume_lettre_df_total_root_length <- data.frame(genotype=names(resume_lettre_total_root_length$Letters),
                                                 letter=resume_lettre_total_root_length$Letters) %>%
  separate(genotype, into= c("genotype","condition"), sep = ":")

stat_all_total_root_length <- stat_total_root_length %>% 
  left_join(resume_lettre_df_total_root_length, by=c("genotype", "condition"))
#graph 
stat_all_total_root_length$condition <- factor(stat_all_total_root_length$condition, levels=c("NoNAA","NAA")) # Choisir ordrer du facte_wrap

position=position_dodge(.9)
ggplot(data=stat_all_total_root_length, mapping=aes(x=genotype, y=moyenne_root_length_cm)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=moyenne_root_length_cm - erreur_std, ymax=moyenne_root_length_cm + erreur_std), position=position, width=.5)+
  geom_text(aes(y=(1.1*moyenne_root_length_cm + erreur_std), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=45, vjust=.5, size=10))+
  ggtitle(label = "Primary root length")+
  xlab(label = "")+
  ylab(label = "length (cm)") +
  facet_wrap(~condition) 











