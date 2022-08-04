library("tidyverse")
library("readxl")
library("readr")
#

source('stat_functions.R')
library(multcompView)


# Data import ------------------------
# Moyenne des Cp r?alis? dans le qPCR_output

plan_plaque_NACL_1 <- read_excel("data/qPCR/in/plan_plaque_NACL_1.xlsx")
X20171226_NACL_Deriv2_Rtech1 <- read_excel("data/qPCR/in/20171226_NACL_Deriv2_Rtech1.xlsx", skip = 1)
X20171226_NACL_Deriv2_Rtech2 <- read_excel("data/qPCR/in/20171226_NACL_Deriv2_Rtech2.xlsx", skip = 1)



#Data fusion & tidying data (ici L1 ler du 72 ?limin? car pas exprim? chez ce genotype)
# + Nom echantillon s?par? en genotype/gene/n?replicat

qPCR_NACL_1 <- X20171226_NACL_Deriv2_Rtech1 %>%
  left_join(plan_plaque_NACL_1) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  select(Cp, genotype, gene) %>%
  separate(genotype, c("genotype", "replicate"), sep="-") %>%
  separate(replicate, c("replicate", "condition"), sep="_") %>%
  rename(Cp1 = Cp)



qPCR_NACL_2 <- X20171226_NACL_Deriv2_Rtech2 %>%
  left_join(plan_plaque_NACL_1) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  select(Cp, genotype, gene) %>%
  separate(genotype, c("genotype", "replicate"), sep="-") %>%
  separate(replicate, c("replicate", "condition"), sep="_") %>%
  rename(Cp2 = Cp)  

qPCR_NACL_1.2 <- qPCR_NACL_1 %>% 
  left_join(qPCR_NACL_2) %>%
  mutate(Cp=(Cp1+Cp2)/2) %>% 
  mutate(CpDiff=(Cp1-Cp2)) %>%
  select(Cp, genotype, gene, replicate, condition)

#Choisir son gene de r?f?rence (choisir celui qui est le moins dispers?; le plus lin?aire, ici ref2)
ggplot(data=qPCR_NACL_1.2, aes(x=genotype, y=Cp, color=replicate)) + geom_point() + facet_wrap(~gene)


#Selection du gene de ref & retirer le gene de ref du dataframe (=gen_data ici)

#Cp mean ref
ref_gene <- "ref2"
ref_data <- qPCR_NACL_1.2 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  select(-gene)

gene_data <- qPCR_NACL_1.2 %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")

# Fusion dataframe (sans le gene de reference) avec Cp du gene de reference pour chaque echantillon (=ref_data)
#& ajout colonne dCp (=variation en nb de cycle/gene par rapport au gene de reference)
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

# Ajout dCp/gene par rapport ? la condition controle (ici Col NoABA) = dCt_data_control
dCt_data_control <- dCt_data %>%
  filter(genotype == "Col") %>%
  filter(condition == "NoNaCl") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

# Calcul du ddCt = Ct gN par rapport ? Ct gN ref coupl? ? Ct condition de ref + calcul du Fold Change (FC)
ddCt_data <- dCt_data %>%
  drop_na()%>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 


#### Graph

ggplot(ddCt_data, aes(x=genotype, y=-ddCt, color=condition)) +
  geom_boxplot() +
  facet_wrap(~gene, nrow = 2, scale="free") +
  theme_bw() 





stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, genotype,condition) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "condition")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

position=position_dodge(.90)

stat_qPCR_cDNA$gene <- factor(stat_qPCR_cDNA$gene, levels = c("UpNPC72","NPC72","DwnNPC72"))
stat_qPCR_cDNA$condition <- factor(stat_qPCR_cDNA$condition, levels = c("NoNaCl","NaCl"))


my_graph <- 
  stat_qPCR_cDNA %>% filter(gene !="NPC48")%>%
  ggplot(data=., mapping=aes(x=genotype, y=-mean,fill=condition)) +
  facet_wrap(~gene)+
  geom_bar(stat="identity", position=position) +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem),width=0.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  geom_hline(yintercept = 0, linetype=1)+
  theme_bw()
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/NPC72 salt.pdf", plot=my_graph, width = 7, height = 3.5)









