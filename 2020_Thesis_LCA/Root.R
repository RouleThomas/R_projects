library(readr)


# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/Root/", showWarnings = FALSE, recursive = TRUE)


#Nadapt -----


germination_LCA_N_adapt <- read_excel("data/Root/Nadapt/germination_LCA_N_adapt.xlsx") %>%
  dplyr::select(-replicate_RNA)
D12 <- read_csv("data/Root/Nadapt/J12.csv") %>%
  dplyr::select(-X6) %>%
  separate(Tag, into = c("plate_nb","plant_nb"), sep = ":") %>%
  mutate(plate_nb=as.numeric(gsub("Nadapt_", "", plate_nb)),
         plant_nb=as.numeric(plant_nb)) %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  left_join(germination_LCA_N_adapt) %>%
  filter(genotype %in% c("col","lca")) %>%
  dplyr::select(genotype, condition,total_root_cm, primary_root_cm,lateral_root_cm,lateral_root_count,lateral_root_sum,density)
  




#Pr ABA -----

stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "primary_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Primary root")

my_graph

# Save 
ggsave(filename="out/Root/Primary_root_Nadapt.pdf", plot=my_graph, width = 7, height = 5)


#LR count -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_count", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)


my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Lateral root count")

my_graph


# Save 
ggsave(filename="out/Root/LR_count_Nadapt.pdf", plot=my_graph, width = 7, height = 5)





#LR lenght -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Mean lenght LR")

my_graph



# Save 
ggsave(filename="out/Root/LR_mean_lenght_Nadapt.pdf", plot=my_graph, width = 7, height = 5)


# density ------


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "density", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Density")

my_graph



# Save 
ggsave(filename="out/Root/density_Nadapt.pdf", plot=my_graph, width = 7, height = 5)


# total root -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "total_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Total root")

my_graph


# Save 
ggsave(filename="out/Root/total_Nadapt.pdf", plot=my_graph, width = 7, height = 5)


# lateral root lenght total -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_sum", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  facet_wrap(~condition, nrow=1)+
  ggtitle("Lateral root lenght")

my_graph



# Save 
ggsave(filename="out/Root/LR_lenght_total_Nadapt.pdf", plot=my_graph, width = 7, height = 5)



# control low light and High light-------

germination_LCA_2 <- read_excel("data/Root/Lowlight/germination_LCA_2.xlsx")

D12 <- read_csv("data/Root/Lowlight/J12-90.csv") %>%
  dplyr::select(-X6) %>%
  separate(Tag, into = c("plate_nb","plant_nb"), sep = ":") %>%
  mutate(plate_nb=as.numeric(gsub("LCA_J12-90-", "", plate_nb)),
         plant_nb=as.numeric(plant_nb)) %>%
  mutate(total_root_cm =`Total Length` / 98, 
         primary_root_cm = `Average Length - Primary roots` / 98, 
         lateral_root_cm =`Average Length - Lateral roots` / 98,
         lateral_root_count = `Lateral Root Count`,
         lateral_root_sum = total_root_cm - primary_root_cm,
         density = lateral_root_count / primary_root_cm) %>%
  filter(primary_root_cm>1) %>%
  left_join(germination_LCA_2) %>%
  filter(genotype %in% c("col","lca")) %>%
  filter(condition =="LowLight") %>%
  dplyr::select(genotype, condition,total_root_cm, primary_root_cm,lateral_root_cm,lateral_root_count,lateral_root_sum,density)





#Pr ABA -----

stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "primary_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Primary root")

my_graph

# Save 
ggsave(filename="out/Root/Primary_root_control.pdf", plot=my_graph, width = 7, height = 5)


#LR count -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_count", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)


my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Lateral root count")

my_graph


# Save 
ggsave(filename="out/Root/LR_count_control.pdf", plot=my_graph, width = 7, height = 5)





#LR lenght -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Mean lenght LR")

my_graph



# Save 
ggsave(filename="out/Root/LR_mean_lenght_control.pdf", plot=my_graph, width = 7, height = 5)


# density ------


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "density", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Density")

my_graph



# Save 
ggsave(filename="out/Root/density_control.pdf", plot=my_graph, width = 7, height = 5)


# total root -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "total_root_cm", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Total root")

my_graph


# Save 
ggsave(filename="out/Root/total_control.pdf", plot=my_graph, width = 7, height = 5)


# lateral root lenght total -----


stat_total_root_length <- 
  D12 %>%
  group_by(condition) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "lateral_root_sum", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

position=position_dodge(.9)

my_graph <- 
  stat_total_root_length %>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - sem, ymax=mean + sem), position=position, width=.5)+
  geom_text(aes(y=(1.1*mean + sem), label=letter), position=position)+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) + 
  theme(legend.position="none")+
  ylab("") +
  xlab("")+
  ggtitle("Lateral root lenght")

my_graph



# Save 
ggsave(filename="out/Root/LR_lenght_total_control.pdf", plot=my_graph, width = 7, height = 5)






