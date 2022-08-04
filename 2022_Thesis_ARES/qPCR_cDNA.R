# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

library("multcompView")

# Out directory if it does not exist
dir.create("out/qPCR cDNA/", showWarnings = FALSE, recursive = TRUE)

#exp1 SALK_140401 ----

#Data import
qPCR_output_exp1 <- read_excel("data/qPCR_cDNA/qPCR_output_exp1.xlsx")


ref_gene <- "ref2"
ref_data <- qPCR_output_exp1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_exp1 %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

ddCt_data_exp1 <- ddCt_data %>% filter(genotype %in% c("col","SALK140401"))

# Graph
ddCt_data %>% filter(genotype %in% c("col","SALK140401")) %>%
ggplot(., aes(x=gene, y=-ddCt)) +
  geom_boxplot(aes(color=genotype)) +
  geom_jitter(aes(color=genotype), width=0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.3))


stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 

position=position_dodge(1)

my_graph <- 
  stat_qPCR_cDNA %>% filter(genotype %in% c("col","SALK140401")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


#exp3 SALK_143569 ----
qPCR_output_exp3 <- read_excel("data/qPCR_cDNA/qPCR_output_exp3.xlsx")


ref_gene <- "ref1"
ref_data <- qPCR_output_exp3 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_exp3 %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

ddCt_data_exp3 <- ddCt_data


# Graph
stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 

position=position_dodge(1)

my_graph <- 
  stat_qPCR_cDNA %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


#exp4 SALK_113294 ----
qPCR_output_exp4 <- read_excel("data/qPCR_cDNA/qPCR_output_exp4.xlsx")


ref_gene <- "ref2"
ref_data <- qPCR_output_exp4 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_exp4 %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

ddCt_data_exp4 <- ddCt_data


# Graph
stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 

position=position_dodge(1)

my_graph <- 
  stat_qPCR_cDNA %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph

#Combine exp 1 3 and 4 -----

ddCt_data_134 <- ddCt_data_exp3 %>% 
  bind_rows(ddCt_data_exp4) %>% 
  filter(genotype !="col") %>%
  bind_rows(ddCt_data_exp1) %>%
  filter(! gene %in% c("ref1","ref2"))


stat_qPCR_cDNA <- 
  ddCt_data_134 %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 

position=position_dodge(.8)

my_graph <- 
  stat_qPCR_cDNA %>% filter(genotype !="SALK143569") %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity", width=.8) +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph

#Combine exp1 and 4 (cause exp 3 is useless) ------

#Combine exp 1 3 and 4 -----

ddCt_data_14 <- ddCt_data_exp4 %>%
  filter(genotype !="col") %>%
  bind_rows(ddCt_data_exp1) %>%
  filter(! gene %in% c("ref1","ref2"))


stat_qPCR_cDNA <- 
  ddCt_data_14 %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 

position=position_dodge(.8)

my_graph <- 
  stat_qPCR_cDNA %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity", width=.8) +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/insertional mutant control.pdf", plot=my_graph, width = 5, height = 3)




# Exp1 Nutt AIA concentration Col -------

#Data import
qPCR_output_NAA_EXP1 <- read_excel("data/qPCR_cDNA/qPCR_output_NAA_EXP1.xlsx") %>%
  mutate(Cp=as.numeric(Cp), replicate=as.character(replicate))

#Data processing

ref_gene <- "ref2"
ref_data <- qPCR_output_NAA_EXP1 %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_NAA_EXP1 %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(time == 0, concentration == 10) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

#norm conc
ddCt_data_0 <- ddCt_data %>% filter(concentration=="0")
ddCt_data_0.1 <- ddCt_data %>% filter(concentration=="0.1")
ddCt_data_0.5 <- ddCt_data %>% filter(concentration=="0.5")
ddCt_data_1 <- ddCt_data %>% filter(concentration=="1")
ddCt_data_10 <- ddCt_data %>% filter(concentration=="10")

ddCt_data_all <- ddCt_data_0 %>%
  bind_rows(ddCt_data_0.1,ddCt_data_0.5,ddCt_data_1,ddCt_data_10)

stat_ddCt_data <- ddCt_data_all %>%
  group_by(gene,time,concentration) %>%
  summarise(mean=mean(-ddCt), 
            (median=median(-ddCt)),
            ecart_type=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

my_graph <- 
stat_ddCt_data %>% filter(concentration %in% c("10"), ! gene %in% c("PID","ref1")) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/col_solinc_NAA.pdf", plot=my_graph, width = 5, height = 3)


my_graph <- 
  stat_ddCt_data %>% filter(concentration %in% c("10"), ! gene %in% c("PID","ref1","AT4G14560")) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  xlab(label = "")+
  ylab(label = "log2(FC)")
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/col_solinc_IA14only_NAA.pdf", plot=my_graph, width = 5, height = 3)





# Mutant epigenetique roots -----
qPCR_output_epigenet_root <- read_excel("data/qPCR_cDNA/qPCR_output_epigenet_root.xlsx")




ref_gene <- "ref2"
ref_data <- qPCR_output_epigenet_root %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_epigenet_root %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  filter(genotype %in% c("col","nrpd2a","nrpd1a4","dcl3"))

# Graph
stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 
stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("col","nrpd2a","nrpd1a4", "rdr2","dcl2.6", "dcl3","dcl2/3","nrpe1")) # Choisir ordrer du facte_wrap



position=position_dodge(0.9)

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene !="ref1", genotype %in% c("col","nrpd2a","nrpd1a4","nrpe1","dcl2.6", "dcl3","dcl2/3", "rdr2")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph

# Save 
ggsave(filename="out/qPCR cDNA/RdDM_mutant.pdf", plot=my_graph, width = 5, height = 3)

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene =="AT4G14548", genotype %in% c("col","nrpd2a","nrpd1a4","dcl3")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph

# Save 
ggsave(filename="out/qPCR cDNA/RdDM_mutant_ARES.pdf", plot=my_graph, width = 3, height = 3)

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene !="ref1", genotype %in% c("col","nrpd2a","nrpd1a4","nrpe1","dcl2.6", "dcl3","dcl2/3", "rdr2")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=gene)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  facet_wrap(~genotype, nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/RdDM_mutant_gene.pdf", plot=my_graph, width = 8, height = 3)



my_graph <- 
stat_qPCR_cDNA %>% filter(gene !="ref1") %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



# RNAi control -----

qPCR_output_RNAi_control <- read_excel("data/qPCR_cDNA/qPCR_output_RNAi_control.xlsx")



ref_gene <- "ref2"
ref_data <- qPCR_output_RNAi_control %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_RNAi_control %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

ddCt_data_RNAi <- ddCt_data

stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 
position=position_dodge(0.9)

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene !="ref1") %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/RNAi_control.pdf", plot=my_graph, width = 5, height = 3)


my_graph <- 
  stat_qPCR_cDNA %>% filter(gene =="AT4G14548") %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/RNAi_control_solinc.pdf", plot=my_graph, width = 3, height = 3)



my_graph <- 
  stat_qPCR_cDNA %>% filter(gene %in%c("AT4G14540","AT4G14550"),
                            genotype%in%c("col","RNAi 1.2")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))
my_graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene %in%c("AT4G14540","AT4G14550"),
                            genotype%in%c("col","RNAi 1.2")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  ylim(-2,2)
my_graph



# Save 
ggsave(filename="out/qPCR cDNA/RNAi_control_nfyiaa14.pdf", plot=my_graph, width = 3, height = 3)



# Combine SALK exp1 and exp4 with the RNAi -------

 

ddCt_data_combine <- ddCt_data_14 %>%
  bind_rows(ddCt_data_RNAi)

stat_qPCR_cDNA <- 
  ddCt_data_combine %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#Anova 
position=position_dodge(0.9)

my_graph <- 
  stat_qPCR_cDNA %>% filter(gene !="ref1", genotype != "SALK140401") %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean, fill=genotype)) +
  geom_bar(position=position, stat="identity") +
  geom_errorbar(position=position, mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(position=position, aes(label=letter), size=4, vjust=-1, hjust=1)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  theme_bw()
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/RNAi and SALK control.pdf", plot=my_graph, width = 7, height = 3)






# kin NAA SALK14 and RNAi 1.2 -----
qPCR_output_NAA_kin <- read_excel("data/qPCR_cDNA/qPCR_output_NAA_kin.xlsx") %>%
  mutate(Cp=as.numeric(Cp))


#Data processing

ref_gene <- "ref2"
ref_data <- qPCR_output_NAA_kin %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_NAA_kin %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(time == 0, genotype =="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

#norm conc
ddCt_data_col <- ddCt_data %>% filter(genotype=="col")
ddCt_data_rnai <- ddCt_data %>% filter(genotype=="RNAi 1.2")
ddCt_data_salk <- ddCt_data %>% filter(genotype=="SALK_113294")

ddCt_data_all <- ddCt_data_col %>%
  bind_rows(ddCt_data_rnai, ddCt_data_salk)
#

stat_ddCt_data <- ddCt_data %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            (median=median(-ddCt)),
            ecart_type=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

my_graph <- 
stat_ddCt_data %>% filter(! gene %in% c("PID","ref1")) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype,nrow=1)
my_graph

# Save 
ggsave(filename="out/qPCR cDNA/NAA kin_1.pdf", plot=my_graph, width = 5, height = 3)


# NAA kin RNAi only ------

qPCR_output_1_SOLINC_NAA <- read_excel("data/kin_NAA/qPCR_output_1_SOLINC_NAA.xlsx") %>%
  dplyr::select(-Pos)

#Data processing

ref_gene <- "ref2"
ref_data <- qPCR_output_1_SOLINC_NAA %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_1_SOLINC_NAA %>%
  filter(gene != ref_gene) 


dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(time == 0, genotype =="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

#norm conc
ddCt_data_col <- ddCt_data %>% filter(genotype=="col")
ddCt_data_rnai12 <- ddCt_data %>% filter(genotype=="RNAi12")
ddCt_data_rnai32 <- ddCt_data %>% filter(genotype=="RNAi32")
ddCt_data_rnai62 <- ddCt_data %>% filter(genotype=="RNAi62")

ddCt_data_all <- ddCt_data_col %>%
  bind_rows(ddCt_data_rnai12, ddCt_data_rnai32,ddCt_data_rnai62)
#

stat_ddCt_data <- ddCt_data %>%
  group_by(gene,time,genotype) %>%
  summarise(mean=mean(-ddCt), 
            (median=median(-ddCt)),
            ecart_type=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart

my_graph <- 
stat_ddCt_data%>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype,nrow=1)+
  ylim(-2.4,4.8)
my_graph

stat_ddCt_data_corr <- stat_ddCt_data %>%
  filter(genotype =="col") %>%
  bind_rows(stat_ddCt_data %>%
              filter(gene !="AT4G14548"))
  
my_graph <- stat_ddCt_data_corr %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype,nrow=1)+
  scale_color_viridis_d()
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/NAA kin_RNAi_withSOLINC.pdf", plot=my_graph, width = 7, height = 3)



my_graph <- stat_ddCt_data_corr %>%
  filter(gene %in% c("AT4G14540","AT4G14550"),
         genotype%in% c("col","RNAi12"))%>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=0,size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype,nrow=1)+
  scale_color_viridis_d()
my_graph


# Save 
ggsave(filename="out/qPCR cDNA/NAA kin_RNAi_withSOLINC_nfyIAA14.pdf", plot=my_graph, width = 7, height = 3)






