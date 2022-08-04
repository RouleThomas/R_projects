
# Out directory if it does not exist
dir.create("out/qPCR_cDNA/", showWarnings = FALSE, recursive = TRUE)


# Pi kin_exp1 ------
#kin

cinetique_Pi_qPCR <- read_excel("data/qPCR/cinetique Pi_qPCR.xlsx",sheet=2) %>%
  dplyr::select(-Pos) %>%
  filter(!is.na(Cp))
  

ref_gene <- "ref2"
ref_data <- cinetique_Pi_qPCR %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- cinetique_Pi_qPCR %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "HighP",
         time == "0",
         genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, time, condition,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

gene_list <- c("CLE14","CAL") %>%
  as.tibble() %>% rename(gene=value)


my_graph <- 
stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(genotype~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/CAL_Pi.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("MRN1","MHAL","CYP705A12","CYP71A16") %>%
  as.tibble() %>% rename(gene=value)


my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(genotype~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/MARS_Pi.pdf", plot=my_graph, width = 5, height = 2)




# Pearson corelation -----



tidy_pearson <- 
  ddCt_data %>%
  filter(!is.na(ddCt), genotype =="ler") %>%
  dplyr::select(replicate, time, condition, gene, ddCt) %>%
  add_column(test="A") %>%
  unite(test,replicate, time,condition) %>%
  spread(key=test , value=ddCt) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("gene") 



View(cor(t(tidy_pearson), t(tidy_pearson), method="pearson",use="complete.obs") %>% 
       as.data.frame  %>% # Matrice de correlation
       rownames_to_column("gene1") %>% 
       gather(gene, pearson_corr, -gene1)  )



# exp _2 Pi starvation kin ------



cinetique_Pi_qPCR <- read_excel("data/qPCR/cinetique Pi_qPCR_2.xlsx",sheet=1) %>%
  filter(!is.na(Cp)) %>%
  mutate(Cp=as.numeric(Cp))


ref_gene <- "ref2"
ref_data <- cinetique_Pi_qPCR %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- cinetique_Pi_qPCR %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "HighP",
         time == "0",
         genotype =="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, time, condition,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	





gene_list <- c("CLE14","CAL") %>%
  as.tibble() %>% rename(gene=value)


my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(genotype~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/CAL_Pi.pdf", plot=my_graph, width = 5, height = 2)











# PAN ------
all_stress_PAN <- read_excel("data/qPCR/all_stress_PAN.xlsx") %>%
  dplyr::select(-Pos) %>%
  filter(!is.na(Cp))

#Pi ----
Pi_PAN <- all_stress_PAN %>%
  filter(manip =="PhosphateAdaptation") %>%
  dplyr::select(-manip)




ref_gene <- "ref2"
ref_data <- Pi_PAN %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- Pi_PAN %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "HighP") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

gene_list <- c("BAP","ncBAP") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/BAP_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("IAA","ncIAA") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/IAA_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)





#Salt ----
Salt_PAN <- all_stress_PAN %>%
  filter(manip =="SaltStress") %>%
  dplyr::select(-manip)




ref_gene <- "ref2"
ref_data <- Salt_PAN %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- Salt_PAN %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "NoNaCl") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoNaCl","NaCl"))


gene_list <- c("BAP","ncBAP") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/BAP_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("IAA","ncIAA") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/IAA_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("BAP","ncBAP") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/IAA_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


#ABA ----
ABA_PAN <- all_stress_PAN %>%
  filter(manip =="ABAStress") %>%
  dplyr::select(-manip)




ref_gene <- "ref2"
ref_data <- ABA_PAN %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- ABA_PAN %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "NoABA") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoABA","ABA"))

gene_list <- c("BAP","ncBAP") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/BAP_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("IAA","ncIAA") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/IAA_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


# candidats et stress -------
#ABA----
candidat_et_stress_folder <- read_excel("data/qPCR/candidat et stress folder.xlsx",sheet=1)


ref_gene <- "ref2"
ref_data <- candidat_et_stress_folder %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- candidat_et_stress_folder %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "NoABA",genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoABA","ABA"))

gene_list <- c("MRN1","MARS","CYP71A16","CYP705A12") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/MARS_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("CAL","CLE14","AT1G63240") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/CAL_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("PAL","PGX1","AT3G26616","AT3G26618") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))
stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoABA","ABA"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))





my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/PAL_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("AT5G38000","LCA","AT5G38010","AT5G38020") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))
stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoABA","ABA"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))


my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

write.csv(stat_Pi %>% inner_join(gene_list), file="ABA_lca.csv")

# Save 
ggsave(filename="out/qPCR_cDNA/LCA_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


#Pi----
candidat_et_stress_folder <- read_excel("data/qPCR/candidat et stress folder.xlsx",sheet=2)


ref_gene <- "ref2"
ref_data <- candidat_et_stress_folder %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- candidat_et_stress_folder %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "HighP",genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_Pi$condition <- factor(stat_Pi$condition, levels = c("HighP","LowP"))

gene_list <- c("MRN1","MARS","CYP71A16","CYP705A12") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/MARS_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("CAL","CLE14","AT1G63240") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/CAL_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("PAL","PGX1","AT3G26616","AT3G26618") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))





my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/PAL_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("AT5G38000","LCA","AT5G38010","AT5G38020") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))


my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

write.csv(  stat_Pi %>% inner_join(gene_list), file="Pi_adapt_lca.csv")

# Save 
ggsave(filename="out/qPCR_cDNA/LCA_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)


#Nacl----
candidat_et_stress_folder <- read_excel("data/qPCR/candidat et stress folder.xlsx",sheet=3)


ref_gene <- "ref1"
ref_data <- candidat_et_stress_folder %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- candidat_et_stress_folder %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(condition == "NoNACL",genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_Pi$condition <- factor(stat_Pi$condition, levels = c("NoNACL","NACL"))

gene_list <- c("MRN1","MARS","CYP71A16","CYP705A12") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/MARS_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)

gene_list <- c("CAL","CLE14","AT1G63240") %>%
  as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/CAL_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("PAL","PGX1","AT3G26616","AT3G26618") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("PGX1","PAL","AT3G26616","AT3G26618"))





my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/PAL_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


gene_list <- c("AT5G38000","LCA","AT5G38010","AT5G38020") %>%
  as.tibble() %>% rename(gene=value) 

gene_list$gene <- factor(gene_list$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))
stat_Pi$gene <- factor(stat_Pi$gene, levels = c("AT5G38000","LCA","AT5G38010","AT5G38020"))


my_graph <- 
  stat_Pi %>% inner_join(gene_list) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean,fill= genotype)) +
  geom_bar(stat='identity', position="dodge") +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2, position=position_dodge(0.9))+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_grid(~gene)
my_graph

write.csv( stat_Pi %>% inner_join(gene_list) , file="salt_lca.csv")

# Save 
ggsave(filename="out/qPCR_cDNA/LCA_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)
















