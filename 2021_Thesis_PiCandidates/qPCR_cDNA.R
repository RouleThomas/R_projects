
# Out directory if it does not exist
dir.create("out/qPCR_cDNA/", showWarnings = FALSE, recursive = TRUE)


# Pi kin_exp1 ------
#kin

cinetique_Pi_qPCR <- read_excel("data/qPCR/cinetique Pi_qPCR.xlsx") %>%
  filter(genotype =="col") %>%
  dplyr::select(-Pos,-genotype)
  
cinetique_Pi_qPCR <- read_excel("data/qPCR/cinetique Pi_qPCR.xlsx",sheet=2) %>%
  filter(genotype =="col") %>%
  dplyr::select(-Pos,-genotype) %>%
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
         time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, time, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	


 ddCt_data %>% inner_join(AT1G55525)

my_graph <- 
stat_Pi %>%	inner_join(AT1G55525) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/AT1G55525_Pi.pdf", plot=my_graph, width = 5, height = 2)




my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002093) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002093_Pi.pdf", plot=my_graph, width = 5, height = 2)




my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002900) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002900_Pi.pdf", plot=my_graph, width = 5, height = 2)


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_005603) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_005603_Pi.pdf", plot=my_graph, width = 5, height = 2)


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_005603) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_005603_Pi.pdf", plot=my_graph, width = 5, height = 2)



# Pearson corelation -----
tidy_pearson <- 
ddCt_data %>%
  filter(!is.na(ddCt)) %>%
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
  gather(gene, pearson_corr, -gene1) %>% 
  filter(gene1 %in% row.names(nc_gene_count), 
         gene %in% row.names(coding_gene_count)) )



# exp _2 Pi starvation kin ------



cinetique_Pi_qPCR <- read_excel("data/qPCR/cinetique Pi_qPCR_2.xlsx",sheet=1) %>%
  filter(genotype =="col") %>%
  dplyr::select(-genotype) %>%
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
         time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat_Pi <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, time, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	


Pi_gene <- c("SPX3","IPS1","CLE14") %>% as.tibble() %>% rename(gene=value)

my_graph <- 
  stat_Pi %>%	inner_join(Pi_gene) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean)) +
  geom_line(aes(color= condition), size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/genes_Pi.pdf", plot=my_graph, width = 5, height = 2)



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



my_graph <- 
  stat_Pi %>%	inner_join(AT1G55525) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/AT1G55525_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002093 <- XLOC_002093 %>% rbind("AT1G08925")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002093) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002093_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002900 <- XLOC_002900 %>% rbind("AT1G08687")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002900) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002900_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)

XLOC_005603 <- XLOC_005603 %>% rbind("AT3G03315")

my_graph <- 
  stat_Pi %>%	inner_join(XLOC_005603) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph




# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_005603_Pi_adapt.pdf", plot=my_graph, width = 5, height = 2)



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

my_graph <- 
  stat_Pi %>%	inner_join(AT1G55525) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/AT1G55525_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002093 <- XLOC_002093 %>% rbind("AT1G08925")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002093) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002093_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002900 <- XLOC_002900 %>% rbind("AT1G08687")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002900) %>% unique() %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002900_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)

XLOC_005603 <- XLOC_005603 %>% rbind("AT3G03315")

my_graph <- 
  stat_Pi %>%	inner_join(XLOC_005603) %>% unique() %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph




# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_005603_salt_adapt.pdf", plot=my_graph, width = 5, height = 2)


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

my_graph <- 
  stat_Pi %>%	inner_join(AT1G55525) %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/AT1G55525_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002093 <- XLOC_002093 %>% rbind("AT1G08925")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002093) %>% unique() %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002093_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)


XLOC_002900 <- XLOC_002900 %>% rbind("AT1G08687")


my_graph <- 
  stat_Pi %>%	inner_join(XLOC_002900) %>% unique() %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_002900_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)

XLOC_005603 <- XLOC_005603 %>% rbind("AT3G03315")

my_graph <- 
  stat_Pi %>%	inner_join(XLOC_005603) %>% unique() %>%
  ggplot(data=.,	 mapping=aes(x=condition,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= condition)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)
my_graph




# Save 
ggsave(filename="out/qPCR_cDNA/XLOC_005603_ABA_adapt.pdf", plot=my_graph, width = 5, height = 2)

























