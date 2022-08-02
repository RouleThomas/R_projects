library(tidyverse)
library(readxl)
library(grid)



# Out directory if it does not exist
dir.create("out/qPCR_cDNA/", showWarnings = FALSE, recursive = TRUE)


#env genet ------

env_genet_lca <- read_excel("data/qPCR_cDNA/env_genet_lca.xlsx") %>%
  dplyr::select(-Pos)

ref_gene <- "ref2"
ref_data <- env_genet_lca %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- env_genet_lca %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

selected_genes_LCA <-
  tibble(gene=c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040"),
         gene_color=c("#fdbb84","#fdbb84","#fdbb84","#fdbb84", "#e34a33", "#fdbb84",
                      "#fdbb84","#fdbb84","#fdbb84"))

my_graph <- 
stat %>% left_join(selected_genes_LCA) %>% filter(gene !="ref1",genotype!="collca") %>%
  ggplot(data=.,	 mapping=aes(x=gene,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= gene)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~genotype,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    values=selected_genes_LCA$gene_color)
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/insertional_mutant_control.pdf", plot=my_graph, width = 8, height = 3)


my_graph <- 
  stat %>% left_join(selected_genes_LCA) %>% filter(gene !="ref1",genotype!="collca") %>%
  ggplot(data=.,	 mapping=aes(x=genotype,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= genotype)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/insertional_mutant_control_all_genes.pdf", plot=my_graph, width = 8, height = 3)

my_graph <- 
  stat %>% left_join(selected_genes_LCA) %>% filter(gene !="ref1",genotype%in%c("col","SAIL_1215_A04")) %>%
  ggplot(data=.,	 mapping=aes(x=genotype,	 y=mean)) +
  geom_bar(stat='identity', aes(fill= genotype)) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/insertional_mutant_control_genes.pdf", plot=my_graph, width = 8, height = 3)


#IAA -----
lca_aia <- read_excel("data/qPCR_cDNA/lca_aia.xlsx")%>%
  dplyr::select(-Pos) 

ref_gene <- "ref1"
ref_data <- lca_aia %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- lca_aia %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         condition =="NoAIA") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene,genotype,condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

selected_genes_LCA <-
  tibble(gene=c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040"),
         gene_color=c("#fdbb84","#fdbb84","#fdbb84","#fdbb84", "#e34a33", "#fdbb84",
                      "#fdbb84","#fdbb84","#fdbb84"))

position=position_dodge(0.9)

stat$condition <- factor(stat$condition, levels = c("NoAIA", "AIA"))

my_graph <- 
  stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")) %>% left_join(selected_genes_LCA) %>%
  ggplot(data=.,	 mapping=aes(x=genotype, y=mean, fill=condition)) +
  geom_bar(stat='identity',position=position) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2,position=position)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph

write.csv(stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")), 
          "aia.csv")

# Save 
ggsave(filename="out/qPCR_cDNA/insertional_mutant_aia.pdf", plot=my_graph, width = 8, height = 3)



#Nadapt ------


lca_nadapt <- read_excel("data/qPCR_cDNA/lca_nadapt.xlsx")%>%
  dplyr::select(-Pos) 


ref_gene <- "ref2"
ref_data <- lca_nadapt %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- lca_nadapt %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col",
         condition =="HighN") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene,genotype,condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

selected_genes_LCA <-
  tibble(gene=c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040"),
         gene_color=c("#fdbb84","#fdbb84","#fdbb84","#fdbb84", "#e34a33", "#fdbb84",
                      "#fdbb84","#fdbb84","#fdbb84"))

position=position_dodge(0.9)


my_graph <- 
  stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")) %>% left_join(selected_genes_LCA) %>%
  ggplot(data=.,	 mapping=aes(x=genotype, y=mean, fill=condition)) +
  geom_bar(stat='identity',position=position) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2,position=position)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph

write.csv(stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")), 
          "nadapt.csv")


# Save 
ggsave(filename="out/qPCR_cDNA/insertional_mutant_Nadapt.pdf", plot=my_graph, width = 8, height = 3)


#ABA kin with LCA -----
ABA_kin_lca <- read_excel("data/qPCR_cDNA/qPCR_ABA_lca.xlsx")



ref_gene <- "ref2"
ref_data <- ABA_kin_lca %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- ABA_kin_lca %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(genotype =="col",
  time==0,
  condition =="NoABA") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene, condition,time,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	



my_graph <- 
  stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +			
  facet_grid(condition~genotype)
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/ABA_kinetic.pdf", plot=my_graph, width = 7, height = 5)


#ANOVA kin aba ----



anova_AT5G37970 <- ddCt_data %>% 
  filter(gene == "AT5G37970") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G37970") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G37980 <- ddCt_data %>% 
  filter(gene == "AT5G37980") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G37980") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G37990 <- ddCt_data %>% 
  filter(gene == "AT5G37990") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G37990") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38000 <- ddCt_data %>% 
  filter(gene == "AT5G38000") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38000") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38005 <- ddCt_data %>% 
  filter(gene == "AT5G38005") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38005") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38010 <- ddCt_data %>% 
  filter(gene == "AT5G38010") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38010") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38020 <- ddCt_data %>% 
  filter(gene == "AT5G38020") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38020") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38030 <- ddCt_data %>% 
  filter(gene == "AT5G38030") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38030") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_AT5G38040 <- ddCt_data %>% 
  filter(gene == "AT5G38040") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="AT5G38040") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_AT5G38040 %>%
  bind_rows(anova_AT5G38030,anova_AT5G38020,anova_AT5G38010,anova_AT5G38005,anova_AT5G38000,anova_AT5G37990,anova_AT5G37980,anova_AT5G37970) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble()




anova_all_gene$genotype <- factor(anova_all_gene$genotype, c("rnai92","rnai38","rnai18","pro"))


my_graph <- 
  anova_all_gene %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.1)), angle=45) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/anova_all_ABA_kin_lca.pdf", plot=my_graph, width = 3.5, height = 3.5)


#mutant voie ABA lca -----
mutant_voie_ABA_lca <- read_excel("data/qPCR_cDNA/mutant_voie_ABA_lca.xlsx") %>%
  filter(condition =="NoABA") %>%
  dplyr::select(-Pos,-condition)

ref_gene <- "ref2"
ref_data <- mutant_voie_ABA_lca %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- mutant_voie_ABA_lca %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref) %>%
  dplyr::select(-Cp)

dCt_data_control <- dCt_data %>%
  filter(genotype=="col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 



stat <- ddCt_data %>%
  filter(!is.na(ddCt)) %>%
  group_by(gene,genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

selected_genes_LCA <-
  tibble(gene=c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040"),
         gene_color=c("#fdbb84","#fdbb84","#fdbb84","#fdbb84", "#e34a33", "#fdbb84",
                      "#fdbb84","#fdbb84","#fdbb84"))

position=position_dodge(0.9)


my_graph <- 
  stat %>% filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")) %>% left_join(selected_genes_LCA) %>%
  ggplot(data=.,	 mapping=aes(x=genotype, y=mean)) +
  geom_bar(stat='identity',position=position) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2,position=position)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/epigenetic_control.pdf", plot=my_graph, width = 8, height = 3)


# LCA genes under all the stress ------

LCA_genes_all_conditions <- read_excel("data/qPCR_cDNA/LCA_genes_all_conditions.xlsx") %>%
  filter(genotype=="col", gene %in% c("AT5G38000","AT5G38005","AT5G38010"), condition %in% c("ABA","AIA","LowN","NACL","LowP")) 

#LCA_genes_all_conditions$condition <- factor(LCA_genes_all_conditions$condition, levels = c("NoABA", "AIA"))



my_graph <- 
LCA_genes_all_conditions  %>%
  ggplot(data=.,	 mapping=aes(x=condition, y=mean)) +
  geom_bar(stat='identity',position=position) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2,position=position)+
  geom_hline(yintercept=0) +			
  geom_hline(yintercept=1,linetype=2)+
  geom_hline(yintercept=-1,linetype=2)+
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/LCA_all_stress.pdf", plot=my_graph, width = 8, height = 3)

#Root methylome FPKM ------

lca_root_tissue_FPKM <- 
  read_excel("data/qPCR_cDNA/lca_root_tissue_FPKM.xlsx") %>%
  as.data.frame() %>% 
  pivot_longer(c(`epidermis`,`cortex`, `endodermis1` , `endodermis2`, `stele` , `root cap`))

lca_root_tissue_FPKM$name <- factor(lca_root_tissue_FPKM$name, levels = c("epidermis","cortex","endodermis1","endodermis2","stele","root cap"))

#pearson methylome-----

lca_root_tissue_FPKM_pearson <- read_excel("data/qPCR_cDNA/lca_root_tissue_FPKM.xlsx") 



#lincRNA
AT5G38005 <- lca_root_tissue_FPKM_pearson %>%
  filter(gene=="AT5G38005") %>% 
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("gene") 





#coding genes

all_gene_count <-
  lca_root_tissue_FPKM_pearson %>%
  filter(gene !="AT5G38005") %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("gene")


pearson_correlation_df <- 
  cor(t(AT5G38005), t(all_gene_count), method="pearson") %>% 
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>%
  dplyr::select(-ncRNA) %>%
  as_tibble()

my_graph <- 
  lca_root_tissue_FPKM  %>% 
  filter(gene %in% c("AT5G37970","AT5G37980","AT5G37990","AT5G38000","AT5G38005","AT5G38010", "AT5G38020","AT5G38030","AT5G38040")) %>%
  left_join(pearson_correlation_df) %>%
  ggplot(data=.,	 mapping=aes(x=name, y=value, label=pearson_corr%>%formatC(digits = 2))) +
  geom_bar(stat='identity',position=position) +		
  geom_hline(yintercept=0) +
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  geom_text(x=4, y=27,color="red")+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/root_methylome_lca.pdf", plot=my_graph, width = 8, height = 3)



solinc_root_tissue_FPKM <- 
  read_excel("data/qPCR_cDNA/solinc_root_tissue_FPKM.xlsx") %>%
  as.data.frame() %>% 
  pivot_longer(c(`epidermis`,`cortex`, `endodermis1` , `endodermis2`, `stele` , `root cap`))

solinc_root_tissue_FPKM$name <- factor(solinc_root_tissue_FPKM$name, levels = c("epidermis","cortex","endodermis1","endodermis2","stele","root cap"))


#lincRNA
AT4G14548 <- read_excel("data/qPCR_cDNA/solinc_root_tissue_FPKM.xlsx") %>%
  filter(gene=="AT4G14548") %>% 
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("gene")



#coding genes

all_gene_count <-
  read_excel("data/qPCR_cDNA/solinc_root_tissue_FPKM.xlsx")  %>%
  filter(gene !="AT4G14548") %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("gene")


pearson_correlation_df <- 
  cor(t(AT4G14548), t(all_gene_count), method="pearson") %>% 
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>%
  dplyr::select(-ncRNA) %>%
  as_tibble()

my_graph <- 
  solinc_root_tissue_FPKM  %>% 
  filter(gene %in% c(c("AT4G14530","AT4G14540","AT4G06195","AT4G14548","AT4G14550","AT4G06200","AT4G14560","AT4G14570"))) %>%
  left_join(pearson_correlation_df) %>%
  ggplot(data=.,	 mapping=aes(x=name, y=value, label=pearson_corr%>%formatC(digits = 2))) +
  geom_bar(stat='identity',position=position) +		
  geom_hline(yintercept=0) +
  theme_bw() +	
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  facet_wrap(~gene,	 nrow=1)+
  geom_text(x=4, y=152,color="red")+
  theme(strip.text = element_text(size = rel(1)))
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/root_methylome_solinc.pdf", plot=my_graph, width = 8, height = 3)

















