# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/qPCR_cDNA_2/", showWarnings = FALSE, recursive = TRUE)


position=position_dodge(0.9)


# MRN1 kin exp1 -----

col_35SMRN1_and_mrn1_kinetic <- read_excel("data/qPCR_cDNA/col 35SMRN1 and mrn1 kinetic.xlsx")

ref_gene <- "ref2"
ref_data <- col_35SMRN1_and_mrn1_kinetic %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- col_35SMRN1_and_mrn1_kinetic %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_mrn <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			



my_graph <- 
stat_mrn %>% filter(gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1)

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/MRN1 mutants _ aba kin.pdf", plot=my_graph, width = 8, height = 3.5)


my_graph <- 
stat_mrn %>% filter(gene %in% c("RAB18","RD29B")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )  + theme(legend.position = "none")


my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/MRN1 mutants _ aba kin_aba genes.pdf", plot=my_graph, width = 8, height = 3.5)


# ANOVA ------

anova_MRN1 <- ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")
anova_MHAL <- ddCt_data %>% 
  filter(gene == "MHAL") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="MHAL") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")
anova_CYP705A12 <- ddCt_data %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")
anova_CYP71A16 <- ddCt_data %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")
anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  mutate(genotype = c("35smrn1","no","mrn1")) %>%
  dplyr::select(-"comp") %>%
  filter(genotype != "no")

anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_MHAL, anova_CYP705A12,anova_CYP71A16,anova_RAB18,anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble() 

anova_all_gene_35smrn1 <- anova_all_gene %>%
  filter(genotype=="35smrn1") %>%
  mutate(diff=-diff)
anova_all_gene_mrn1 <- anova_all_gene %>%
  filter(genotype=="mrn1") 

anova_all_gene <- anova_all_gene_35smrn1 %>%
  bind_rows(anova_all_gene_mrn1)

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene




my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.0001))) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_mrn_cluster.pdf", plot=my_graph, width = 3.5, height = 3.5)



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_mrn_aba.pdf", plot=my_graph, width = 3.5, height = 3.5)



  


# mro kin ------
#kin
qPCR_cDNA_mro_kin <- read_excel("data/qPCR_cDNA/qPCR_cDNA_mro_kin.xlsx")



ref_gene <- "ref"
ref_data <- qPCR_cDNA_mro_kin %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_cDNA_mro_kin %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_mro <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			





my_graph <- 
  stat_mro %>% filter(gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) + theme(legend.position = "none")+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/mro mutants _ aba kin.pdf", plot=my_graph, width = 5, height = 3.5)


my_graph <- 
  stat_mro %>% filter(gene %in% c("RAB18","RD29B")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) + theme(legend.position = "none")+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/mro mutants _ aba kin_aba genes.pdf", plot=my_graph, width = 5, height = 3.5)


# ANOVA ------

anova_MRN1 <- ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_MHAL <- ddCt_data %>% 
  filter(gene == "MHAL") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="MHAL") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP705A12 <- ddCt_data %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP71A16 <- ddCt_data %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_MHAL, anova_CYP705A12,anova_CYP71A16,anova_RAB18,anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble() 

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene


my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.0001))) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_mro_cluster.pdf", plot=my_graph, width = 3.5, height = 3.5)



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff - 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_mro_aba.pdf", plot=my_graph, width = 3.5, height = 3.5)





# water kin ------
#kin

qPCR_output_water <- read_excel("data/qPCR_cDNA/qPCR_water_control.xlsx", sheet=1) 


ref_gene <- "ref"
ref_data <- qPCR_output_water %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_output_water %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_water <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

my_graph <- 
stat_water %>% filter(gene %in% c("RAB18","RD29B")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) + theme(legend.position = "none")+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/water_aba resp.pdf", plot=my_graph, width = 9, height = 3.5)



my_graph <- 
stat_water %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MHAL")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) + theme(legend.position = "none")+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/water_cluster resp.pdf", plot=my_graph, width = 9, height = 3.5)


#ANOVA -----

anova_MRN1 <- ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_MHAL <- ddCt_data %>% 
  filter(gene == "MHAL") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MHAL") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP705A12 <- ddCt_data %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP71A16 <- ddCt_data %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_MHAL, anova_CYP705A12,anova_CYP71A16,anova_RAB18,anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble()

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene


anova_all_gene$genotype <- factor(anova_all_gene$genotype, c("rnai92","rnai38","rnai18","pro"))


my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.0001))) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_water_cluster.pdf", plot=my_graph, width = 3.5, height = 3.5)



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff - 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_water_aba.pdf", plot=my_graph, width = 3.5, height = 3.5)




# High Treatment ABA 100uM -------

qPCR_output_100uM <- read_excel("data/qPCR_cDNA/qPCR_High_ABA_100uM.xlsx",sheet=1) %>%
  mutate(Cp=as.numeric(Cp))

ref_gene <- "ref"
ref_data <- qPCR_output_100uM %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_output_100uM %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

#not normalized


stat_High <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

my_graph <- 
stat_High %>% filter(gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +		
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1 )+
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/high_cluster resp.pdf", plot=my_graph, width = 9, height = 3.5)


my_graph <- 
  stat_High %>% filter(gene %in% c("RAB18","RD29B")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1 )+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/high_cluster resp_aba resp.pdf", plot=my_graph, width = 9, height = 3.5)


# ANOVA ------

anova_MRN1 <- ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_MHAL <- ddCt_data %>% 
  filter(gene == "MHAL") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MHAL") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP705A12 <- ddCt_data %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP71A16 <- ddCt_data %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_MHAL, anova_CYP705A12,anova_CYP71A16,anova_RAB18,anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble() 

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene


anova_all_gene$genotype <- factor(anova_all_gene$genotype, c("rnai92","rnai38","rnai18","pro"))


my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.0001))) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_high_cluster.pdf", plot=my_graph, width = 3.5, height = 3.5)



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff - 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_high_aba.pdf", plot=my_graph, width = 3.5, height = 3.5)








# RNAi1.8, line 3 kinetic ------

qPCR_cDNA_RNAi_line3 <- read_excel("data/qPCR_cDNA/qPCR_cDNA_RNAi_line3.xlsx")


ref_gene <- "ref2"
ref_data <- qPCR_cDNA_RNAi_line3 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_cDNA_RNAi_line3 %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

stat_rnai <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1","MHAL")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/kin RNAi 3 marneral genes.pdf", plot=my_graph, width = 5, height = 3.5)


my_graph <- 
  stat_rnai %>%	filter(gene%in% c("RAB18","RD29B")) %>% 
  left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/kin RNAi 3_ABA genes.pdf", plot=my_graph, width = 5, height = 3.5)


#re all qPCR output ------
re_all_qPCR_output <- read_excel("data/qPCR_cDNA/re_all_qPCR_output.xlsx",sheet=1)


ref_gene <- "ref"
ref_data <- re_all_qPCR_output %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- re_all_qPCR_output %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


ddCt_data_rnai <- ddCt_data

col_ABA <- ddCt_data %>%
  filter(genotype=="col") %>%
  add_column(condition="ABA stimulus") 

stat_rnai <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

#ANOVA------

ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  summary


anova_MRN1 <- ddCt_data %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_MHAL <- ddCt_data %>% 
  filter(gene == "MHAL") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="MHAL") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP705A12 <- ddCt_data %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP71A16 <- ddCt_data %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col"))  %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  filter(rowname %in% c("pro-col","rnai18-col","rnai38-col","rnai92-col")) %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_MHAL, anova_CYP705A12,anova_CYP71A16,anova_RAB18,anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble()

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene


anova_all_gene$genotype <- factor(anova_all_gene$genotype, c("rnai92","rnai38","rnai18","pro"))


my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff + sign(diff + 0.0001)), angle=45) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_cluster.pdf", plot=my_graph, width = 3.5, height = 3.5)



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=diff - 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )+
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_all_aba.pdf", plot=my_graph, width = 3.5, height = 3.5)


###############




my_graph <- 
  stat_rnai %>%	filter(gene%in% c("MHAL")) %>% left_join(genotype_metadata)%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/MARS_only_RNAi_kin.pdf", plot=my_graph, width = 9, height = 4)






my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1","MHAL")) %>% left_join(genotype_metadata)%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral.pdf", plot=my_graph, width = 9, height = 4)





my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1"),
                       genotype%in%c("col","rnai92")) %>% left_join(genotype_metadata)%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_withoutMARS.pdf", plot=my_graph, width = 9, height = 3)


my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1"),
                       genotype%in%c("col","rnai92"),time%in%c(0,0.25,0.5,1)) %>% left_join(genotype_metadata)%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.04)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_withoutMARS_zoom.pdf", plot=my_graph, width = 6, height = 3)





#  norm t0 each genotype ##### ------

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_col <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(genotype=="col")

dCt_data_control <- dCt_data %>%
  filter(genotype == "pro") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_pro <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(genotype=="pro")

dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai18") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_rnai18 <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(genotype=="rnai18")

dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai38") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_rnai38 <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(genotype=="rnai38")

dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai92") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_rnai92 <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(genotype=="rnai92")

dCt_data <- ddCt_data_rnai92 %>%
  bind_rows(ddCt_data_rnai38) %>%
  bind_rows(ddCt_data_rnai18) %>%
  bind_rows(ddCt_data_col) %>%
  bind_rows(ddCt_data_pro)


######

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_rnai <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	


# Save 

my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1")) %>% left_join(genotype_metadata)%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_norm.pdf", plot=my_graph, width = 9, height = 4)

my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1")) %>% left_join(genotype_metadata) %>% filter(time<1.5) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_norm_zoom.pdf", plot=my_graph, width = 9, height = 4)


my_graph <- 
  stat_rnai %>%	filter(gene%in% c("RAB18","RD29B")) %>% 
  left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_ABA.pdf", plot=my_graph, width = 9, height = 3)

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_ABA_norm.pdf", plot=my_graph, width = 5, height = 3.5)


#time 0 ----

stat_qPCR_cDNA <- ddCt_data %>%
  filter(time =="0") %>%
  group_by(gene, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(time =="0") %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/marneral cluster control condition_gene.pdf", plot=my_graph, width = 7, height = 3.5)




my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(!gene %in% c("RAB18","RD29B"), genotype_name != "SALK_133089") %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/marneral cluster control condition_gene_RNAi_only.pdf", plot=my_graph, width = 7, height = 3.5)




my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(!gene %in% c("RAB18","RD29B"),
                                                             genotype%in%c("rnai92")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=gene)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~genotype, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )+
  ylim(-10,5)
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/marneral cluster control condition_gene RNAi 1.pdf", plot=my_graph, width = 4, height = 3)






#All stress Col-0 ------


qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=1)


##Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_nitrate <- ddCt_data %>%
  add_column(condition="NO3- starvation")

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time %in% c("0","12","24")) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral nitrogen starvation.pdf", plot=my_graph, width = 5, height = 3)


# Pi starvation -----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=2) %>%
  clean_genotype_code

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_Pi <- ddCt_data %>%
  add_column(condition="Pi starvation") %>%
  filter(time %in% c("0","2","1", "12","24"))

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time %in% c("0","2","1", "12","24")) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral Pi starvation.pdf", plot=my_graph, width = 5, height = 3)


# NAA treatment ----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=3)

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_naa <- ddCt_data %>%
  add_column(condition="NAA stimulus") %>%
  filter(time != 24)

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time != 24) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral NAA treatment.pdf", plot=my_graph, width = 5, height = 3)

#Heat Stress -----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=4)

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) 
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_hs <- ddCt_data %>%
  add_column(condition="Heat stress")

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral Heat Stress treatment.pdf", plot=my_graph, width = 5, height = 3)


Col_all_stress <- col_nitrate %>% bind_rows(col_ABA, col_Pi,col_naa, col_hs) %>%
  dplyr::select(time,gene,ddCt,condition)


stat_qPCR_cDNA <- Col_all_stress %>%
  group_by(gene, time,condition) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA$condition <- factor(stat_qPCR_cDNA$condition, c("Pi starvation","NO3- starvation","NAA stimulus", "Heat stress", "ABA stimulus"))


my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) +
  facet_wrap(~condition, nrow=1, scales = "free_x") +
  xlab("time")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/Col-0 all treatment.pdf", plot=my_graph, width = 10, height = 3)






















# high res ------

high_res_kin <- read_excel("data/qPCR_cDNA/high_res_kin.xlsx",sheet=2)


ref_gene <- "ref1"
ref_data <- high_res_kin %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- high_res_kin %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_high <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 		

 
stat_high %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MHAL"), genotype %in% c("col","pro")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=genotype)) +
  geom_line(aes(color = genotype),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1) 


my_graph <- 
  stat_high %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MHAL")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 
my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/rnai_cluster resp.pdf", plot=my_graph, width = 7, height = 3.5)




my_graph <- 
  stat_high %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/rnai_aba resp.pdf", plot=my_graph, width = 7, height = 3.5)


# pro kin ------ de la merde, car marneral cluster trop haut!! utiliser pro exp high res
montrer expe tdna intergenic!! ou a lancer!!!
NON: mettre le rnai.18 avec le salk time 2!!!!!!!!!!!!, imposible car salk time 2 jai pas 
passé les RAB18 et RD29B!!





output_pro_ABA_1 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_1.xlsx")
plate_plan_pro_ABA_1 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_1.xlsx")
output_pro_ABA_2 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_2.xlsx")
plate_plan_pro_ABA_2 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_2.xlsx")
output_pro_ABA_3 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_3.xlsx")
plate_plan_pro_ABA_3 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_3.xlsx")

output_pro_ABA_4 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_4.xlsx")
plate_plan_pro_ABA_4_Zenhua1 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_4_Zenhua1.xlsx")

output_pro_ABA_5 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_5.xlsx")
plate_plan_pro_ABA_5_Zenhua2 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_5_Zenhua2.xlsx")

output_pro_ABA_6 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_6.xlsx")
plate_plan_pro_ABA_6_Zenhua3 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_6_Zenhua3.xlsx")

output_pro_ABA_7 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_7.xlsx")
plate_plan_pro_ABA_7 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_7.xlsx")

output_pro_ABA_8 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_8.xlsx")
plate_plan_pro_ABA_8 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_8.xlsx")

output_pro_ABA_9 <- read_excel("data/qPCR_cDNA/in/output_pro_ABA_9.xlsx")
plate_plan_pro_ABA_9 <- read_excel("data/qPCR_cDNA/in/plate_plan_pro_ABA_9.xlsx")



qPCR_proABA_1 <- output_pro_ABA_1 %>%
  left_join(plate_plan_pro_ABA_1) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp)) %>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)


qPCR_proABA_2 <- output_pro_ABA_2 %>%
  left_join(plate_plan_pro_ABA_2) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_proABA_3 <- output_pro_ABA_3 %>%
  left_join(plate_plan_pro_ABA_3) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_Zenhua_1 <- output_pro_ABA_4 %>%
  left_join(plate_plan_pro_ABA_4_Zenhua1) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_Zenhua_2 <- output_pro_ABA_5 %>%
  left_join(plate_plan_pro_ABA_5_Zenhua2) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_Zenhua_3 <- output_pro_ABA_6 %>%
  left_join(plate_plan_pro_ABA_6_Zenhua3) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_proABA_7 <- output_pro_ABA_7 %>%
  left_join(plate_plan_pro_ABA_7) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)


qPCR_proABA_8 <- output_pro_ABA_8 %>%
  left_join(plate_plan_pro_ABA_8) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)

qPCR_proABA_9 <- output_pro_ABA_9 %>%
  left_join(plate_plan_pro_ABA_9) %>%
  separate(condition, into = c("genotype", "gene"), sep = "/") %>%
  separate(genotype, c("genotype", "condition"), sep="-") %>%
  separate(genotype, c("genotype", "time"), sep="_") %>%
  separate(genotype, c("genotype", "replicate"), sep="&") %>%
  filter(!is.na(Cp))%>%
  dplyr::select(genotype, time, replicate, condition, gene, Cp)



qPCR_proABA <- qPCR_proABA_1 %>% bind_rows(qPCR_proABA_2, qPCR_proABA_3) 
qPCR_proABA <- qPCR_proABA_1 %>% bind_rows(qPCR_Zenhua_1, qPCR_Zenhua_2, qPCR_Zenhua_3) 
qPCR_proABA <- qPCR_proABA_1 %>% bind_rows(qPCR_proABA_9) 
qPCR_proABA <- qPCR_proABA_1 %>% bind_rows(qPCR_proABA_2, qPCR_proABA_3, qPCR_proABA_9) 

#Choisir son gene de rÃ©fÃ©rence (choisir celui qui est le moins dispersÃ©; le plus linÃ©aire, ici ref2)
ggplot(data=qPCR_proABA, aes(x=time, y=Cp)) + geom_jitter(aes(color=condition), width = 0.1) + facet_wrap(~gene, scale = "free")



#Selection du gene de ref & retirer le gene de ref du dataframe (=gen_data ici)

#Cp mean ref
ref_gene <- "ref2"
ref_data <- qPCR_proABA %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_proABA %>%
  filter(gene != ref_gene) 

# Fusion dataframe (sans le gene de reference) avec Cp du gene de reference pour chaque echantillon (=ref_data)
#& ajout colonne dCp (=variation en nb de cycle/gene par rapport au gene de reference)
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



# Ajout dCp/gene par rapport Ã  la condition controle (ici Col NoABA) = dCt_data_control
dCt_data_control <- dCt_data %>%
  filter(time == 0, genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))




# Calcul du ddCt = Ct gN par rapport Ã  Ct gN ref couplÃ© Ã  Ct condition de ref + calcul du Fold Change (FC)
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

ddCt_data_time0 <- ddCt_data %>%
  filter(time == 0) %>%
  dplyr::select(-condition)

ddCt_data_times <- ddCt_data %>%
  filter(time != 0) %>%
  filter(condition =="ABA")%>%
  dplyr::select(-condition)

ddCt_data <- ddCt_data_time0 %>%
  bind_rows(ddCt_data_times)


ddCt_data_to_change_time <- read_excel("data/qPCR_cDNA/in/ddCt_data_to change time.xlsx",sheet=2) %>%
  mutate(ddCt=as.numeric(ddCt),
         time=as.numeric(time))
  



stat_pro <- ddCt_data_to_change_time %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	


my_graph <- 
  stat_pro %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MHAL"), genotype != "nrpd2a") %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 
my_graph



my_graph <- 
  stat_pro %>% filter(gene %in% c("RD29B","RAB18","ABI5"), genotype != "nrpd2a") %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 
my_graph


# ying powder ------


qPCR_output_ying <- read_excel("data/qPCR_cDNA/MARS CRISPR/qPCR_output_ying.xlsx") %>%
  select(-Well) %>%
  drop_na()



#Cp mean ref
ref_gene <- "ref2"
ref_data <- qPCR_output_ying %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- qPCR_output_ying %>%
  filter(gene != ref_gene) 

# Fusion dataframe (sans le gene de reference) avec Cp du gene de reference pour chaque echantillon (=ref_data)
#& ajout colonne dCp (=variation en nb de cycle/gene par rapport au gene de reference)
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



# Ajout dCp/gene par rapport Ã  la condition controle (ici Col NoABA) = dCt_data_control
dCt_data_control <- dCt_data %>%
  filter(time == "t0", genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))




# Calcul du ddCt = Ct gN par rapport Ã  Ct gN ref couplÃ© Ã  Ct condition de ref + calcul du Fold Change (FC)
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))


stat_ying <- ddCt_data %>%
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	


stat_ying %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +	
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 


ddCt_data %>%
  mutate(replicate=as.character(replicate))%>%
  ggplot(., aes(time, -ddCt))+
  geom_jitter(position = position_jitter(height = .2, width = .2), aes(color=replicate))+
  facet_grid(genotype~gene)


# 





