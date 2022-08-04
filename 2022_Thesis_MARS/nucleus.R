# Theme loading
source("lib/mytheme.R")

# Out directory if it does not exist
dir.create("out/nucleus/", showWarnings = FALSE, recursive = TRUE)

#Data import
nucleus_qPCR_output <- read_excel("data/nucleus/nucleus_qPCR_output.xlsx")%>%
  mutate(Cp=as.numeric(Cp))

#Data processing
ref_gene <- "ref2"
ref_data <- nucleus_qPCR_output %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- nucleus_qPCR_output %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(condition == "total") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

# NEW UPDATED with isofmrs -----

ddCt_data$gene <- factor(ddCt_data$gene, levels=c("ref1","APOLO","ASCO","MARS.isoform1","MARS.isoform2","MARS.isoform3", "MHAL","U6"))

my_graph <- 
ddCt_data %>% 
  mutate(ddCt=-ddCt) %>%
  filter(condition=="nucleus")%>%
  ggbarplot(., x = "gene", y = "ddCt",add = "mean_se") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "APOLO"), c("ref1", "ASCO"), c("ref1", "MARS.isoform1"), 
                                        c("ref1", "MARS.isoform2"), c("ref1", "MARS.isoform3"), c("ref1", "U6") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() + 
  theme(legend.position="none")


my_graph

# Save 
ggsave(filename="out/nucleus/nucleus_isoforms.pdf", plot=my_graph, width = 7, height = 4) 




# t.test

ddCt_data <- ddCt_data %>% 
  filter(condition == "nucleus", gene %in% c("U6","MHAL", "ASCO", "APOLO", "ref1")) %>% 
  dplyr::select(gene, ddCt) %>% 
  mutate(log2FC=-ddCt) %>% 
  left_join(selected_gene_nucleus)

ddCt_data$gene_name <- factor(ddCt_data$gene_name, levels=c("ref1","APOLO","ASCO","MHAL.5'","U6"))

my_graph <- 
ddCt_data %>% 
ggbarplot(., x = "gene_name", y = "log2FC",add = "mean_se", fill="gene_name") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "APOLO"), c("ref1", "ASCO"), c("ref1", "MHAL.5'"), c("ref1", "U6") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() +
  scale_fill_manual(breaks=selected_gene_nucleus$gene,
                     values=selected_gene_nucleus$gene_color,
                     labels=selected_gene_nucleus$gene_name,) + 
  theme(legend.position="none")

my_graph

# Save 
ggsave(filename="out/nucleus/nucleus.pdf", plot=my_graph, width = 7, height = 4) 


ddCt_data %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="gene") +
  theme_bw() 


## nucleus iso MARS ---------

nucleurs_iso_qPCR <- read_excel("data/nucleus/nucleurs_iso_qPCR.xlsx")%>%
  mutate(Cp=as.numeric(Cp)) %>%
  select(-Pos)


ref_gene <- "ref"
ref_data <- nucleurs_iso_qPCR %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)

gene_data <- nucleurs_iso_qPCR %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(condition == "total") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 




stat_iso <- ddCt_data %>%		
  filter(replicate %in% c(1,2,5))%>%
  group_by(gene, condition) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


my_graph <-  
  stat_iso %>% 
  filter(condition =="nucleus")%>%
  ggplot(data=.,	 mapping=aes(x=gene,	 y=mean,	 group=gene)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=gene)) +		
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw()

my_graph

#######################can play with the ref or pooler with the previous manip plutot!
# Save 
ggsave(filename="out/XXX.pdf", plot=my_graph, width = 7, height = 2)





