# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/meDIP/", showWarnings = FALSE, recursive = TRUE)



meDIP_qPCR_output <- read_excel("data/meDIP/meDIP_qPCR_output.xlsx") 


meDIP_output_qPCR_input <- meDIP_qPCR_output %>%
  filter(condition =="input") %>%
  mutate(Cp_input=Cp) %>%
  mutate(Cp_input=Cp_input-log2(10)) %>%
  dplyr::select(-condition, -Cp)


meDIP_output_qPCR_sample <- meDIP_qPCR_output %>%
  filter(condition != "input") %>%
  left_join(meDIP_output_qPCR_input) %>%
  mutate(Percent=2^-(Cp-Cp_input)*100)




stat_meDIP <- meDIP_output_qPCR_sample %>%
  group_by(gene, condition, genotype) %>%
  summarise(mean=mean(Percent), 
            median= median(Percent),
            SD=sd(Percent), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_meDIP <- 
  meDIP_output_qPCR_sample %>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "Percent", "condition")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_meDIP$genotype <- factor(stat_meDIP$genotype, levels=c("col","rnai92","rnai38"))
stat_meDIP$gene <- factor(stat_meDIP$gene, levels=c("APOLO.5'","MARS.5'","MARS.3'"))


#ANOVa GOOD ------
label_y_shift <-
  max(stat_meDIP$mean + stat_meDIP$sem) * 0.05


my_graph <- 
  stat_meDIP %>% filter(genotype %in% c("col","rnai38","rnai92")) %>%
  ggplot(data=., mapping=aes(x=gene, y=mean,fill=condition)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  geom_errorbar(position=position_dodge(width=0.9), mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y =0.025, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "% input")+
  facet_wrap(~genotype, nrow=1, scale="free") +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1))) +
  scale_fill_grey(start=0.8, end=0.2)


my_graph


# Save 
ggsave(filename="out/meDIP/meDIP_APOLOMARS_input.pdf", plot=my_graph, width = 5, height = 4) 



#Enrichment ----


meDIP_output_qPCR_sample_IGG <- meDIP_output_qPCR_sample %>%
  filter(condition =="IGG") %>%
  mutate(Percent_IGG=Percent) %>%
  dplyr::select(-Percent, -condition, -Cp, -Cp_input) 

meDIP_output_qPCR_sample_mC <- meDIP_output_qPCR_sample %>%
  filter(condition =="mC") %>%
  mutate(Percent_mC=Percent) %>%
  dplyr::select(-Percent, -condition, -Cp, -Cp_input) %>%
  left_join(meDIP_output_qPCR_sample_IGG) %>%
  mutate(Enrichment=Percent_mC/Percent_IGG)




stat_Enrichment <- meDIP_output_qPCR_sample_mC %>%
  group_by(gene, genotype) %>%
  summarise(mean=mean(Enrichment), 
            median= median(Enrichment),
            SD=sd(Enrichment), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_Enrichment <- 
  meDIP_output_qPCR_sample_mC %>%
  group_by(genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "Enrichment", "gene")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

stat_Enrichment$genotype <- factor(stat_Enrichment$genotype, levels=c("col","rnai92","rnai38"))
stat_Enrichment$gene <- factor(stat_Enrichment$gene, levels=c("APOLO.5'","MARS.5'","MARS.3'"))


my_graph <- 
  stat_Enrichment %>% filter(genotype %in% c("col","rnai38","rnai92","rnai18")) %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=gene, y=mean)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=gene)) +
  geom_errorbar(position=position_dodge(width=0.9), mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y =2, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "Enrichment")+
  facet_wrap(~genotype_name, nrow=1, scale="free") +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  theme(strip.text = element_text(size = rel(1))) +
  scale_fill_manual(values=selected_gene_meDIP$gene_color,)

my_graph

# Save 
ggsave(filename="out/meDIP/meDIP_APOLOMARS_Enrichment.pdf", plot=my_graph, width = 5, height = 3) 




