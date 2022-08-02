# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/ChIRP/", showWarnings = FALSE, recursive = TRUE)


# ChIRP RNA ----
# Data import
qPCR_output_all_ChIRPRNA <- read_excel("data/ChIRP/qPCR_output_all_ChIRPRNA.xlsx")

#Data processing

qPCR_chirp_odd <- qPCR_output_all_ChIRPRNA %>% filter(condition == "odd")
qPCR_chirp_even <- qPCR_output_all_ChIRPRNA %>% filter(condition == "even")
qPCR_chirp_lacZ <- qPCR_output_all_ChIRPRNA %>% filter(condition == "lacZ")

qPCR_chirp_input <- qPCR_output_all_ChIRPRNA %>%
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(10))

qPCR_IP_input_chirp <- qPCR_chirp_odd %>%
  bind_rows(qPCR_chirp_even, qPCR_chirp_lacZ, qPCR_chirp_input)


qPCR_input_chirp <- qPCR_IP_input_chirp %>%
  filter(condition == "input") %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_IP_chirp <- qPCR_IP_input_chirp %>%
  filter(condition %in% c("odd","even","lacZ")) %>%
  mutate(Cp_sample=Cp) %>%
  dplyr::select(-Cp)


qPCR_percent_chirp <- qPCR_IP_chirp %>%
  full_join(qPCR_input_chirp) %>%
  mutate(Percent=2^-(Cp_sample-Cp_input)*100)


ggplot(data=qPCR_percent_chirp, aes(x=condition, y=Percent)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free")




qPCR_percent_chirp_IP <- qPCR_percent_chirp %>% filter(condition %in% c("odd","even")) %>%
  dplyr::select(-Cp_sample,-Cp_input) %>%
  rename(Percent_IP=Percent)
qPCR_percent_chirp_lacZ <- qPCR_percent_chirp %>% filter(condition %in% c("lacZ")) %>%
  dplyr::select(-Cp_sample,-Cp_input, -condition)

qPCR_percent_chirp_IP_lacZ <- qPCR_percent_chirp_IP %>%
  left_join(qPCR_percent_chirp_lacZ) %>%
  mutate(Enrichment=Percent_IP/Percent)

qPCR_percent_chirp_IP_lacZ %>%
  ggplot(data=., aes(x=condition, y=Enrichment)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free")


qPCR_percent_chirp_IP_lacZ %>% filter(gene %in% c("ref2", "MARS.3'")) %>%
  ggplot(data=., aes(x=condition, y=Enrichment)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~gene, nrow=1) +
  theme_bw()



# Statistics and Graph
# Anova
qPCR_percent_chirp_IP_lacZ$gene <- factor(qPCR_percent_chirp_IP_lacZ$gene, levels=c("ref2","MARS.3'"))


my_graph <- 
  qPCR_percent_chirp_IP_lacZ %>% filter(gene %in% c("ref2", "MARS.3'")) %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref2","MARS.3'") ), 
                     label = "p.label", bracket.size= 0.5) +
  scale_fill_grey() +
  theme_bw() +
  ylab("Enrichent") +
  xlab("") + 
  theme(legend.position="none")

my_graph


# Save 
ggsave(filename="out/ChIRP/ChIRPRNA.pdf", plot=my_graph, width = 5, height = 4) 


# ChIRPDNA ------
#Data import
qPCR_output_all_ChIRPDNA <- read_excel("data/ChIRP/qPCR_output_all_ChIRPDNA.xlsx")

#Data processing

qPCR_chirp_odd <- qPCR_output_all_ChIRPDNA %>% filter(condition == "ODD")
qPCR_chirp_even <- qPCR_output_all_ChIRPDNA %>% filter(condition == "EVEN")
qPCR_chirp_lacZ <- qPCR_output_all_ChIRPDNA %>% filter(condition == "LacZ")

qPCR_chirp_input <- qPCR_output_all_ChIRPDNA %>%
  filter(condition == "input") %>%
  mutate(Cp=Cp-log2(0.5))

qPCR_IP_input_chirp <- qPCR_chirp_odd %>%
  bind_rows(qPCR_chirp_even, qPCR_chirp_lacZ, qPCR_chirp_input)


qPCR_input_chirp <- qPCR_IP_input_chirp %>%
  filter(condition == "input") %>%
  mutate(Cp_input=Cp) %>%
  dplyr::select(-Cp, -condition)

qPCR_IP_chirp <- qPCR_IP_input_chirp %>%
  filter(condition %in% c("ODD","EVEN","LacZ")) %>%
  mutate(Cp_sample=Cp) %>%
  dplyr::select(-Cp)


qPCR_percent_chirp <- qPCR_IP_chirp %>%
  full_join(qPCR_input_chirp) %>%
  mutate(Percent=2^-(Cp_sample-Cp_input)*100)


stat_CHIRP <- qPCR_percent_chirp %>%
  group_by(condition, gene) %>%
  summarise(mean=mean(Percent),
            median= median(Percent),
            SD=sd(Percent), #ecart-type
            n=n(), #nombre d'?chantillon
            erreur_std=SD/sqrt(n))

position = position_dodge(0.8)

stat_CHIRP %>%
  ggplot(data=., mapping=aes(x=condition, y=mean, group=gene)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5),position=position)+
  xlab(label = "")+
  ylab(label = "")+
  facet_wrap(~gene, scale="free")+
  theme_bw()


qPCR_percent_chirp_IP <- qPCR_percent_chirp %>% filter(condition %in% c("ODD","EVEN")) %>%
  dplyr::select(-Cp_sample,-Cp_input) %>%
  rename(Percent_IP=Percent)
qPCR_percent_chirp_lacZ <- qPCR_percent_chirp %>% filter(condition %in% c("LacZ")) %>%
  dplyr::select(-Cp_sample,-Cp_input, -condition)

qPCR_percent_chirp_IP_lacZ <- qPCR_percent_chirp_IP %>%
  left_join(qPCR_percent_chirp_lacZ) %>%
  mutate(Enrichment=Percent_IP/Percent)

qPCR_percent_chirp_IP_lacZ %>%
  ggplot(data=., aes(x=condition, y=Enrichment)) +
  geom_boxplot() +
  facet_wrap(~gene, scale="free")


stat_CHIRP <- qPCR_percent_chirp_IP_lacZ %>%
  group_by(condition, gene) %>%
  summarise(mean=mean(Enrichment),
            median= median(Enrichment),
            SD=sd(Enrichment), #ecart-type
            n=n(), #nombre d'?chantillon
            erreur_std=SD/sqrt(n))

stat_CHIRP %>% filter(gene %in% c("intergenic1","intergenic2","MRN1.5'","MRN1.3'", "MRN1.promoter", "H3.3")) %>%
  ggplot(data=., mapping=aes(x=condition, y=mean, group=gene)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5),position=position)+
  xlab(label = "")+
  ylab(label = "")+
  facet_wrap(~gene, scale="free")+
  geom_hline(yintercept=1,linetype=2)+
  theme_bw()


# Anova
qPCR_percent_chirp_IP_lacZ$gene <- factor(qPCR_percent_chirp_IP_lacZ$gene, levels=c("H3.3", "intergenic1", "intergenic2","MRN1.promoter","MRN1.5'","MRN1.3'"))


my_graph <- 
  qPCR_percent_chirp_IP_lacZ %>% filter(gene %in% c("intergenic1","intergenic2","MRN1.5'","MRN1.3'", "MRN1.promoter", "H3.3")) %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("H3.3","intergenic1"), c("H3.3","intergenic2"), c("H3.3","MRN1.promoter"), c("H3.3","MRN1.5'"), c("H3.3","MRN1.3'") ), 
                     label = "p.label", bracket.size= 0.5) +
  scale_fill_grey() +
  theme_bw() +
  ylab("Enrichent") +
  xlab("") +
  geom_hline(yintercept=0.42,linetype=2)+
  theme(legend.position="none")

my_graph


# Save 
ggsave(filename="out/ChIRP/ChIRPDNA.pdf", plot=my_graph, width = 5, height = 4) 


my_graph <- 
  qPCR_percent_chirp_IP_lacZ %>% filter(gene %in% c("intergenic1","intergenic2","MRN1.5'","MRN1.3'", "MRN1.promoter", "H3.3")) %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("H3.3","intergenic1"), c("H3.3","intergenic2"), c("H3.3","MRN1.promoter"), c("H3.3","MRN1.5'"), c("H3.3","MRN1.3'") ), 
                     label = "p.label", bracket.size= 0.5) +
  scale_fill_grey() +
  theme_bw() +
  ylab("Enrichent") +
  xlab("") +
  geom_hline(yintercept=0.42,linetype=2)+
  theme(legend.position="none")

my_graph


# Save 
ggsave(filename="out/ChIRP/ChIRPDNA.pdf", plot=my_graph, width = 5, height = 4) 




my_graph <- 
  qPCR_percent_chirp_IP_lacZ %>% filter(gene %in% c("intergenic1","intergenic2","MRN1.5'","MRN1.3'", "MRN1.promoter")) %>%
  ggbarplot(., x = "gene", y = "Enrichment",add = "mean_se", fill="gene") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("intergenic1","intergenic2"), c("intergenic1","MRN1.promoter"), c("intergenic1","MRN1.5'"), c("intergenic1","MRN1.3'") ), 
                     label = "p.label", bracket.size= 0.5) +
  scale_fill_grey() +
  theme_bw() +
  ylab("Enrichment") +
  xlab("") +
  geom_hline(yintercept=0.42,linetype=2)+
  theme(legend.position="none")

my_graph


# Save 
ggsave(filename="out/ChIRP/ChIRPDNA_withoutH3.pdf", plot=my_graph, width = 5, height = 4) 



