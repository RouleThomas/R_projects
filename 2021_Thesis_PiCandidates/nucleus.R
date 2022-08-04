




# Out directory if it does not exist
dir.create("out/nucleus/", showWarnings = FALSE, recursive = TRUE)

#Data import
nucleus_qPCR_output <- read_excel("data/nucleus/nucleus_qPCR.xlsx")


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
  filter(condition == "Total") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 


# t.test

ddCt_data <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "AT1G55525","XLOC_002093",
                                             "XLOC_002900","XLOC_005603", "XLOC_002755")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data$gene <- factor(ddCt_data$gene, levels=c("ref1", "XLOC_002755","XLOC_005603","XLOC_002900",
                                                       "XLOC_002093","AT1G55525","U6"))

my_graph <- 
ddCt_data %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                   comparisons = list(c("ref1", "XLOC_002755"), c("ref1", "XLOC_005603"), c("ref1", "XLOC_002900"), c("ref1", "XLOC_002093"), 
                                      c("ref1", "AT1G55525", "ref1","U6")), 
                   label = "p.label", bracket.size= 0.5)
  
my_graph


# Save 
ggsave(filename="out/nucleus/nucleus.pdf", plot=my_graph, width = 7, height = 4) 

#CAL -----

ddCt_data1 <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "CAL")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data1$gene <- factor(ddCt_data1$gene, levels=c("ref1", "CAL","U6"))


my_graph <- 
  ddCt_data1 %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "CAL"), c("ref1", "U6")), 
                     label = "p.label", bracket.size= 0.5)

my_graph


# Save 
ggsave(filename="out/nucleus/nucleus_CAL.pdf", plot=my_graph, width = 7, height = 4) 

#LCA ------


ddCt_data1 <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "LCA")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data1$gene <- factor(ddCt_data1$gene, levels=c("ref1", "LCA","U6"))


my_graph <- 
  ddCt_data1 %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "LCA"), c("ref1", "U6")), 
                     label = "p.label", bracket.size= 0.5)

my_graph


# Save 
ggsave(filename="out/nucleus/nucleus_LCA.pdf", plot=my_graph, width = 7, height = 4) 

#BAP ------


ddCt_data1 <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "BAP")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data1$gene <- factor(ddCt_data1$gene, levels=c("ref1", "BAP","U6"))


my_graph <- 
  ddCt_data1 %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "BAP"), c("ref1", "U6")), 
                     label = "p.label", bracket.size= 0.5)

my_graph


# Save 
ggsave(filename="out/nucleus/nucleus_BAP.pdf", plot=my_graph, width = 7, height = 4) 

#BAP ------


ddCt_data1 <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "NPC48")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data1$gene <- factor(ddCt_data1$gene, levels=c("ref1", "NPC48","U6"))


my_graph <- 
  ddCt_data1 %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "NPC48"), c("ref1", "U6")), 
                     label = "p.label", bracket.size= 0.5)

my_graph


# Save 
ggsave(filename="out/nucleus/nucleus_NPC48.pdf", plot=my_graph, width = 7, height = 4) 

#IAA ------


ddCt_data1 <- ddCt_data %>% 
  filter(condition == "Noyau", gene %in% c("U6","ref1", "IAA")) %>% 
  mutate(log2FC=-ddCt) %>% 
  dplyr::select(gene, log2FC)

ddCt_data1$gene <- factor(ddCt_data1$gene, levels=c("ref1", "IAA","U6"))


my_graph <- 
  ddCt_data1 %>% 
  ggbarplot(., x = "gene", y = "log2FC",add = "mean_se", fill="grey") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("ref1", "IAA"), c("ref1", "U6")), 
                     label = "p.label", bracket.size= 0.5)

my_graph


# Save 
ggsave(filename="out/nucleus/nucleus_IAA.pdf", plot=my_graph, width = 7, height = 4) 









