# Out directory if it does not exist
dir.create("out/RNAseq_graph/", showWarnings = FALSE, recursive = TRUE)



library(tidyverse)
library(corrplot)
library(readr)
library(readxl)
theme_set(theme_bw())
position=position_dodge(.10)






normcounts <- read.table("./in/normcounts.txt") %>%
  rownames_to_column("ID")



#Pour graphique RNAseq
normcounts_tidy <- normcounts %>% 
  gather(col_t0_A:ler_t2_C, key = "condition", value = "counts", na.rm = TRUE) %>%
  separate(condition, c("ecotype", "time", "replicate"), sep="_")
stat_normcounts_tidy <- normcounts_tidy %>%
  group_by(ID, ecotype, time) %>%
  summarise(mean_counts=mean(counts), 
            median_counts=median(counts),
            ecart_type=sd(counts),
            n=n(), 
            erreur_std=ecart_type/sqrt(n))


#AT1G55525 -----
gene_list <- c("AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220")


 
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/AT1G55525_trans.pdf", plot=my_graph, width = 7, height = 5)

gene_list <- c("AT1G55525")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/AT1G55525.pdf", plot=my_graph, width = 7, height = 5)

#AT3G61198 ----------



gene_list <- c("AT4G08920", "AT5G03150", "AT1G72150", "AT1G22530")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/AT3G61198_trans.pdf", plot=my_graph, width = 7, height = 5)


gene_list <- c("AT3G61198")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/AT3G61198.pdf", plot=my_graph, width = 7, height = 5)



#AT4G13495------

gene_list <- c("AT4G13495")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph
# Save 
ggsave(filename="out/RNAseq_graph/AT4G13495.pdf", plot=my_graph, width = 7, height = 5)


gene_list <- c("AT1G28560")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph
# Save 
ggsave(filename="out/RNAseq_graph/AT4G13495_trans.pdf", plot=my_graph, width = 7, height = 5)


#XLOC_002093/AT1G08925 -----

gene_list <- c("AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

my_graph
# Save 
ggsave(filename="out/RNAseq_graph/XLOC_002093_trans.pdf", plot=my_graph, width = 7, height = 5)

gene_list <- c("XLOC_002093")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

my_graph
# Save 
ggsave(filename="out/RNAseq_graph/XLOC_002093.pdf", plot=my_graph, width = 7, height = 5)


#XLOC_002755/AT1G08173

gene_list <- c("AT2G23430", "AT1G72160", "AT5G58010", "AT1G67710", "AT3G04630", "AT1G13260", "AT2G34680", "AT1G72150", "AT3G22400", "AT4G33880", "AT4G34580", "AT2G31090")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph
# Save 
ggsave(filename="out/RNAseq_graph/XLOC_002755_trans.pdf", plot=my_graph, width = 7, height = 5)

gene_list <- c("XLOC_002755")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph
# Save 
ggsave(filename="out/RNAseq_graph/XLOC_002755.pdf", plot=my_graph, width = 7, height = 5)


gene_list <- c("XLOC_002755")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
my_graph
# Save 
ggsave(filename="out/RNAseq_graph/XLOC_002755.pdf", plot=my_graph, width = 7, height = 5)


# UCR ------------



gene_list <- c("Col_NEW_RNA_F_29664","AT5G21940")
gene_list <- c("AT2G18328","AT2G18323")
gene_list <- c("AT2G05170","AT2G05171")
gene_list <- c("AT4G22550","AT4G22545")



normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 
my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")+
  theme_bw()
my_graph








