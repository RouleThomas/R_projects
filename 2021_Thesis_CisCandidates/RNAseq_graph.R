# Out directory if it does not exist
dir.create("out/RNAseq_graph/", showWarnings = FALSE, recursive = TRUE)


library(RColorBrewer)
library(tidyverse)
library(corrplot)
library(readr)
library(readxl)
theme_set(theme_bw())
position=position_dodge(.10)








normcounts <- read.table("data/rnaseq/normcounts.txt") %>%
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


#CLE14 -----
gene_list <- c("AT1G63230","AT1G63240","Col_NEW_RNA_R_6757","AT1G63245", "AT1G63250")


 
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

normcounts_current$ID <- factor(normcounts_current$ID, levels = c("AT1G63230","AT1G63240","Col_NEW_RNA_R_6757","AT1G63245", "AT1G63250"))
stat_normcounts_current$ID <- factor(stat_normcounts_current$ID, levels = c("AT1G63230","AT1G63240","Col_NEW_RNA_R_6757","AT1G63245", "AT1G63250"))


my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free", nrow=1)

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/CAL.pdf", plot=my_graph, width = 7, height = 5)



#PGX1 -----
gene_list <- c("AT3G26600","AT3G26610","AT3G26612", "AT3G26618", "AT3G26620")



normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

normcounts_current$ID <- factor(normcounts_current$ID, levels = c("AT3G26600","AT3G26610","AT3G26612", "AT3G26618", "AT3G26620"))


my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free",nrow=1)

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/PGX.pdf", plot=my_graph, width = 7, height = 5)


#BAP_AT3G61190 -----
gene_list <- c("AT3G61180","AT3G61190","AT3G61198", "AT3G61200","AT3G61210")



normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

normcounts_current$ID <- factor(normcounts_current$ID, levels = c("AT3G61180","AT3G61190","AT3G61198", "AT3G61200","AT3G61210"))


my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free",nrow=1)

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/BAP.pdf", plot=my_graph, width = 7, height = 5)


#IAA -----
gene_list <- c("AT4G14530","AT4G14540","AT4G06195","AT4G14548","AT4G14550","AT4G06200","AT4G14560","AT4G14570")



normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

normcounts_current$ID <- factor(normcounts_current$ID, levels = c("AT4G14530","AT4G14540","AT4G06195","AT4G14548","AT4G14550","AT4G06200","AT4G14560","AT4G14570"))


my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free",nrow=1)

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/IAA.pdf", plot=my_graph, width = 7, height = 5)




#LCA -----
gene_list <- c("AT5G38000","AT5G38005","AT5G38010","AT5G38020")



normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

normcounts_current$ID <- factor(normcounts_current$ID, levels = c("AT5G38000","AT5G38005","AT5G38010","AT5G38020"))


my_graph <- 
  ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free",nrow=1)

my_graph 
# Save 
ggsave(filename="out/RNAseq_graph/LCA.pdf", plot=my_graph, width = 7, height = 5)


#Pearson Pi RNAseq -------


normcounts_col <- normcounts %>% 
  select(ID, starts_with("c"))


#lincRNA
AT5G38005 <- normcounts_col %>%
  filter(ID=="AT5G38005") %>% 
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") 
  
  



#coding gene avec GO

all_gene_count <-
  normcounts_col %>%
  filter(ID !="AT5G38005") %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID")



pearson_correlation_df <- 
  cor(t(AT5G38005), t(all_gene_count), method="pearson") %>% 
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>%
  filter(gene %in% c("AT5G37830" ,"AT5G37840","AT5G37850","AT5G37860","AT5G05555","AT5G37870","AT5G37872","AT5G37875","AT5G37880","AT5G37890","AT5G37900","AT5G37910","AT5G05565",
                    "AT5G37920","AT5G37930","AT5G37940","
                     AT5G05585","AT5G37950","AT5G37960","AT5G37970","AT5G37980","AT5G37990","AT5G05605","
                     AT5G38000","AT5G38005","AT5G38010","AT5G38020","AT5G38030","AT5G38035","AT5G38037","AT5G38040",
                     "AT5G38050","AT5G38060","AT5G38070","AT5G38080","AT5G38090","AT5G38096","AT5G38100"))
#   14 - 13      
         
         










#pearson Araport11 -------
Araport11_all_counts <- 
  read_tsv("data/pearson/Araport11_all_gene_counts.tsv.gz") %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "normalised_count")
#CLE14-----
cluster <- c("AT1G63230","AT1G63240","AT1G09745","AT1G63245", "AT1G63250")


cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/CLE14_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


#PGX1-----
cluster <- c("AT3G26600","AT3G26610","AT3G26612", "AT3G26618", "AT3G26620")


cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/PGX1_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()




#BAP-----
cluster <- c("AT3G61180","AT3G61190","AT3G61198", "AT3G61200","AT3G61210")


cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/BAP_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()

#IAA-----
cluster <- c("AT4G14530","AT4G14540","AT4G06195","AT4G14548","AT4G14550","AT4G06200","AT4G14560","AT4G14570")
cluster <- c("AT4G14530","AT4G14540","AT4G14548","AT4G14550","AT4G14560","AT4G14570")

cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/IAA_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


#LCA-----
cluster <- c("AT5G38000","AT5G38005","AT5G38010","AT5G38020")


cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/LCA_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


#LCA ------

cluster <- c("AT5G37830" ,"AT5G37840","AT5G37850","AT5G37860","AT5G05555","AT5G37870","AT5G37872","AT5G37875","AT5G37880","AT5G37890","AT5G37900","AT5G37910","AT5G05565",
             "AT5G37920","AT5G37930","AT5G37940","
                     AT5G05585","AT5G37950","AT5G37960","AT5G37970","AT5G37980","AT5G37990","AT5G05605","
                     AT5G38000","AT5G38005","AT5G38010","AT5G38020","AT5G38030","AT5G38035","AT5G38037","AT5G38040",
             "AT5G38050","AT5G38060","AT5G38070","AT5G38080","AT5G38090","AT5G38096","AT5G38100")



cluster_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/LCA_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()




# TRAVA correlation ------
TRAVA_all_counts <- 
  read_tsv("data/pearson/TraVA_normcounts.tsv.gz") %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "normalised_count")


#IAA-----
cluster <- c("AT4G14530","AT4G14540","AT4G14548","AT4G14550","AT4G14560","AT4G14570")


cluster_expr <-
  TRAVA_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/IAA_TRAVA.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


#other candidiq-----
cluster <- c("AT1G48610", "AT1G48620","AT1G48625","AT1G48630" ,"AT1G48635", "AT1G48640")


cluster_expr <-
  TRAVA_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% cluster) %>%
  mutate(gene_id = factor(gene_id, levels = cluster))

cluster_corr <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

cluster_sig <-
  cluster_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/Pearson/RACK1B_TRAVA.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()
