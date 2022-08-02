library('tidyverse')
library('stringr')
library('readxl')
library('readr')
library('Hmisc')


# mhal aba GTF ------

#normcounts <- read_excel("P:/marneral/RNAseq_mhal_all/Analysis_with_GTF/Analysis_with_GTF/normcounts.xlsx")

normcounts <- read_csv("data/RNAseq_mhal/BH0.01/out/normcounts.csv")
normcounts <- read_delim("data/RNAseq_mhal/BH0.01_MARS corrected/normcounts.csv", 
                         ";", escape_double = FALSE, trim_ws = TRUE)  


#normcounts <- read_excel("data/RNAseq_mhal/normcounts.xlsx")


#clustering ------

clustering_DE <- read_excel("data/RNAseq_mhal/clustering_DE.xlsx")
#clustering_DE <- read_excel("data/RNAseq_mhal/clustering_DE.xlsx",sheet=2)


diffcounts <- normcounts %>% inner_join(clustering_DE) %>% 
  mutate(col0=rowMeans(diffcounts[,c(2:4)]))%>%
  mutate(col4=rowMeans(diffcounts[,c(5:7)]))%>%
  mutate(RNAi0=rowMeans(diffcounts[,c(8:10)]))%>%
  mutate(RNAi4=rowMeans(diffcounts[,c(11:13)]))%>%
  dplyr::select(col0,RNAi0,col4,RNAi4, ID)%>%
  filter(!duplicated(ID))
  
  
diffcounts <- as.data.frame(diffcounts) 

row.names(diffcounts)=diffcounts$ID
diffcounts$ID=NULL


library(pheatmap)

heatmap=pheatmap(diffcounts,scale = "row",cluster_cols = F )

saveRDS(heatmap, "heatmap.rds")

cluster_gene=cbind(diffcounts,
                   cluster = cutree(heatmap$tree_row,
                                    k = 4))%>%add_rownames(var="ID")%>%arrange(desc(cluster))

scaled_cluster_gene=cbind(diffcounts %>% as.matrix() %>% t %>% scale() %>% t,
                          cluster = cutree(heatmap$tree_row, k = 4))%>%
              as.data.frame %>% add_rownames(var="ID")%>%arrange(desc(cluster))

# collect genes from the cluster
cluster_genes_mhal <- cluster_gene %>%
  rename("code"="ID") %>%
  mutate(code=as.numeric(code))%>%
  left_join(clustering_DE)%>%
  dplyr::select(ID,cluster)%>%
  arrange(cluster)

write.csv(cluster_genes_mhal, file="cluster_genes_mhal.csv")




hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows = 4,
                      annotation_row =data.frame(
                        row.names = cluster_gene$ID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))

my_graph

# Save 
ggsave(filename="out/RNAseq_mhal/cluster_4_heatmap.pdf", plot=hm, width =4 , height = 3) 




# optimal nb of cluster------

pkgs <- c("factoextra",  "NbClust")
#install.packages(pkgs)

library(factoextra)
library(NbClust)

df <- scale(diffcounts)


# Elbow method --> Give 4
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method --> Give 2
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic --> Give 1
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")


# facet line -----



cluster_gene_tidy <- cluster_gene %>%
  pivot_longer(!c(ID,cluster), names_to="type", values_to="count") %>%
  group_by(ID) %>%
  mutate(count=(count-min(count))/(max(count)-min(count)),cluster=as.character(cluster))


cluster_gene_tidy$type <- factor(cluster_gene_tidy$type, c("col0","RNAi0","col4","RNAi4"))

cluster_gene_tidy_stat <- cluster_gene_tidy %>%
  group_by(cluster,type)%>%
  summarise(mean=mean(count),n=n()) %>%
  select(cluster,n) %>% unique

title_cluster <- data.frame (cluster  = c("1","2","3","4"),
                  cluster_clear = c("cluster 1 (2,446 genes)", "cluster 2 (2,867 genes)", "cluster 3 (27 genes)", "cluster 4 (3 genes)")
)

my_graph <- 
cluster_gene_tidy %>% 
  left_join(title_cluster) %>%
  ggplot(data=., mapping=aes(x=type, y=count, group=ID)) +
  geom_line(size=0.75, colour="grey") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "normcount")+
  facet_wrap(~cluster_clear,scales = "free",nrow=1)+
  stat_summary(aes(y = count,group=1), fun.y=mean, colour="red", geom="line",group=1)
my_graph



# Save 
ggsave(filename="out/RNAseq_mhal/cluster_4_profile.pdf", plot=my_graph, width = 7, height = 2) 

# corected facet_line -----------

scaled_cluster_gene



scaled_cluster_gene_tidy <- scaled_cluster_gene %>%
  pivot_longer(!c(ID,cluster), names_to="type", values_to="count") %>%
  mutate(cluster=as.character(cluster))




scaled_cluster_gene_tidy$type <- factor(scaled_cluster_gene_tidy$type, c("col0","RNAi0","col4","RNAi4"))

scaled_cluster_gene_stat <- scaled_cluster_gene_tidy %>%
  group_by(cluster,type)%>%
  summarise(mean=mean(count),n=n()) %>%
  select(cluster,n) %>% unique

title_cluster <- data.frame (cluster  = c("1","2","3","4"),
                             cluster_clear = c("cluster 1 (2,446 genes)", "cluster 2 (2,867 genes)", "cluster 3 (27 genes)", "cluster 4 (3 genes)")
)

my_graph <- 
  scaled_cluster_gene_tidy %>% 
  left_join(title_cluster) %>%
  ggplot(data=., mapping=aes(x=type, y=count, group=ID)) +
  geom_line(size=0.75, colour="grey") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "normcount")+
  facet_wrap(~cluster_clear,scales = "free",nrow=1)+
  stat_summary(aes(y = count,group=1), fun.y=mean, colour="red", geom="line",group=1)
print(my_graph)



# Save 
ggsave(filename="out/RNAseq_mhal/cluster_4_profile.pdf", plot=my_graph, width = 7, height = 2) 


#number of DE genes BH0.01 original ------------

Col_t0h_vs_Col_t4h_diff_expressed <- read_delim("data/RNAseq_mhal/BH0.01/out/Col_t0h_vs_Col_t4h_diff_expressed.tsv", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed <- read_delim("data/RNAseq_mhal/BH0.01/out/RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed.tsv", 
                                                          "\t", escape_double = FALSE, trim_ws = TRUE)


# individualize the up and down in each genotypes
Col_up <- Col_t0h_vs_Col_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange>0) %>%
  add_column(genotype ="Col",
             DE="up") %>%
  unique()
Col_down <- Col_t0h_vs_Col_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange<0) %>%
  add_column(genotype ="Col",
             DE="down")%>%
  unique()

RNAi_up <- RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange>0) %>%
  add_column(genotype ="RNAi",
             DE="up")%>%
  unique()
RNAi_down <- RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange<0) %>%
  add_column(genotype ="RNAi",
             DE="down")%>%
  unique()

all_de_RNAi_col <- Col_up %>%
  bind_rows(Col_down) %>%
  bind_rows(RNAi_up) %>%
  bind_rows(RNAi_down) %>%
  group_by(genotype,DE) %>%
  summarise(n=n())


write.csv(all_de_RNAi_col, file="all_de_RNAi_col.csv")

my_graph <- 
ggplot(all_de_RNAi_col, aes(genotype,n, label=n))+
  geom_col(aes(fill=DE),position=position_dodge(0.9))+
  xlab("")+
  ylab("count")
my_graph

# Save 
ggsave(filename="out/RNAseq_mhal/numberDEgenes.pdf", plot=my_graph, width = 3, height = 2.5) 

my_graph <- 
ggplot(all_de_RNAi_col, aes(DE,n, label=n))+
  geom_col(aes(fill=genotype),position=position_dodge(0.9))+
  xlab("")+
  ylab("count")

my_graph
# Save 
ggsave(filename="out/RNAseq_mhal/numberDEgenesV2.pdf", plot=my_graph, width = 3, height = 2.5) 






#number of DE genes FDR0.01 new ------------

Col_t0h_vs_Col_t4h_diff_expressed <- read_delim("data/RNAseq_mhal/FDR0.01_MARS corrected/Col_t0h_vs_Col_t4h_stat_expressed.tsv", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed <- read_delim("data/RNAseq_mhal/FDR0.01_MARS corrected/RNAiMHAL_t0h_vs_RNAiMHAL_t4h_stat_expressed.tsv", 
                                                          "\t", escape_double = FALSE, trim_ws = TRUE)


# individualize the up and down in each genotypes
Col_up <- Col_t0h_vs_Col_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange>0, padj<0.01) %>%
  add_column(genotype ="Col",
             DE="up") %>%
  unique()
Col_down <- Col_t0h_vs_Col_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange<0, padj<0.01) %>%
  add_column(genotype ="Col",
             DE="down")%>%
  unique()

RNAi_up <- RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange>0, padj<0.01) %>%
  add_column(genotype ="RNAi",
             DE="up")%>%
  unique()
RNAi_down <- RNAiMHAL_t0h_vs_RNAiMHAL_t4h_diff_expressed %>%
  dplyr::select(gene_id, log2FoldChange,padj) %>%
  filter(log2FoldChange<0, padj<0.01) %>%
  add_column(genotype ="RNAi",
             DE="down")%>%
  unique()

all_de_RNAi_col <- Col_up %>%
  bind_rows(Col_down) %>%
  bind_rows(RNAi_up) %>%
  bind_rows(RNAi_down) %>%
  group_by(genotype,DE) %>%
  summarise(n=n())


write.csv(all_de_RNAi_col, file="all_de_RNAi_col.csv")

my_graph <- 
  ggplot(all_de_RNAi_col, aes(genotype,n, label=n))+
  geom_col(aes(fill=DE),position=position_dodge(0.9))+
  xlab("")+
  ylab("count")
my_graph

# Save 
ggsave(filename="out/RNAseq_mhal/numberDEgenes.pdf", plot=my_graph, width = 3, height = 2.5) 

my_graph <- 
  ggplot(all_de_RNAi_col, aes(DE,n, label=n))+
  geom_col(aes(fill=genotype),position=position_dodge(0.9))+
  xlab("")+
  ylab("count")

my_graph
# Save 
ggsave(filename="out/RNAseq_mhal/numberDEgenesV2_FDR0.01MARS corr.pdf", plot=my_graph, width = 3, height = 2.5) 




# graph per gene-------




normcounts_tidy <- normcounts %>% 
  gather("Col_t0h_1":"RNAiMHAL_t4h_3", key = "condition", value = "counts", na.rm = TRUE) %>%
  separate(condition, c("genotype", "time", "replicate"), sep="_") %>%
  dplyr::select(ID, genotype, replicate, counts, time) %>%
  filter(replicate %in% c(3,2))


stat_normcounts_tidy <- normcounts_tidy %>%
  group_by(ID, genotype, time) %>%
  summarise(mean_counts=mean(counts), 
            median_counts=median(counts),
            ecart_type=sd(counts), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart


# List gene  -------------

gene_list <- c("AT5G42600", "AT5G42580", "AT5G42590", "AT5G00580")
gene_list <- c("AT5G42560", "AT5G42570", "AT5G42620", "AT5G42610")
gene_list <- c("AT3G30720","AT5G65080","AT1G70580","AT1G79600",
               "AT2G15620","AT2G32290","AT2G42600","AT3G21670",
               "AT3G61220","AT4G34480","AT4G34550","AT5G15500")

gene_list <- c("AT5G66400", "AT5G52300","AT5G42600", "AT5G42580", "AT5G42590", "AT5G00580")


   


gene_list <- c("AT3G30720", "AT1G70580", "AT1G79600","AT2G15620","AT2G42600")
gene_list <- c("AT1G67105", "AT3G13277", "AT4G08035","AT5G18255","AT1G46554","AT5G03615","AT1G05213","AT4G13495")



  filter(ID %in% gene_list)
normcounts_current <- normcounts_tidy %>%


normcounts_current %>%
  filter(time=="t4h")%>%
  ggbarplot(., x = "genotype", y = "counts",add = "mean_se", fill="ID") +
  stat_compare_means(method="t.test", 
                     comparisons = list(c("Col","RNAiMHAL") ), 
                     label = "p.label", bracket.size= 0.5)+
  facet_wrap(~ID,nrow=1,scale="free")



stat_normcounts_current <-
  stat_normcounts_tidy %>%
  filter(ID %in% gene_list) # Stat que sur les genes en ?tudes



ggplot(data=stat_normcounts_current, aes(x=time, color=genotype)) +
  geom_jitter(data=normcounts_current, aes(y=counts, shape=replicate)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free") +
  theme_classic()



my_graph <- 
ggplot(data=stat_normcounts_current, aes(x=time, fill=genotype)) +
  geom_col(aes(y=mean_counts), position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=0.25, position =position_dodge(0.9)) +
  facet_wrap(~ID, scale="free",nrow=1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.3))

my_graph
# Save 
ggsave(filename="out/RNAseq_mhal/norm_count_neigb.pdf", plot=my_graph, width = 8, height = 2.5) 

ggsave(filename="out/RNAseq_mhal/norm_count_target2.pdf", plot=my_graph, width = 7, height = 6) 
ggsave(filename="out/RNAseq_mhal/norm_count_marneral cluster.pdf", plot=my_graph, width = 7, height = 6) 
ggsave(filename="out/RNAseq_mhal/norm_count_marneral cluster ABA_corr.pdf", plot=my_graph, width = 8, height = 2.25) 

ggsave(filename="out/RNAseq_mhal/norm_count_marneral cluster_corr.pdf", plot=my_graph, width = 7, height = 6) 



# Matt Prior cluster --------
#####SP1------------
SP_counts <- read_delim("D:/Fulbright_UCR/data/Matt/SP_counts.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ID=X1)


diffcounts <- SP_counts %>%
  drop_na() %>%
  mutate(Total = select(., NRTLEAFTRA1:`35SLEAFTOT2`) %>% rowSums(na.rm = TRUE)) %>%
  filter(Total > 10) %>%
  unique() %>%
  select(-Total) %>%
  select("ID","35SLEAFTOT1", "35SLEAFTOT2","35SLEAFTRA1","gSWT13LEAFTOT70L912","gSWT13LEAFTOT100L912","gSWT13LEAFTOT3L12","gSWT13LEAFTRA1",
"NRTLEAFTOT1L74","NRTLEAFTOT2L74","NRTLEAFTRA1","GSD1LEAFTOT1L31","GSD1LEAFTOT2L31","GSD1LEAFTRA2L31","GSD1LEAFTRA1L31",
"NRTSCEPTOT1L2374","gSWT11LEAFRECL16","gSWT11SCEPTOT1","gSWT13STEMTOT1","35SSTEMPTOT1")



diffcounts <- as.data.frame(diffcounts) 

row.names(diffcounts)=diffcounts$ID
diffcounts$ID=NULL


library(pheatmap)

heatmap=pheatmap(diffcounts,scale = "row",cluster_cols = F )


saveRDS(heatmap, "heatmap.rds")

cluster_gene=cbind(diffcounts,
                   cluster = cutree(heatmap$tree_row,
                                    k = 8))%>%add_rownames(var="ID")%>%arrange(desc(cluster))


hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows =8,
                      annotation_row =data.frame(
                        row.names = cluster_gene$ID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))


write.csv(cluster_gene, file="cluster_gene.csv")

# Save 
ggsave(filename="out/ChIRP/1_8cl.pdf", plot=hm, width = 5, height = 4) 


# optimal nb of cluster------

pkgs <- c("factoextra",  "NbClust")
#install.packages(pkgs)

library(factoextra)
library(NbClust)

df <- scale(diffcounts)


# Elbow method --> Give ??
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method --> Give 2
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic --> Give 1
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")




#SP2-------------


SP_counts <- read_delim("D:/Fulbright_UCR/data/Matt/SP_counts2.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ID=X1)


diffcounts <- SP_counts %>%
  drop_na() %>%
  mutate(Total = select(., `35SLEAFTOT1`:GSD1LEAFTRA2L31) %>% rowSums(na.rm = TRUE)) %>%
  filter(Total > 10) %>%
  unique() %>%
  select(-Total) %>%
  select("ID","35SLEAFTOT1","35SLEAFTRA1","NRTLEAFTOT1L74","NRTLEAFTRA1","gSWT13LEAFTOT70L912","gSWT13LEAFTRA1","GSD1LEAFTOT1L31",
         "GSD1LEAFTRA1L31","GSD1LEAFTRA2L31")



diffcounts <- as.data.frame(diffcounts) 

row.names(diffcounts)=diffcounts$ID
diffcounts$ID=NULL


library(pheatmap)

heatmap=pheatmap(diffcounts,scale = "row",cluster_cols = F )


saveRDS(heatmap, "heatmap.rds")

cluster_gene=cbind(diffcounts,
                   cluster = cutree(heatmap$tree_row,
                                    k = 8))%>%add_rownames(var="ID")%>%arrange(desc(cluster))


hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows = 8,
                      annotation_row =data.frame(
                        row.names = cluster_gene$ID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))



# Save 
ggsave(filename="out/ChIRP/2_8cl_without2222.pdf", plot=hm, width = 5, height = 4) 




# ARES ---------

# AFFIMETRIX ------

# LogFC -----
ARES_DEG_filteredonlyAGIinAffimetrixATH1_nodup_merged_NAA_dataset <- 
  read_csv("C:/Users/roule/Box/Ongoing projects/ARES/R clustering/in/ARES_DEG_filteredonlyAGIinAffimetrixATH1_nodup_merged_NAA_dataset.csv") %>%
  select(geneID,FoldChange_Col0_2h,FoldChange_arf19_2h,FoldChange_arf7_2h, FoldChange_arf19_arf7_2h)




diffcounts <- as.data.frame(ARES_DEG_filteredonlyAGIinAffimetrixATH1_nodup_merged_NAA_dataset) 

row.names(diffcounts)=diffcounts$geneID
diffcounts$geneID=NULL


library(pheatmap)

heatmap=pheatmap(diffcounts,scale = "row",cluster_cols = F )




saveRDS(heatmap, "heatmap.rds")

cluster_gene=cbind(diffcounts,
                   cluster = cutree(heatmap$tree_row,
                                    k = 9))%>%add_rownames(var="geneID")%>%arrange(desc(cluster))


#SCALE
scaled_cluster_gene=cbind(diffcounts %>% as.matrix() %>% t %>% scale() %>% t,
                          cluster = cutree(heatmap$tree_row, k = 9))%>%
  as.data.frame %>% add_rownames(var="geneID")%>%arrange(desc(cluster))



# collect genes from the cluster
cluster_genes_scale <- cluster_gene

write.csv(cluster_genes_scale, file="cluster_genes_scale_4.csv")

##

#Replace scale by none or row for unscale or scale:
hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "none",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows = 9,
                      annotation_row =data.frame(
                        row.names = cluster_gene$geneID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))

my_graph

# Save 
ggsave(filename="out/RNAseq_mhal/ARES_cluster_4_heatmap_unscale.pdf", plot=hm, width =4 , height = 10) 




# optimal nb of cluster------

pkgs <- c("factoextra",  "NbClust")
#install.packages(pkgs)

library(factoextra)
library(NbClust)

df <- scale(diffcounts)


# Elbow method --> Give 6
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method --> Give 9
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic --> Give 9
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")











# corected facet_line -----------

scaled_cluster_gene



scaled_cluster_gene_tidy <- scaled_cluster_gene %>%
  pivot_longer(!c(ID,cluster), names_to="type", values_to="count") %>%
  mutate(cluster=as.character(cluster))




scaled_cluster_gene_tidy$type <- factor(scaled_cluster_gene_tidy$type, c("col0","RNAi0","col4","RNAi4"))

scaled_cluster_gene_stat <- scaled_cluster_gene_tidy %>%
  group_by(cluster,type)%>%
  summarise(mean=mean(count),n=n()) %>%
  select(cluster,n) %>% unique

title_cluster <- data.frame (cluster  = c("1","2","3","4"),
                             cluster_clear = c("cluster 1 (2,446 genes)", "cluster 2 (2,867 genes)", "cluster 3 (27 genes)", "cluster 4 (3 genes)")
)

my_graph <- 
  scaled_cluster_gene_tidy %>% 
  left_join(title_cluster) %>%
  ggplot(data=., mapping=aes(x=type, y=count, group=ID)) +
  geom_line(size=0.75, colour="grey") +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "normcount")+
  facet_wrap(~cluster_clear,scales = "free",nrow=1)+
  stat_summary(aes(y = count,group=1), fun.y=mean, colour="red", geom="line",group=1)
print(my_graph)



# Save 
ggsave(filename="out/RNAseq_mhal/cluster_4_profile.pdf", plot=my_graph, width = 7, height = 2) 







