normcounts <- read_delim("heatmap data/normcounts.csv", 
                         ";", escape_double = FALSE, trim_ws = TRUE)  


#clustering ------

clustering_DE <- read_excel("heatmap data/DE genes.xlsx")


diffcounts <- normcounts %>% inner_join(clustering_DE)
diffcounts <- diffcounts %>% 
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



hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows = 4,
                      annotation_row =data.frame(
                        row.names = cluster_gene$ID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))



## optimal number of cluster ----


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



