# Out directory if it does not exist
dir.create("out/Pearson/", showWarnings = FALSE, recursive = TRUE)


library("RColorBrewer")
library(tidyverse)
library(corrplot)
library(readr)
library(readxl)
library(ggpubr)
theme_set(theme_bw())

normcounts <- read.table("./raw/in/normcounts.txt") %>%
  rownames_to_column("ID")
GO_1 <- read_tsv("./raw/in/GO root.tsv")
GO_3 %>% unique()
GO_2 <- read_tsv("./raw/in/GO meristem.tsv")
GO_3 <- read_tsv("./raw/in/GO cell division.tsv")
all_info_V201810415 <- read_excel("./raw/in/all_info_V201810415.xlsx")

GO_all <- GO_1 %>%
  bind_rows(GO_2) %>%
  bind_rows(GO_3) %>%
  dplyr::select(-"Gene.primaryIdentifier") %>%
  unique()

write.csv(GO_all, file="./GO_all.csv")

#Identification des lincRNA -----



lincRNA <- all_info_V201810415 %>%
  filter(NAT == "FALSE", 
         width > 200, 
         type == "noncoding") %>%
  select(ID, Araport11_ID)
lincRNA %>% unique()


#Nombre de lncRNA avec voisins GO HORS SUJET ------


lincRNA_num <- all_info_V201810415 %>%
  filter(NAT == "FALSE", 
         width > 200, 
         type == "noncoding") %>%
  select(ID, Araport11_ID, upstream_ID, downstream_ID)

GO_root_genes <- GO_1 %>% 
  select(Araport11_ID = Gene.primaryIdentifier) %>%
  mutate(GeneOntology = "root")
GO_meristem_genes <- GO_2 %>% 
  select(Araport11_ID = Gene.primaryIdentifier) %>%
  mutate(GeneOntology = "meristem")
GO_cell_genes <- GO_3 %>% 
  select(Araport11_ID = Gene.primaryIdentifier) %>%
  mutate(GeneOntology = "cell")
#write.csv(GO_cell_genes, "./VennDiagramGO/GO_cell_genes.csv")

GO_all_genes <-
  bind_rows(GO_root_genes, GO_meristem_genes, GO_cell_genes) %>%
  unique %>%
  rename(coding = Araport11_ID)



lincRNA_number_up <- lincRNA_num %>% 
  select(upstream_ID, Araport11_ID, ID) %>%
  rename(coding = upstream_ID) %>%
  inner_join(GO_all_genes)
lincRNA_number_down_up <- lincRNA_num %>% 
  select(downstream_ID, Araport11_ID, ID) %>%
  rename(coding = downstream_ID) %>%
  inner_join(GO_all_genes) %>% 
  bind_rows(lincRNA_number_up)


###################









# Identification des genes avec GO ------


GO_root_genes <- GO_1 %>% 
  select(Araport11_ID = Gene.primaryIdentifier)
GO_meristem_genes <- GO_2 %>% 
  select(Araport11_ID = Gene.primaryIdentifier)
GO_cell_genes <- GO_3 %>% 
  select(Araport11_ID = Gene.primaryIdentifier)
GO_all_genes <-
  bind_rows(GO_root_genes, GO_meristem_genes, GO_cell_genes) %>%
  unique 

coding_GO <-
  all_info_V201810415 %>%
  filter(!is.na(Araport11_ID), # Non s?lection des Araport11_ID vide
         Araport11_ID %in% GO_all_genes$Araport11_ID) %>% # dans colonne Araport_11 ID, selection des ID avec GO
  select(ID, Araport11_ID) %>%
  unique 



#S?lection des g?nes diff?rentiellement exprim? ------


gene_diff <- all_info_V201810415 %>%
  filter(t1t2 != "NA" | t0t1 != "NA" | t0t2 != "NA" | Col_Ler != "NA") %>%
  select(ID, Araport11_ID)


gene_diff_Pi <- all_info_V201810415 %>%
  filter(t1t2 != "NA" | t0t1 != "NA" | t0t2 != "NA") %>%
  select(ID, Araport11_ID)

gene_diff_ecotype <- all_info_V201810415 %>%
  filter(Col_Ler != "NA") %>%
  select(ID, Araport11_ID)


#lincRNA
nc_gene_count <-
  normcounts %>%
  filter(ID %in% lincRNA$ID,
         ID %in% gene_diff$ID) %>% 
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") 


#coding gene avec GO

coding_gene_count <-
  normcounts %>%
  filter(ID %in% coding_GO$ID,
         ID %in% gene_diff$ID) %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID")


# Correlation ------


pearson_correlation_df <- 
  cor(t(nc_gene_count), t(coding_gene_count), method="pearson") %>% 
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>% 
  filter(ncRNA %in% row.names(nc_gene_count), 
         gene %in% row.names(coding_gene_count))  




my_graph <- 
ggplot(pearson_correlation_df, aes(x=NA, y=pearson_corr)) + geom_violin()

my_graph

# Save 
ggsave(filename="out/Pearson/col and ler pearson.pdf", plot=my_graph, width = 7, height = 5)


#les 100 1er couples positifs/négatifs sont-ils Pi responsive? - HORHS SIJET -----------
#codant
pearson_correlation_df %>% 
  arrange(desc(pearson_corr)) %>% 
  head(100) %>% 
  select(gene) %>% 
  rename(Araport11_ID=gene) %>%
  inner_join(gene_diff_Pi) %>% 
  unique() %>% 
  nrow()

pearson_correlation_df %>% 
  arrange(desc(pearson_corr)) %>% 
  head(100) %>% 
  select(gene) %>% 
  rename(Araport11_ID=gene) %>%
  inner_join(gene_diff_ecotype) %>% 
  unique()%>% 
  nrow()
#ncRNA
pearson_correlation_df %>% 
  arrange(desc(pearson_corr)) %>% 
  head(100) %>% 
  select(ncRNA) %>% 
  rename(Araport11_ID=ncRNA) %>%
  inner_join(gene_diff_Pi) %>% 
  unique()%>% 
  nrow()

pearson_correlation_df %>% 
  arrange(desc(pearson_corr)) %>% 
  head(100) %>% 
  select(ncRNA) %>% 
  rename(Araport11_ID=ncRNA) %>%
  inner_join(gene_diff_ecotype) %>% 
  unique()%>% 
  nrow()


#######



#Distribution des couples de corr?lations (ncRNA et coding gene with GO)**
#Trop de corr?lation car forte diff?rence entre comptage Col et comptage Ler --> cr?ation de faux positifs

nrow(pearson_correlation_df) # couples de corr?lations


## S?paration des comptages issus de Col et issus de Ler -------

nc_gene_count_Col <-
  normcounts %>%
  filter(ID %in% lincRNA$ID,
         ID %in% gene_diff$ID) %>% # Dans colonne ID, selection des ID des lncRNA diff exprim? (, = et)
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") %>%
  select(col_t0_A, col_t0_B,  col_t0_C, col_t1_A, col_t1_B, col_t1_C, col_t2_A, col_t2_B, col_t2_C)  # select que les count Col
nc_gene_count_Ler <-
  normcounts %>%
  filter(ID %in% lincRNA$ID,
         ID %in% gene_diff$ID) %>% # Dans colonne ID, selection des ID des lncRNA diff exprim? (, = et)
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") %>%
  select(ler_t0_A, ler_t0_B,  ler_t0_C, ler_t1_A, ler_t1_B, ler_t1_C, ler_t2_A, ler_t2_B, ler_t2_C) 

coding_gene_count_Col <-
  normcounts %>%
  filter(ID %in% coding_GO$ID,
         ID %in% gene_diff$ID) %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") %>%
  select(col_t0_A, col_t0_B,  col_t0_C, col_t1_A, col_t1_B, col_t1_C, col_t2_A, col_t2_B, col_t2_C)  # select que les count Col
coding_gene_count_Ler <-
  normcounts %>%
  filter(ID %in% coding_GO$ID,
         ID %in% gene_diff$ID) %>% # Dans colonne ID, selection des ID des coding_GO diff exprim? (, = et)
  unique() %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ID") %>%
  select(ler_t0_A, ler_t0_B,  ler_t0_C, ler_t1_A, ler_t1_B, ler_t1_C, ler_t2_A, ler_t2_B, ler_t2_C) 

correlation_Col <- 
  cor(t(nc_gene_count_Col), t(coding_gene_count_Col), method="pearson") %>% 
  as.data.frame  %>% 
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>% 
  filter(ncRNA %in% row.names(nc_gene_count_Col), 
         gene %in% row.names(coding_gene_count_Col)) %>%
  rename(pearson_corr_col = pearson_corr)
correlation_Ler <- # correlation des count ncRNA et count coding gene
  cor(t(nc_gene_count_Ler), t(coding_gene_count_Ler), method="pearson") %>% 
  as.data.frame  %>% 
  rownames_to_column("ncRNA") %>% 
  gather(gene, pearson_corr, -ncRNA) %>% 
  filter(ncRNA %in% row.names(nc_gene_count_Ler), 
         gene %in% row.names(coding_gene_count_Ler)) %>%
  rename(pearson_corr_ler = pearson_corr)

#
#les 100 1er couples positifs/négatifs sont-ils Pi responsive? - HORS SUJET -------

correlation_Col %>% 
  arrange(desc(pearson_corr_col)) %>% 
  head(100) %>% 
  select(gene) %>% 
  rename(Araport11_ID=gene) %>%
  inner_join(gene_diff_Pi) %>% 
  unique()%>% nrow()

correlation_Col %>% 
  arrange(desc(pearson_corr_col)) %>% 
  head(100) %>% 
  select(gene) %>% 
  rename(Araport11_ID=gene) %>%
  inner_join(gene_diff_ecotype) %>% 
  unique()%>% nrow()
#ncRNA
correlation_Col %>% 
  arrange(desc(pearson_corr_col)) %>% 
  head(100) %>% 
  select(ncRNA) %>% 
  rename(Araport11_ID=ncRNA) %>%
  inner_join(gene_diff_Pi) %>% 
  unique()%>% nrow()

correlation_Col %>% 
  arrange(desc(pearson_corr_col)) %>% 
  head(100) %>% 
  select(ncRNA) %>% 
  rename(Araport11_ID=ncRNA) %>%
  inner_join(gene_diff_ecotype) %>% 
  unique()%>% nrow()




# corelation comptage séparé ------



my_graph <- 
ggplot(correlation_Col, aes(x=NA, y=pearson_corr_col)) + geom_violin()
my_graph
# Save 
ggsave(filename="out/Pearson/Pearson_col.pdf", plot=my_graph, width = 7, height = 5)


my_graph <- 
ggplot(correlation_Ler, aes(x=NA, y=pearson_corr_ler)) + geom_violin()
my_graph 
# Save 
ggsave(filename="out/Pearson/Pearson_ler.pdf", plot=my_graph, width = 7, height = 5)


#Bonne distribution des corr?lations


# Selection de candidats --------


#*Recherche de cluster parmis les 100 1eres correlations les plus forte :
#Sur les comptage de Col_Positives -----
  


correlation_col_100 <- correlation_Col %>% arrange(desc(pearson_corr_col)) %>% filter(row_number() <= 100) %>% as.data.frame()

#correlation_ler_100 <- correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 100) %>% as.data.frame() 
#correlation_all <- correlation_Col %>% full_join(correlation_Ler)
#correlation_col_100 %>% group_by(ncRNA) %>% summarise(mean=mean(pearson_corr_col),n=n()) %>% as.data.frame()



correlation_col_100 %>% 
  group_by(ncRNA) %>% 
  summarise(n=n(),mean=mean(pearson_corr_col)) %>% 
  as.data.frame()%>% 
  arrange(desc(n)) %>% 
  filter(n>=5) %>%
  add_column(Corr="Pi")



correlation_col_100 %>% filter(ncRNA == "AT5G09710") %>% select(gene)
correlation_col_100 %>% filter(ncRNA == "AT1G55525") %>% select(gene)
correlation_col_100 %>% filter(ncRNA == "XLOC_002093") %>% select(gene)
correlation_col_100 %>% filter(ncRNA == "XLOC_002755") %>% select(gene)
correlation_col_100 %>% filter(ncRNA == "XLOC_002900") %>% select(gene)
correlation_col_100 %>% filter(ncRNA == "XLOC_005603") %>% select(gene)




















































#AT1G55525 lncRNA


correlation_col_100 %>% filter(ncRNA == "AT1G55525") %>% select(gene) %>% as.list()
gene_list <- c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")


C1 <- correlation_all %>% 
  filter(ncRNA == "AT1G55525", gene %in% c("AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220")) %>% 
  add_column(type = "Col_count_Pos_100")
C1
``` 
*Araport
Unknown gene
GO AIA, root

AT1G19220	ARF19	auxin response factor 19
AT1G50460	HKL1	hexokinase-like 1
AT1G55580	LAS	GRAS family transcription factor
AT1G58340	ZF14	MATE efflux family protein
AT1G78240	TSD2	S-adenosyl-L-methionine-dependent methyltransferases superfamily protein
AT2G34650	PID	Protein kinase superfamily protein
AT3G03660	WOX11	WUSCHEL related homeobox 11
AT3G14370	WAG2	Protein kinase superfamily protein
AT3G60630	HAM2	GRAS family transcription factor
AT4G16780	HB-2	homeobox protein 2
AT4G17460	HAT1	Homeobox-leucine zipper protein 4 (HB-4) / HD-ZIP protein
AT4G34410	RRTF1	redox responsive transcription factor 1
AT5G52050	MATE efflux family protein

*GV
L?g?re corr?lation, surtout avec WOX11


*clustering GLM
``` {r, echo=FALSE}
gene_list <- c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220") %>% as.data.frame()
``` 

1  AT1G55525
2  AT5G52050 9
3  AT4G34410 9
4  AT2G34650 6
5  AT1G58340 6
6  AT1G55580
7  AT3G03660 9
8  AT1G78240 2
9  AT3G60630 2
10 AT4G17460 6
11 AT4G16780 6
12 AT3G14370 6
13 AT1G50460 2
14 AT1G19220 2


*Gem2Net
Pas le lncNRA



**Oui**
  
  ``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220")) %>% arrange(desc(type))
``` 
R?pondent ? cinetique Pi
Bcp de GO root dev


*AT3G52670 lncRNA*
  ```{r, include=FALSE}
correlation_col_100 %>% filter(ncRNA == "AT3G52670") %>% select(gene) %>% as.list()
gene_list <- c("AT3G52670", "AT5G43080", "AT4G33495") # On vire At5G43080 car correlation moche
gene_list <- c("AT3G52670", "AT4G33495")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT3G52670", gene == "AT4G33495")
``` 
*Araport
Pseudogene
FBD%2C F-box%2C Skp2-like and Leucine Rich Repeat domains containing protein  


Risque que ?a soit une prot?ine --> OUT
**Non**
  
#AT3G61198 lncRNA ------


correlation_col_100 %>% filter(ncRNA == "AT3G61198") %>% select(gene) %>% as.list()
gene_list <- c("AT3G61198", "AT4G08920", "AT5G03150", "AT1G72150", "AT1G22530")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")



#Belle corr?lation mais d?j? s?lectionner pour son action en CIS sur sa cible (BAP1). Niveau d'expression tr?s fort des TRANS candidats en plus...

#AT4G13495 lncRNA


correlation_col_100 %>% filter(ncRNA == "AT4G13495") %>% select(gene) %>% as.list()
gene_list <- c("AT4G13495", "AT1G28560")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")






correlation_all %>% filter(ncRNA == "AT4G13495", gene == "AT1G28560")

Positivement corr?l? chez Col
Negativement corr?l? chez Ler




**Non**


#AT5G09710 lncRNA

  
correlation_col_100 %>% filter(ncRNA == "AT5G09710") %>% select(gene) %>% as.list()
gene_list <- c("AT5G09710", "AT3G60630", "AT3G07390", "AT1G70560", "AT1G78240", "AT4G34410", "AT5G52050", "AT3G03660", "AT1G17110", "AT3G14370", "AT1G19220", "AT1G58340", "AT2G36400", "AT5G10720", "AT4G17460", "AT2G34650", "AT5G59030", "AT1G70940", "AT1G55580", "AT4G16780", "AT1G50460", "AT4G09510")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  


correlation_all %>% filter(ncRNA == "AT5G09710", gene %in% c("AT3G60630", "AT3G07390", "AT1G70560", "AT1G78240", "AT4G34410", "AT5G52050", "AT3G03660", "AT1G17110", "AT3G14370", "AT1G19220", "AT1G58340", "AT2G36400", "AT5G10720", "AT4G17460", "AT2G34650", "AT5G59030", "AT1G70940", "AT1G55580", "AT4G16780", "AT1G50460", "AT4G09510"))


Belle corr?lation. Similaire entre Col et Ler

*Araport
Pseudogene
Magnesium transporter CorA-like family protein; BEST Arabidopsis thaliana protein match is: magnesium transporter 7 (TAIR:AT5G09690.2)

Probablement une prot?ine ou pseudogene : trop risqu? --> OUT

**Non**





#XLOC_002093/AT1G08925 lncRNA*

correlation_col_100 %>% filter(ncRNA == "XLOC_002093") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")



C2 <- correlation_all %>% filter(ncRNA == "XLOC_002093", gene %in% c("AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")) %>% add_column(type = "Col_count_Pos_100")
C2
``` 
Belle corr?lation. Similaire entre Col et Ler

*Araport
novel transcribed region

AT1G17110	UBP15	ubiquitin-specific protease 15
AT1G19220	ARF19	auxin response factor 19
AT1G55580	LAS	GRAS family transcription factor
AT1G58340	ZF14	MATE efflux family protein
AT1G70560	TAA1	tryptophan aminotransferase of Arabidopsis 1
AT1G78240	TSD2	S-adenosyl-L-methionine-dependent methyltransferases superfamily protein
AT3G03660	WOX11	WUSCHEL related homeobox 11
AT3G07390	AIR12	auxin-induced in root cultures-like protein
AT3G60630	HAM2	GRAS family transcription factor
AT4G14550	IAA14	indole-3-acetic acid inducible 14
AT4G34410	RRTF1/ERF	redox responsive transcription factor 1
AT5G52050	AT5G52050	MATE efflux family protein
AT5G59030	COPT1	copper transporter 1

GO root


*GV
AT4G34410
AT5G52050

bien cor?l? (cluster 9)


*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

1  XLOC_002093
2    AT1G78240  2
3    AT3G07390  6
4    AT3G03660  9
5    AT1G70560  6
6    AT4G34410  9
7    AT1G58340  6
8    AT3G60630  2
9    AT1G17110  10
10   AT1G19220  2
11   AT5G52050  9
12   AT5G59030  6
13   AT1G55580
14   AT4G14550  2

*Gem2Net
Pas le lncNRA

**Oui**, Limiter la corr?lation aux 2 cluster 9?


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")) %>% arrange(desc(type))
``` 
R?pondent ? cinetique Pi,
Bcp de GO root dev






#XLOC_002755/AT1G08173 lncRNA*
```{r, echo=FALSE}
correlation_col_100 %>% filter(ncRNA == "XLOC_002755") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002755", "AT2G23430", "AT1G72160", "AT5G58010", "AT1G67710", "AT3G04630", "AT1G13260", "AT2G34680", "AT1G72150", "AT3G22400", "AT4G33880", "AT4G34580", "AT2G31090")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  


C3 <- correlation_all %>% filter(ncRNA == "XLOC_002755", gene %in% c("AT2G23430", "AT1G72160", "AT5G58010", "AT1G67710", "AT3G04630", "AT1G13260", "AT2G34680", "AT1G72150", "AT3G22400", "AT4G33880", "AT4G34580", "AT2G31090")) %>% add_column(type = "Col_count_Pos_100")
C3


Belle corr?lation. Similaire entre Col et Ler

*Araport
Rien

AT1G08173	""	""
AT1G13260	RAV1	related to ABI3/VP1 1
AT1G67710	ARR11	response regulator 11
AT1G72150	PATL1	PATELLIN 1
AT1G72160	AT1G72160	Sec14p-like phosphatidylinositol transfer family protein
AT2G23430	ICK1	Cyclin-dependent kinase inhibitor family protein
AT2G31090	AT2G31090	transmembrane protein
AT2G34680	AIR9	Outer arm dynein light chain 1 protein
AT3G04630	WDL1	WVD2-like 1
AT3G22400	LOX5	PLAT/LH2 domain-containing lipoxygenase family protein
AT4G33880	RSL2	ROOT HAIR DEFECTIVE 6-LIKE 2
AT4G34580	COW1	Sec14p-like phosphatidylinositol transfer family protein
AT5G58010	LRL3	LJRHL1-like 3


GO root

*GV
pas le lncRNA
bof

*clustering GLM
```{r, include=FALSE}
gene_list <- c("AT1G08173", "AT2G23430", "AT1G72160", "AT5G58010", "AT1G67710", "AT3G04630", "AT1G13260", "AT2G34680", "AT1G72150", "AT3G22400", "AT4G33880", "AT4G34580", "AT2G31090") %>% as.data.frame()
``` 

1  AT1G08173  4
2  AT2G23430  4
3  AT1G72160  7
4  AT5G58010  4
5  AT1G67710  4
6  AT3G04630  4
7  AT1G13260  4
8  AT2G34680  7
9  AT1G72150  5
10 AT3G22400  7
11 AT4G33880  1
12 AT4G34580  4
13 AT2G31090  7

*Gem2Net
Pas le lncRNA

**Oui**, Limiter la corr?lation aux cluster 4?


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")) %>% arrange(desc(type))
``` 
Repondent a cinetique Pi
Bcp de GO root dev





*XLOC_002900/AT1G08687 lncRNA*
```{r, include = FALSE}
correlation_col_100 %>% filter(ncRNA == "XLOC_002900") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002900", "AT4G17500", "AT3G04630", "AT5G61600", "AT1G67710", "AT1G56010", "AT2G26670")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
```
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "XLOC_002900", gene %in% c("AT4G17500", "AT3G04630", "AT5G61600", "AT1G67710", "AT1G56010", "AT2G26670"))
``` 
Belle corr?lation. Similaire entre Col et Ler

*Araport
novel transcribed region

AT1G56010	NAC1	NAC domain containing protein 1
AT1G67710	ARR11	response regulator 11
AT2G26670	TED4	Plant heme oxygenase (decyclizing) family protein
AT3G04630	WDL1	WVD2-like 1
AT4G17500	ERF-1	ethylene responsive element binding factor 1
AT5G61600	ERF104	ethylene response factor 104


*GV
Pas le lncRNA

plutot cor?l?, surtout les 2 ERFs



*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002900", "AT4G17500", "AT3G04630", "AT5G61600", "AT1G67710", "AT1G56010", "AT2G26670") %>% as.data.frame()
``` 

1 XLOC_002900 8
2   AT4G17500 4
3   AT3G04630 4
4   AT5G61600 1
5   AT1G67710 4
6   AT1G56010 4
7   AT2G26670 4

*Gem2Net
Pas le lncRNA

**Non**, Corr?lation pas si belle




*XLOC_005603/AT3G03315 lncRNA*
```{r, echo=FALSE}
correlation_col_100 %>% filter(ncRNA == "XLOC_005603") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_005603", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
```
```{r, echo=FALSE}
C4 <- correlation_all %>% filter(ncRNA == "XLOC_005603", gene %in% c("XLOC_005603", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430")) %>% add_column(type = "Col_count_Pos_100")
C4
``` 
Belle corr?altion Col et Ler


*Araport
novel transcribed region

AT1G14740	TTA1	class I heat shock protein, putative (DUF1423)
AT1G16510	AT1G16510	SAUR-like auxin-responsive protein family
AT1G25220	ASB1	anthranilate synthase beta subunit 1
AT2G42430	LBD16	lateral organ boundaries-domain 16
AT3G03315	""	novel transcribed region; detected in root, root apical meristem, dark-grown seedling, light-grown seedling, aerial
AT5G12330	LRP1	Lateral root primordium (LRP) protein-like protein

GO root


*GV
Bien corr?l?


*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_005603", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430") %>% as.data.frame()
``` 

1 XLOC_005603
2   AT1G25220 2
3   AT1G14740 2
4   AT5G12330 6
5   AT1G16510 6
6   AT2G42430 6


*Gem2Net
Pas le lncRNA

**Oui**


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_005603", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430")) %>% arrange(desc(type))
``` 
Repondent ? cinetique Pi


*XLOC_007011/AT3G00400 lncRNA*
```{r, include=FALSE}
correlation_col_100 %>% filter(ncRNA == "XLOC_007011") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_007011", "AT5G39610")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
```
```{r}
correlation_all %>% filter(ncRNA == "XLOC_007011", gene == "AT5G39610")
``` 
Corr?lation col et ler, profil diff?rent


*Araport
novel transcribed region

Coding = Encodes a NAC-domain transcription factor. Positively regulates aging-induced cell death and senescence in leaves. This gene is upregulated in response to salt stress in wildtype as well as NTHK1 transgenic lines although in the latter case the induction was drastically reduced. It was also upregulated by ABA, ACC and NAA treatment, although in the latter two cases, the induction occurred relatively late when compared with NaCl or ABA treatments. Note: this protein (AtNAC6) on occasion has also been referred to as AtNAC2, not to be confused with the AtNAC2 found at locus AT3G15510

*GV


*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

XLOC_007011 8
AT5G39610   9

*Gem2Net
not lncRNA

**Non**, coding pas sexy





*XLOC_007825/AT4G04595 lncRNA*
```{r, echo=FALSE}
correlation_col_100 %>% filter(ncRNA == "XLOC_007825") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_007825", "AT1G05630")
normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
```
```{r, echo=FALSE}
C5 <- correlation_all %>% filter(ncRNA == "XLOC_007825", gene == "AT1G05630") %>% add_column(type = "Col_count_Pos_100")
C5
``` 
Correlation positive chez Col
Negative chez Ler


*Araport
novel transcribed region

The WD40 repeat region of 5PTase13 interacts with the sucrose nonfermenting-1-related kinase, and loss of function in 5PTase13 leads to nutrient level-dependent reduction of root growth, along with abscisic acid (ABA) and sugar insensitivity. 
Arabidopsis 5PTase13 functions in root gravitropism by regulating vesicle trafficking and gravity-induced auxin redistribution.

*GV
Pas le lncRNA

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

AT1G05630 7

*Gem2Net
Pas le lncRNA

**Oui**, anticor?lation et coding sexy



``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_007825", "AT1G05630")) %>% arrange(desc(type))
``` 
Differentiel Col Ler
GO root dev




```{r}
correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 100) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
corr_count_col <- correlation_Col %>% arrange(desc(pearson_corr_col)) %>% filter(row_number() <= 100) %>% as.data.frame()
corr_count_ler <- correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 100) %>% as.data.frame()
corr_count_col_ler <- corr_count_col %>% inner_join(corr_count_ler) %>% select(ncRNA, gene) %>% group_by(ncRNA) %>% as.data.frame()
#Aucun couple retrouv? chez comptage Col et Ler parmis les 100 1eres corr?lations
```

### Negatives

```{r, echo=FALSE}
correlation_Col %>% arrange(pearson_corr_col) %>% filter(row_number() <= 100) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
```

```{r, include=FALSE}
correlation_col_100_neg <- correlation_Col %>% arrange(pearson_corr_col) %>% filter(row_number() <= 100) %>% as.data.frame()
correlation_ler_100_neg <- correlation_Ler %>% arrange(pearson_corr_ler) %>% filter(row_number() <= 100) %>% as.data.frame() 
```

*AT2G31585 lncRNA*
```{r, include=FALSE}
correlation_col_100_neg %>% filter(ncRNA == "AT2G31585") %>% select(gene) %>% as.list()
gene_list <- c("AT2G31585", "AT1G02330")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT1G02330"))
``` 
Anticorr?lation chez Col et Ler avec 2 profils diff?rents

*Araport
Unknown gene

*GV


*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 



*Gem2Net


**Non**, coding peu int?ressant


*AT4G34881 lncRNA*
```{r, echo=FALSE}
correlation_col_100_neg %>% filter(ncRNA == "AT4G34881") %>% select(gene) %>% as.list()
gene_list <- c("AT4G34881", "AT5G05560")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C6 <- correlation_all %>% filter(ncRNA == "AT4G34881", gene %in% c("AT5G05560")) %>% add_column(type = "Col_count_Neg_100")
C6
``` 
Anticorr?lation chez Col, non exprim? chez Ler --> target r?agit diff?rement



*Araport
transmembrane protein

Coding = Encodes a subunit of the Arabidopsis thaliana E3 ubiquitin ligase complex that plays a synergistic role with APC4 both in female gametogenesis and in embryogenesis.  --> GO AIA homeostasi

*GV
Non cor?l? (bien)

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 



*Gem2Net
Pas le lncRNA

**Oui**


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT4G34881", "AT5G05560")) %>% arrange(desc(type))
```
Col_Ler diff


*AT5G09710 lncRNA*
```{r, include=FALSE}
correlation_col_100_neg %>% filter(ncRNA == "AT5G09710") %>% select(gene) %>% as.list()
gene_list <- c("AT5G09710", "AT5G54670")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT5G09710", gene %in% c("AT5G54670"))
``` 
Anticorr?lation similaire chez Col et Ler


Probablement protein --> **Non**



*Col_NEW_RNA_F_30663/ lncRNA*
```{r, include = FALSE}
correlation_col_100_neg %>% filter(ncRNA == "Col_NEW_RNA_F_30663") %>% select(gene) %>% as.list()
gene_list <- c("Col_NEW_RNA_F_30663", "AT3G47690", "AT1G55540")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "Col_NEW_RNA_F_30663", gene %in% c("AT3G47690", "AT1G55540"))
``` 
Anticorr?lation chez Col. Non exprim? chez Ler --> Target r?agit diff?remment


*Araport
AT3G47690   encodes a homolog of animal microtubule-end-binding protein. There are two other members of this family. EB1 forms foci at regions where the minus ends of microtubules are gathered during mitosis and early cytokinesis
Data find that antibodies directed against EB1 proteins colocalize with microtubules in roots, and that mutants with reduced expression from EB1A genes have roots that deviate toward the left on vertical or inclined plates.

AT1G55540   LNO1 is required for seed viability in Arabidopsis.



*GV
les 2 coding non cor?l?

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

AT3G47690 10
AT1G55540 10


*Gem2Net
Pas le lncRNA, les coding pas retrouv? dans meme cluster

**Non**, anticor?lation pas si belle, car lncRNA expression peu variable en regardant l'?chelle




*Col_NEW_RNA_F_30663/AT1G09765 lncRNA*
  ```{r, include =FALSE}
correlation_col_100_neg %>% filter(ncRNA == "Col_NEW_RNA_F_7035") %>% select(gene) %>% as.list()
gene_list <- c("Col_NEW_RNA_F_7035", "AT5G05560")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "Col_NEW_RNA_F_7035", gene %in% c("AT5G05560"))
```
Anticorr?lation chez Col. Non exprim? chez Ler --> Target r?agit diff?remment


*Araport
novel transcribed region

AT5G05560 Encodes a subunit of the Arabidopsis thaliana E3 ubiquitin ligase complex that plays a synergistic role with APC4 both in female gametogenesis and in embryogenesis. 
GO AIA homesostasis


*GV
PAs le lncRNA

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

Col_NEW_RNA_F_7035  3
AT5G05560           10

Plutot diff?rent (bien)

*Gem2Net
PAs le lncRNA


**Non**, coding pas sexy





*XLOC_002093/AT2G00680 lncRNA*
  ```{r, include=FALSE}
correlation_col_100_neg %>% filter(ncRNA == "XLOC_002093") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002093", "AT5G54670", "AT4G33650", "AT1G30690", "AT2G25180")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "Col_NEW_RNA_R_10878", gene %in% c("AT5G43070"))
```
Anticorr?lation chez Col et Ler similaire

*Araport
novel transcribed region

AT1G30690	AT1G30690	Sec14p-like phosphatidylinositol transfer family protein
AT2G00680	""	novel transcribed region
AT2G25180	RR12	response regulator 12
AT4G33650	DRP3A	dynamin-related protein 3A
AT5G54670	ATK3	kinesin 3

AT2G25180, root CK, transduction


*GV
Pas le lncRNA
Plutot cor?l? (surtout les 2 APTs)



*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 

AT2G00680
AT5G54670 7
AT4G33650 7
AT1G30690 4
AT2G25180 7


*Gem2Net
Pas le lncRNA

AT2G25180, AT4G33650 ds meme cluster

**Non**, codings pas sexy, pas de lien entre eux...





*XLOC_002900/AT1G08687 lncRNA*
  ```{r, include = FALSE}
correlation_col_100_neg %>% filter(ncRNA == "XLOC_002900") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002900", "AT2G39800", "AT3G47690", "AT3G24240", "AT5G05560", "AT1G55540")
gene_list <- c("XLOC_002900", "AT2G39800", "AT3G24240") # Filtre pour select Correlation avec Ler aussi

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "XLOC_002900", gene %in% c("AT2G39800", "AT3G24240"))
```
Anticorr?lation chez Col et Ler avec 2 profils diff?rents




*Araport
novel transcribed region

AT2G39800	encodes a delta1-pyrroline-5-carboxylate synthase that catalyzes the rate-limiting enzyme in the biosynthesis of proline. Gene is expressed in reproductive organs and tissues under non-stress conditions but in the whole plant under water-limiting condition. Expression is also induced by abscisic acid and salt stress in a light-dependent manner. encodes a delta1-pyrroline-5-carboxylate synthase that catalyzes the rate-limiting enzyme in the biosynthesis of proline. Gene is expressed in reproductive organs and tissues under non-stress conditions but in the whole plant under water-limiting condition. Expression is also induced by abscisic acid and salt stress in a light-dependent manner. P5CS1 appears to be involved in salt stress responses related to proline accumulation, including protection from reactive oxidative species. P5CS1 appears to be present in different cells and/or different subcellular locations from P5CS2 in a tissue-dependent manner.  Source: TAIR, Jun 30, 2015.



*GV
Pas le lncRNA
Coding pas cor?l? entre eux

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 
XLOC_002900 8
AT2G39800   2
AT3G24240   2

*Gem2Net
PAs le lncRNA



**Non**, correlation pas si belle et coding pas fou




*XLOC_005603/AT3G03315 lncRNA*
  ```{r, include=FALSE}
correlation_col_100_neg %>% filter(ncRNA == "XLOC_005603") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_005603", "AT4G33650")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "XLOC_005603", gene %in% c("AT4G33650"))
```
Anticorr?lation chez Col et Ler similaire



*Araport
novel transcribed region



*GV
Pas le lncRNA

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 
Non

*Gem2Net
Pas le lncRNA


**Non**, coding pas fou et peu d'infos 





#Sur les comptage de Ler

## Positives
```{r}
correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 100) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame() 
correlation_ler_100 <- correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 100) %>% as.data.frame() 
```


*AT1G55525 lncRNA*
```{r, echo=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT1G55525") %>% select(gene) %>% as.list()
gene_list <- c("AT1G55525", "AT4G34410", "AT2G42430", "AT1G70940", "AT5G26930", "AT1G78240", "AT1G58340", "AT2G34650", "AT1G19220", "AT5G59030", "AT3G60630", "AT1G25220")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C7 <- correlation_all %>% filter(ncRNA == "AT1G55525", gene %in% c("AT4G34410", "AT2G42430", "AT1G70940", "AT5G26930", "AT1G78240", "AT1G58340", "AT2G34650", "AT1G19220", "AT5G59030", "AT3G60630", "AT1G25220")) %>% add_column(type = "Ler_count_Pos_100")
C7
```
Corr?lation chez Col et Ler similaire


*Araport
Unknown gene

AT1G19220	ARF19	auxin response factor 19
AT1G25220	ASB1	anthranilate synthase beta subunit 1
AT1G55525	AT1G55525	""
AT1G58340	ZF14	MATE efflux family protein
AT1G70940	PIN3	Auxin efflux carrier family protein
AT1G78240	TSD2	S-adenosyl-L-methionine-dependent methyltransferases superfamily protein
AT2G34650	PID	Protein kinase superfamily protein
AT2G42430	LBD16	lateral organ boundaries-domain 16
AT3G60630	HAM2	GRAS family transcription factor
AT4G34410	RRTF1	redox responsive transcription factor 1
AT5G26930	GATA23	GATA transcription factor 23
AT5G59030	COPT1	copper transporter 1

Go root

*GV
peu de donn?es avec lncRNA
L?g?re corr?lation

*clustering GLM
```{r, include=FALSE}
gene_list <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.data.frame()
``` 


AT4G34410 9
AT2G42430 6
AT1G70940 2
AT5G26930 9
AT1G78240 2
AT1G58340 6
AT2G34650 6
AT1G19220 2
AT5G59030 6
AT3G60630 2
AT1G25220 2
AT1G55525

6 et 9 se ressemble

*Gem2Net
Pas le lncRNA
AT1G19220, AT5G26930, AT2G42430
AT2G42430, AT4G34410
AT2G42430, AT2G34650
AT5G59030, AT1G78240
AT3G60630, AT5G59030


**Oui**, belle corr?lation avec des genes root sexy

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT1G55525", "AT4G34410", "AT2G42430", "AT1G70940", "AT5G26930", "AT1G78240", "AT1G58340", "AT2G34650", "AT1G19220", "AT5G59030", "AT3G60630", "AT1G25220")) %>% arrange(desc(type))
```
Pi resp
GO root



*AT1G66725 lncRNA*
*microRNA163*
```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT1G66725") %>% select(gene) %>% as.list()
gene_list <- c("AT1G66725", "AT2G26170", "AT1G08090", "AT3G53480", "AT1G10940")
gene_list <- c("AT1G66725", "AT2G26170") #Filtre

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT1G55525", gene %in% c("AT4G34410", "AT2G42430", "AT1G70940", "AT5G26930", "AT1G78240", "AT1G58340", "AT2G34650", "AT1G19220", "AT5G59030", "AT3G60630", "AT1G25220"))
```
Corr?lation positive chez Col et n?gative chez Ler


*Araport
Encodes a protein with similarity to thromboxane-A synthase, member of the CYP711A cytochrome P450 family. MAX1 is a specific repressor of vegetative axillary buds generated by the axillary meristem. Expressed in vascular traces in the rosette stem and axillary buds throughout plant development. Mutants have increased axillary branches. Along with MAX3,4 thought to mediate control of shoot branching via synthesis of a signal molecule which is transported over long distance mediated by MAX2. cDNA supports the existence of the longer transcript predicted for this locus, no cDNA isolated for shorter transcript. MAX1 downregulates 11 genes involved in flavonoid pathway (CHS, CHI, F3H, F3'H, FLS, DFR, ANS, UFGT, RT, AAC and GST).

MAX1 is implicated in synthesis of the carotenoid-derived branch regulator(s) from the root, it likely links long-distance signaling with local control of bud outgrowth

*GV
bof cor?l? (col, bien)

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 


AT1G66725 9
AT2G26170 5

*Gem2Net

**Non**, correlation pas si belle




*AT1G67195 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT1G67195") %>% select(gene) %>% as.list()
gene_list <- c("AT1G67195", "AT3G54870", "AT2G17500", "AT5G54670", "AT4G33650", "AT2G47000")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT1G55525", gene %in% c("AT4G34410", "AT2G42430", "AT1G70940", "AT5G26930", "AT1G78240", "AT1G58340", "AT2G34650", "AT1G19220", "AT5G59030", "AT3G60630", "AT1G25220"))
```
Corr?lation similaire chez Col et Ler


**Non**, MIR414 c un coding !
  
  
  
  
  *AT1G67195 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT2G31585") %>% select(gene) %>% as.list()
gene_list <- c("AT2G31585", "AT4G12550")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT4G12550"))
```
Corr?lation chez Col et ler avec profil diff?rent

**Non**, MIR414 c un coding !
  
  
  
  
  *AT2G32795 lncRNA*
  ```{r, include = FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT2G32795") %>% select(gene) %>% as.list()
gene_list <- c("AT2G32795", "AT5G58710")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT2G32795", gene %in% c("AT5G58710"))
```
Corr?lation Positive chez Ler et n?gative chez Col




*Araport
Unknown gene


*GV
Bof cor?l? (bien, col)

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 

AT2G32795
AT5G58710 7

*Gem2Net
inutile comme anticor?l? col


**Non**, coding pas sexy





*AT3G27884 lncRNA*
  ```{r, include = FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT3G27884") %>% select(gene) %>% as.list()
gene_list <- c("AT3G27884", "AT2G17500")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT3G27884", gene %in% c("AT2G17500"))
```
Corr?lation positive chez Ler et n?gative chez Col



*Araport
Unknown gene

coding = PILS5 is a member of the PIN-LIKES (PILS) protein family of putative auxin carriers residing at the endoplasmic reticulum. PILS5 activity decreases intracellular levels of free IAA by increasing the levels of certain auxin conjugates.

*GV

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 


*Gem2Net

**Non**, correlation col pas si belle et coding pas fou




*AT5G09710 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT5G09710") %>% select(gene) %>% as.list()
gene_list <- c("AT5G09710", "AT5G52050", "AT2G34650", "AT5G26930", "AT3G12280", "AT1G19220", "AT1G17110")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT3G27884", gene %in% c("AT2G17500"))
```

pseudogene :
  
  Magnesium transporter CorA-like family protein; BEST Arabidopsis thaliana protein match is: magnesium transporter 7 (TAIR:AT5G09690.2).  

*AT5G24105 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "AT5G24105") %>% select(gene) %>% as.list()
gene_list <- c("AT5G24105", "AT4G24670", "AT2G34680", "AT5G60660", "AT5G03150", "AT5G65800")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT5G24105", gene %in% c("AT5G24105", "AT4G24670", "AT2G34680", "AT5G60660", "AT5G03150", "AT5G65800"))
```
Corr?lation similaire entre Col et Ler
V?rifier que c pas codant

**Non**
  
  
  *Ler_NEW_RNA_F_12534 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "Ler_NEW_RNA_F_12534") %>% select(gene) %>% as.list()
gene_list <- c("Ler_NEW_RNA_F_12534", "AT5G51040", "AT1G71040")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "Ler_NEW_RNA_F_12534", gene %in% c("AT5G51040", "AT1G71040"))
```
Corr?lation chez Ler uniquement


*Araport

AT5G51040 = Encodes succinate dehydrogenase assembly factor 2 (SDHAF2), a low abundance mitochondrial protein needed for assembly and activity of mitochondrial complex II and for normal root elongation

AT1G71040 = Encodes LPR2. Function together with LPR1 (AT1G23010) and a P5-type ATPase (At5g23630/PDR2) in a common pathway that adjusts root meristem activity to Pi (inorganic phosphate) availability.

*GV

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 


*Gem2Net

**Non**, anticor?lation chez col aps clair




*Ler_NEW_RNA_R_34309 lncRNA*
  ```{r, echo=FALSE}
correlation_ler_100 %>% filter(ncRNA == "Ler_NEW_RNA_R_34309") %>% select(gene) %>% as.list()
gene_list <- c("Ler_NEW_RNA_R_34309", "AT4G05530", "AT1G08090")
gene_list <- c("Ler_NEW_RNA_R_34309", "AT4G05530") #Filtre

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, echo=FALSE}
C8 <- correlation_all %>% filter(ncRNA == "Ler_NEW_RNA_R_34309", gene %in% c("AT4G05530")) %>% add_column(type = "Ler_count_Pos_100")
C8
```
Corr?lation positive chez Ler uniquement


*Araport

Encodes a peroxisomal member of the short-chain dehydrogenase/reductase (SDR) family of enzymes. Loss of IBR1 function causes increased resistance to indole-3-butyric acid without affecting plant responses to IAA, NAA, and 2,4-D. This enzyme may be responsible for catalyzing a dehydrogenation step in the beta-oxidation-like conversion of IBA to IAA. 

SDRa is a peroxisomal oxidoreductase-like protein that is required for response to pro-auxins. [SDRa]

*GV

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 

AT4G05530 7

*Gem2Net

**WhyNot**, pur hasard...


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("Ler_NEW_RNA_R_34309", "AT4G05530")) %>% arrange(desc(type))
```
Col_ler diff


*XLOC_000851/AT1G47570 lncRNA*
  ```{r, include=FALSE}
correlation_ler_100 %>% filter(ncRNA == "XLOC_000851") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_000851", "AT5G48630")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "XLOC_000851", gene %in% c("AT5G48630"))
```
Corr?lation Positive chez Ler et n?gative chez Col

*Araport
RING/U-box superfamily protein; FUNCTIONS IN: zinc ion binding; EXPRESSED IN: 21 plant structures; EXPRESSED DURING: 10 growth stages; CONTAINS InterPro DOMAIN/s: Zinc finger, RING-type (InterPro:IPR001841), Zinc finger, C3HC4 RING-type (InterPro:IPR018957
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, surment une proteine
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_001547 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_100 %>% filter(ncRNA == "XLOC_001547") %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              gene_list <- c("XLOC_001547", "AT1G12560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C9 <- correlation_all %>% filter(ncRNA == "XLOC_001547", gene %in% c("AT1G12560")) %>% add_column(type = "Ler_count_Pos_100")
                                                                                                                                                                                                                                              C9
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation chez Ler uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              coding = ? 	Member of Alpha-Expansin Gene Family. Naming convention from the Expansin Working Group (Kende et al, 2004. Plant Mol Bio). Containing a conserved root hair-specific cis-element RHE. Expressed specifically in root hair cell and involved in root hair elongation
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G12560 5
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **WhyNot**, pur hasard
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_001547", "AT1G12560")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              precurseur smRNA24 Ler>Col
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002093/AT1G08925 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_100 %>% filter(ncRNA == "XLOC_002093") %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002093", "AT1G17110", "AT5G52050", "AT1G19220", "AT5G10720")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C10 <- correlation_all %>% filter(ncRNA == "XLOC_002093", gene %in% c("AT1G17110", "AT5G52050", "AT1G19220", "AT5G10720")) %>% add_column(type = "Ler_count_Pos_100")
                                                                                                                                                                                                                                              C10
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Correlation similaire chez Col et Ler
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              novel transcribed region
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G17110 = Encodes a ubiquitin-specific protease, and its activity has been confirmed in an in vitro assay. ubp15 mutants have defects in cell proliferation, and the associated processes of leaf, root, stem, flower, and silique development. UBP15 can be found in the nucleus and cytoplasm in transient assays. Though UBP15 is expressed in many tissues, UBP15 transcript levels are higher in rosette leaves and inflorescences than in other parts of the plant
                                                                                                                                                                                                                                              UBP15 affects leaf shape by controlling cell proliferation, and has ubiquitin-specific protease activity in vitro.
                                                                                                                                                                                                                                              Genetic analyses indicate that UBP15 functions antagonistically in a common pathway with DA1 to influence seed size, but does so independently of DA2 and EOD1.
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT5G52050 = AtDTX50, functions as an abscissic acid efflux transporter. [AtDTX50]
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G19220 = ARF19,	Encodes an auxin response factor that contains the conserved VP1-B3 DNA-binding domain at its N-terminus and the Aux/IAA-like domains III and IV present in most ARFs at its C-terminus. The protein interacts with IAA1 (yeast two hybrid) and other auxin response elements such as ER7 and ER9 (yeast one hybrid). ARF19 protein can complement many aspects of the arf7 mutant phenotype and , together with ARF7, is involved in the response to ethylene. In the arf7 arf19 double mutant, several auxin-responsive genes (e.g. IAA5, LBD16, LBD29 and LBD33) are no longer upregulated by auxin.
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT5G10720, ROS..
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              GO root (les 4)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              Plutot cor?l?
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              Not lncRNA
                                                                                                                                                                                                                                              AT1G17110 10
                                                                                                                                                                                                                                              AT5G52050 9
                                                                                                                                                                                                                                              AT1G19220 2
                                                                                                                                                                                                                                              AT5G10720 2
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              (2 et 10 se ressemble)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              AT1G17110, AT5G10720
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Oui**, belle corr?lation avec gene int?ressant
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002093", "AT1G17110", "AT5G52050", "AT1G19220", "AT5G10720")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              GO root
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002900/AT1G08687 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_100 %>% filter(ncRNA == "XLOC_002900") %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002900", "AT4G01060", "AT4G39990", "AT3G04630")
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002900", "AT3G04630")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C11 <- correlation_all %>% filter(ncRNA == "XLOC_002900", gene %in% c("AT3G04630")) %>% add_column(type = "Ler_count_Pos_100")
                                                                                                                                                                                                                                              C11
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation similaire entre Col et Ler mais profil diff?rents
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              novel transcribed region
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT3G04630 = Member of a small gene family which have a KLEEK domain which may be involved in protein- protein interactions. Over expression of WDL1 results in abnormal root development
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              XLOC_002900 8
                                                                                                                                                                                                                                              AT3G04630 4
                                                                                                                                                                                                                                              (ressemble surtout pr comptage ler)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Oui**, corr?lation avec un coding int?ressant
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002900", "AT3G04630")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_005603/AT3G03315 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_100 %>% filter(ncRNA == "XLOC_005603") %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              gene_list <- c("XLOC_005603", "AT3G11260", "AT4G01540", "AT1G14740", "AT3G03660", "AT5G52310", "AT3G07390")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C12 <- correlation_all %>% filter(ncRNA == "XLOC_005603", gene %in% c("AT3G11260", "AT4G01540", "AT1G14740", "AT3G03660", "AT5G52310", "AT3G07390")) %>% add_column(type = "Ler_count_Pos_100")
                                                                                                                                                                                                                                              C12
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Correlation similaire chez Col et Ler
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              novel transcribed region
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G14740	TTA1	class I heat shock protein, putative (DUF1423)
                                                                                                                                                                                                                                              AT3G03660	WOX11	WUSCHEL related homeobox 11
                                                                                                                                                                                                                                              AT3G07390	AIR12	auxin-induced in root cultures-like protein
                                                                                                                                                                                                                                              AT3G11260	WOX5	WUSCHEL related homeobox 5
                                                                                                                                                                                                                                              AT4G01540	NTM1	NAC with transmembrane motif1
                                                                                                                                                                                                                                              AT5G52310	LTI78	low-temperature-responsive protein 78 (LTI78) / desiccation-responsive protein 29A (RD29A)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              GO root dvlpmnt
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              Bof cor?l?
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G14740 2
                                                                                                                                                                                                                                              AT3G03660 9
                                                                                                                                                                                                                                              AT3G07390 6
                                                                                                                                                                                                                                              AT3G11260 9
                                                                                                                                                                                                                                              AT4G01540 6
                                                                                                                                                                                                                                              AT5G52310 9
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              9, 6 se ressemble
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              AT3G03660, AT4G01540
                                                                                                                                                                                                                                              AT3G07390, AT5G52310
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Bof**, dur de trouver un lien entre Col et Ler expression du lncRNa / a ses potentiels targets
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_005603", "AT3G11260", "AT4G01540", "AT1G14740", "AT3G03660", "AT5G52310", "AT3G07390")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              GO root (les 2 gem2net)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              # Negatives
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_Ler %>% arrange(pearson_corr_ler) %>% filter(row_number() <= 100) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              correlation_ler_100_neg
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT1G55525 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo = FALSE}
                                                                                                                                                                                                                                              correlation_ler_100_neg %>% filter(ncRNA == "AT1G55525") %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              gene_list <- c("AT1G55525", "AT2G01830", "AT5G35750")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo = FALSE}
                                                                                                                                                                                                                                              C13 <- correlation_all %>% filter(ncRNA == "AT1G55525", gene %in% c("AT2G01830", "AT5G35750")) %>% add_column(type = "Ler_count_Neg_100")
                                                                                                                                                                                                                                              C13
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation negative chez Col et Ler
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              Unknown gene
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT2G01830	HK4/WOL	CHASE domain containing histidine kinase protein
                                                                                                                                                                                                                                              Histidine kinase: cytokinin-binding receptor that transduces cytokinin signals across the plasma membrane
                                                                                                                                                                                                                                              CK receptor
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT5G35750	HK2	histidine kinase 2
                                                                                                                                                                                                                                              Encodes histidine kinase AHK2.
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Plutot anticor?l? (peu de donn?es...)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT2G01830 4
                                                                                                                                                                                                                                              AT5G35750 4
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas coding, pas de groupe/2
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Bof**, Pas plus d'infos que correlation...

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT1G55525", "AT2G01830", "AT5G35750")) %>% arrange(desc(type))
```
Pi resp


*AT1G66725/microRNA163 lncRNA*
```{r, include = FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT1G66725") %>% select(gene) %>% as.list()
gene_list <- c("AT1G66725", "AT5G65700")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "XLOC_005603", gene %in% c("AT3G11260", "AT4G01540", "AT1G14740", "AT3G03660", "AT5G52310", "AT3G07390"))
```
Corr?lation n?gative chez Co let Ler

**Non**, mir163



*AT2G31585 lncRNA*
```{r, echo=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT2G31585") %>% select(gene) %>% as.list()
gene_list <- c("AT2G31585", "AT2G37630", "AT4G17490")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C14 <- correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))  %>% add_column(type = "Ler_count_Neg_100")
C14
```
Corr?lation chez Col et Ler avec profil diff?rent

*Araport
Unknown gene

AT2G37630 = Encodes a MYB-domain protein involved in specification of the leaf proximodistal axis. Mutation results in lobed and dissected leaves with a characteristic asymmetry. Homologous to the Antirrhinum PHANTASTICA (PHAN) and maize ROUGH SHEATH2 (RS2) genes Asymmetric placement of auxin response at the distal leaf tip precedes visible asymmetric leaf growth. Acts alongside AXR1 to exclude BP expression in leaves and with PIN1 to repress BP and promote lateral organ growth. Interacts physically with AS2 to form a complex that binds to the BP promoter and silences BP. Also functions as a regulator of the plant immune response
RDR6, SGS3 and AGO7 act in the same pathway, which genetically interacts with the AS1-AS2 pathway for leaf development.
Our data suggest that RS2/AS1 and HIRA mediate the epigenetic silencing of knox genes, possibly by modulating chromatin structure.

AT4G17490 = Encodes a member of the ERF (ethylene response factor) subfamily B-3 of ERF/AP2 transcription factor family (ATERF-6). The protein contains one AP2 domain. There are 18 members in this subfamily including ATERF-1, ATERF-2, AND ATERF-5. It is involved in the response to reactive oxygen species and light stress



*GV
Plutot anticor?l? (surtout avec At2G)


*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 
AT2G31585 2
AT2G37630 7
AT4G17490 9

(2 et 7 semble inverse)

*Gem2Net
Pas le lncRNA, les 2 genes pas ensemble

**Oui**, histoire de chromatin int?ressante avec AT2G, bien anticor?l? avec lui (v?rifier que histoire pas completement ?lucid?...)

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT2G31585", "AT2G37630", "AT4G17490")) %>% arrange(desc(type))
```
Col_ler diff


*AT3G52670 lncRNA*
```{r, include=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT3G52670") %>% select(gene) %>% as.list()
gene_list <- c("AT3G52670", "AT1G78870")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))
```
Corr?lation chez Col et Ler avec profil diff?rent


*Araport
Pseudogene
FBD%2C F-box%2C Skp2-like and Leucine Rich Repeat domains containing protein  
risque que sa soit pseudogene ou protein...

**Non**,



*AT3G61198 lncRNA BAP1*
```{r, include=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT3G61198") %>% select(gene) %>% as.list()
gene_list <- c("AT3G61198", "AT5G63770")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))
```
Corr?lation chez Col et Ler

BAP1 lncRNA

**Non**


*AT5G09710 lncRNA*
```{r, include = FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT5G09710") %>% select(gene) %>% as.list()
gene_list <- c("AT5G09710", "AT5G57740", "AT1G22530", "AT5G35750", "AT1G72160", "AT5G48000")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))
```
Corr?lation chez Col et Ler

Pseudogene, 
 	Magnesium transporter CorA-like family protein
Risque que sa soit pseudogene ou une proteine

**Non**,





*AT5G24105 lncRNA*
```{r, include = FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "AT5G24105") %>% select(gene) %>% as.list()
gene_list <- c("AT5G24105", "AT1G07880", "AT1G27740", "AT3G07390", "AT3G11260", "AT5G52310", "AT4G09510", "AT3G03660", "AT1G12950")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include = FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))
```
Corr?lation chez Col et Ler
Attention AGP..

**Non**

*Ler_NEW_RNA_F_12534/AT2G27440 lncRNA*
```{r, include = FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "Ler_NEW_RNA_F_12534") %>% select(gene) %>% as.list()
gene_list <- c("Ler_NEW_RNA_F_12534", "AT5G05560")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "AT2G31585", gene %in% c("AT2G37630", "AT4G17490"))
```
Corr?lation n?gative chez Ler uniquement


*Araport
 	pseudogene of Rho GTPase activating protein with PAK-box/P21-Rho-binding domain-containing protei
 	
**Non**



*Ler_NEW_RNA_R_14094/AT2G08790 lncRNA*
```{r, echo=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "Ler_NEW_RNA_R_14094") %>% select(gene) %>% as.list()
gene_list <- c("Ler_NEW_RNA_R_14094", "AT2G33790")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C15 <- correlation_all %>% filter(ncRNA == "Ler_NEW_RNA_R_14094", gene %in% c("AT2G33790")) %>% add_column(type = "Ler_count_Neg_100")
C15
```
Corr?lation chez Ler uniquement

*Araport
Gene

AT2G33790 = pollen Ole e 1 allergen protein containing 14.6% proline residues, similar to arabinogalactan protein (Daucus carota) GI:11322245, SP:Q03211 Pistil-specific extensin-like protein precursor (PELP) Nicotiana tabacum; contains Pfam profile PF01190: Pollen proteins Ole e I family 
role in root generation (cell wall expression)

*GV
Pas le lncRNA

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 
AT2G33790 2

*Gem2Net

**Oui**, coding int?resant, belle anticorr?lation chez ler uniquement

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("Ler_NEW_RNA_R_14094", "AT2G33790")) %>% arrange(desc(type))
```
Col_ler diff
GO root


*XLOC_001547 lncRNA*
```{r}
correlation_ler_100_neg %>% filter(ncRNA == "XLOC_001547") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_001547", "AT5G48600")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r}
correlation_all %>% filter(ncRNA == "XLOC_001547", gene %in% c("AT5G48600"))
```
Corr?lation chez Ler uniquement




*Araport
AT5G48600 = structural maintenance of chromosome 3
results reveal a critical role for AtCAP-C during cell division

**Non**, coding aps sexy



*XLOC_002093/AT1G08925 lncRNA*
```{r, echo=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "XLOC_002093") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_002093", "AT5G35750", "AT5G03150", "AT1G30690", "AT1G72160", "AT1G22530", "AT1G72150", "AT5G48000")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C16 <- correlation_all %>% filter(ncRNA == "XLOC_002093", gene %in% c("AT5G35750", "AT5G03150", "AT1G30690", "AT1G72160", "AT1G22530", "AT1G72150", "AT5G48000")) %>% add_column(type = "Ler_count_Neg_100")
C16
```
Corr?lation chez Col et Ler



*Araport
novel transcribed region

AT1G22530	PATL2	PATELLIN 2
AT1G30690	AT1G30690	Sec14p-like phosphatidylinositol transfer family protein
AT1G72150	PATL1	PATELLIN 1
AT1G72160	AT1G72160	Sec14p-like phosphatidylinositol transfer family protein
AT5G03150	JKD	C2H2-like zinc finger protein
AT5G35750	HK2	histidine kinase 2
AT5G48000	CYP708A2	cytochrome P450, family 708, subfamily A, polypeptide 2

GO cell division


AT1G72150 = PATL1 is involved in membrane-trafficking events associated with cell-plate expansion or maturation and point to the involvement of phosphoinositides in cell-plate biogenesis. [Patellin1 (PATL1)]

AT1G72160 = the interaction between Alfalfa mosaic virus (AMV) MP and Arabidopsis Patellin 3 (atPATL3) and Patellin 6 (atPATL6), is reported.


*GV
Pas le lncRNA

Plutot cor?l?

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 
XLOC_002093
AT5G35750 4
AT5G03150 7
AT1G30690 4
AT1G72160 7
AT1G22530 4
AT1G72150 5
AT5G48000 4

se ressemble bof

*Gem2Net
Pas le lncRNA
AT1G22530, AT1G30690, AT1G72150, AT1G72160
dans plusieurs clusters

**Bof**, belle antir?gulation mais pas de role ds root clair des codings


``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002093", "AT5G35750", "AT5G03150", "AT1G30690", "AT1G72160", "AT1G22530", "AT1G72150", "AT5G48000")) %>% arrange(desc(type))
```
Pi resp
GO root


*XLOC_005603/AT3G03315 lncNRA*
```{r, echo=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "XLOC_005603") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_005603", "AT5G40330", "AT1G77740")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C17 <- correlation_all %>% filter(ncRNA == "XLOC_005603", gene %in% c("AT5G40330", "AT1G77740")) %>% add_column(type = "Ler_count_Neg_100")
C17
```
Corr?lation chez Col et Ler


*Araport
novel transcribed region

AT5G40330 = 	Encodes a MYB gene that, when overexpressed ectopically, can induce ectopic trichome formation. It is a member of subgroup 15, together with WER and GL1. Members of this subgroup share a conserved motif of 19 amino acids in the putative transcription activation domain at the C-terminal end. The gene is expressed in leaves, stems, flowers, seeds and roots and quite strongly in trichomes. There is partial functional redundancy between ATMYB23 and GL1. The two proteins are functionally equivalent with respect to the regulation of trichome initiation but not with respect to trichome branching - which is controlled by MYB23 and not GL1.
MYB23 participates in a positive feedback loop to reinforce cell fate decisions and ensure robust establishment of the cell type pattern in the Arabidopsis root epidermis.

AT1G77740 = 	Encodes PIP5K2, a phosphatidylinositol-4-phosphate 5-kinase (PtdIns(4)P 5-kinase 2; or PIP5K2 that is involved in regulating lateral root formation and root gravity response. 
PIP5K2 is involved in regulating lateral root formation and root gravity response, and reveal a critical role of PIP5K2/PtdIns(4,5)P(2) in root development through regulation of PIN proteins, providing direct evidence of crosstalk between the phosphatidylinositol signaling pathway and auxin response, and new insights into the control of polar auxin transport

*GV
Pas le lncRNA

Coding non cor?l?


*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 
AT1G77740 7
AT5G40330 4

se ressemble bof

*Gem2Net
Pas le lncRNA
coding jamais ds 2 meme cluster


**Bof**, rien de plus a part correlation

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_005603", "AT5G40330", "AT1G77740")) %>% arrange(desc(type))
```
Pi resp


*XLOC_005937/AT3G00660 lncRNA*
```{r, include=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "XLOC_005937") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_005937", "AT5G48630")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, include=FALSE}
correlation_all %>% filter(ncRNA == "XLOC_005937", gene %in% c("AT5G48630"))
```
Corr?lation positive chez Col et n?gative chez Ler


*Araport
novel transcribed region

AT5G48630 = Cyclin family protein  


**Non**, coding pas sexy et correlation pas si belle






*XLOC_011131/AT5G09295 lncRNA*
```{r, echo=FALSE}
correlation_ler_100_neg %>% filter(ncRNA == "XLOC_011131") %>% select(gene) %>% as.list()
gene_list <- c("XLOC_011131", "AT5G20490", "AT3G18780")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
  
``` 
```{r, echo=FALSE}
C18 <- correlation_all %>% filter(ncRNA == "XLOC_011131", gene %in% c("AT5G20490", "AT3G18780")) %>% add_column(type = "Ler_count_Neg_100")
C18
```
Corr?lation chez Col et Ler, profil diff?rent


*Araport
Gene

AT5G20490 = Encodes a member of the type XI myosin protein family involved in root hair growth, trichome development, and organelle trafficking. Required for fast root hair growth. This gene appears to be expressed at low levels throughout the plant. 

AT3G18780 = Encodes an actin that is constitutively expressed in vegetative structures but not pollen. ACT2 is involved in tip growth of root hairs


*GV
Pas le lncRNA
correlatin des coding

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 
XLOC_011131 2
AT5G20490 7
AT3G18780 7

2 inverse de 7

*Gem2Net
Pas le lncRNA
Pas de cluster des codings


**Oui**, belle anticor?lation avec profil different

``` {r, echo=FALSE}
all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_011131", "AT5G20490", "AT3G18780")) %>% arrange(desc(type))
```
Col_ler diff


# Correlation diff?rentiel entre ?cotypes

```{r}
nrow(pearson_correlation_df)*0.05
```
Sur les 1000 1eres corr?lations les plus fortes = (environ) 5 % 
 
## 1000 1eres corr?lation positives des comptages Col
```{r}
correlation_Col %>% arrange(desc(pearson_corr_col)) %>% filter(row_number() <= 1000) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
corr_Pos_1000 <- correlation_Col %>% 
  arrange(desc(pearson_corr_col)) %>% 
  filter(row_number() <= 1000) %>% 
  inner_join(correlation_Ler) %>% 
  mutate(pearson_diff = pearson_corr_col-pearson_corr_ler) 
```
Parims ces corr?lations, recherche des 10 couples les plus diff?rentiels entre Col et Ler
```{r, echo=FALSE}
corr_Pos_1000 %>% arrange(desc(pearson_diff)) %>% head(10) #AT1G52347 source de faux-positif
corr_Pos_1000 %>% filter(ncRNA != "AT1G52347") %>% arrange(desc(pearson_diff)) %>% head(10)
```

*XLOC_007825/AT4G04595 lncRNA*
```{r, include=FALSE}
gene_list <- c("XLOC_007825", "AT5G48000", "AT5G48010", "AT5G47990", "At5g47980", "Ler_NEW_RNA_F_34582")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")

```
Corr?lation positive chez Col et n?gative chez Ler
(sortie juste avec AT5g48000)

*Araport
novel transcribed region

AT5g48000 = Encodes a member of the CYP708A family of cytochrome P450 enzymes. THAH appears to add a hydroxyl group to the triterpene thalianol. thah1 mutants have an elevated accumulation of thalianol. thah1-1 mutants have longer roots than wild type plants. Thalian-diol and desaturated thalian-diol are lost from the root extracts of thah1-1 mutants. Overexpression of the sequence from At5g48000.1 rescues the thah1-1 mutant phenotype (Field 2008); it is unknown whether the shorter sequences associated with other gene models would provide functional complementation

Appartient au thalianol cluster :

Data indicate that the thalianol gene cluster consists of four contiguous coexpressed genes, At5g48010, At5g48000, At5g47990, and At5g47980, encoding enzymes that together catalyze the synthesis and elaboration of thalianol.


Pr?sence d'un lncRNA dans ce cluster :  AT5G07035 = Ler_NEW_RNA_F_34582, GO root dev dans notre RNAseq
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              XLOC_007825 8
                                                                                                                                                                                                                                              AT5g48000   4
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              cot? ler p^rofil semblable
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, thalianol cluster risqu? et anticor?lation ler pas si belle
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G13495 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G13495", "AT1G23010")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation positive chez Col et n?gative chez Ler
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C19 <- correlation_all %>% filter(ncRNA == "AT4G13495", gene %in% c("AT1G23010")) %>% add_column(type = "Col_count_Pos_1000_diff")
                                                                                                                                                                                                                                              C19
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              Unknown gene
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G23010 = LPR1 Encodes a protein with multicopper oxidase activity. Located in ER. Function together with LPR2 (AT1G71040) and a P5-type ATPase (At5g23630/PDR2) in a common pathway that adjusts root meristem activity to Pi (inorganic phosphate) availability.
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Plutot anticor?l? ecotype assay (bien)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              AT4G13495 5
                                                                                                                                                                                                                                              AT1G23010 6
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              se ressemble pas
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Oui**, induction chez col positive pr expression LPR1 logique avc reponse differentiel ecotyp a -Pi
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT4G13495", "AT1G23010")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              smRNA 24 chez Col
                                                                                                                                                                                                                                              GO Pi starvation
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G13495 lncRNA*
                                                                                                                                                                                                                                                ```{r, include = FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G13495", "AT5G64630")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation positive chez Col et n?gative chez Ler
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              FAS2 = Chromatin Assembly Factor-1 (CAF-1) p60 subunit. Involved in organization of the shoot and root apical meristems. In Arabidopsis, the three CAF-1 subunits are encoded by FAS1, FAS2 and, most likely, MSI1, respectively. Mutations in FAS1 or FAS2 lead to increased frequency of homologous recombination and T-DNA integration in Arabidopsis.
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, coding pas sexy 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ### Les 20 meilleures corr?lations Col avec ncRNA non exprim? chez Ler (NA ler)
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              corr_Pos_1000 %>% 
                                                                                                                                                                                                                                                mutate(is.na(pearson_corr_ler)) %>% 
                                                                                                                                                                                                                                                filter(is.na(pearson_corr_ler) == "TRUE") %>% 
                                                                                                                                                                                                                                                select(ncRNA, gene, pearson_corr_col) %>% 
                                                                                                                                                                                                                                                arrange(desc(pearson_corr_col)) %>% 
                                                                                                                                                                                                                                                head(20)
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Col_NEW_RNA_F_30071 lncRNA*
                                                                                                                                                                                                                                                ```{r, include = FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Col_NEW_RNA_F_30071", "AT1G52150")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation chez Col uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, correlation pas si belle
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G34881 lncNRA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G34881", "AT1G09560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation chez Col uniquement
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C20 <- correlation_all %>% filter(ncRNA == "AT4G34881", gene %in% c("AT1G09560")) %>% add_column(type = "Col_count_Pos_1000_spe")
                                                                                                                                                                                                                                              C20
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              transmembrane protein
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT1G09560 = Encodes a plasodesmata-located protein involved in regulating primary root growth by controlling phloem-mediated allocation of resources between the primary and lateral root meristems
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Pas de correlation assay ecotype (pas grave)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              AT4G34881 3
                                                                                                                                                                                                                                              AT1G09560 5
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              se ressemble
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas le lncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Bof**, pur hasard 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT4G34881", "AT1G09560")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              GO root dev
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ### 1000 1eres corr?lation negatives des comptages Col
                                                                                                                                                                                                                                              ```{r}
                                                                                                                                                                                                                                              corr_Neg_1000 <- correlation_Col %>% 
                                                                                                                                                                                                                                                arrange(pearson_corr_col) %>% 
                                                                                                                                                                                                                                                filter(row_number() <= 1000) %>% 
                                                                                                                                                                                                                                                inner_join(correlation_Ler) %>% 
                                                                                                                                                                                                                                                mutate(pearson_diff = pearson_corr_col-pearson_corr_ler) 
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Parims ces corr?lations, recherche des 10 couples les plus diff?rentiels entre Col et Ler
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              corr_Neg_1000 %>% arrange(pearson_diff) %>% head(10) 
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT2G32795 lncRNA*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT2G32795", "AT2G33790")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation positive chez ler, negative chez col
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C21 <- correlation_all %>% filter(ncRNA == "AT2G32795", gene %in% c("AT2G33790")) %>% add_column(type = "Col_count_Neg_1000_diff")
                                                                                                                                                                                                                                              C21
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              Unknown gene
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT2G33790 = pollen Ole e 1 allergen protein containing 14.6% proline residues, similar to arabinogalactan protein (Daucus carota) GI:11322245, SP:Q03211 Pistil-specific extensin-like protein precursor (PELP) {Nicotiana tabacum}; contains Pfam profile PF01190: Pollen proteins Ole e I family 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Plutot anticor?l? (bien)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              AT2G33790 2
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              Pas retrouv? ds meme cluster (bien)
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Oui**, anticor?lation soutenu/GV
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT2G32795", "AT2G33790")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              GO root dev
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002384/At1NC041650 lncRNA*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002384", "AT3G01090")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation negative chez col uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT3G01090 = encodes a SNF1-related protein kinase that physically interacts with SCF subunit SKP1/ASK1 and 20S proteosome subunit PAD1. It can also interact with PRL1 DWD-containing protein. Based on in vitro degradation assays and cul4cs and prl1 mutants, there is evidence that AKIN10 is degraded in a proteasome-dependent manner, and that this depends on a CUL4-PRL1 E3 ligase
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, absence lncRNa modifie pas tellement coprealtion coding
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT1G49032 lncRNA*
                                                                                                                                                                                                                                                ```{r, include = FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT1G49032", "AT1G80600")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation n?gative chez col uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ***Non**, correlation moche et risque protein
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ### Les 20 meilleures anticor?lations Col avec ncRNA non exprim? chez Ler (NA ler)
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              corr_Neg_1000 %>% 
                                                                                                                                                                                                                                                mutate(is.na(pearson_corr_ler)) %>% 
                                                                                                                                                                                                                                                filter(is.na(pearson_corr_ler) == "TRUE") %>% 
                                                                                                                                                                                                                                                select(ncRNA, gene, pearson_corr_col) %>% 
                                                                                                                                                                                                                                                arrange(pearson_corr_col) %>% 
                                                                                                                                                                                                                                                head(20)
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G34881 lncRNA*
                                                                                                                                                                                                                                                ```{r, include = FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G34881", "AT5G05560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation negative chez col uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              transmembrane protein
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT5G05560 = Encodes a subunit of the Arabidopsis thaliana E3 ubiquitin ligase complex that plays a synergistic role with APC4 both in female gametogenesis and in embryogenesis.  
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              Bof anticor?l?
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              AT4G34881 3
                                                                                                                                                                                                                                              AT5G05560 10
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              plutot diff?rent
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, coding pas sexy
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Col_NEW_RNA_F_30663 lncRNA*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Col_NEW_RNA_F_30663", "AT3G47690")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?altion negative chez col uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT3G47690 = encodes a homolog of animal microtubule-end-binding protein. There are two other members of this family. EB1 forms foci at regions where the minus ends of microtubules are gathered during mitosis and early cytokinesi
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, anticor?lation pas flagrante et coding pas si sexy
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ### Les 20 meilleures corr?lations Ler avec ncRNA non exprim? chez Col (NA ler)
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_1000 <- correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 1000) %>% as.data.frame() %>% inner_join(correlation_Col)
                                                                                                                                                                                                                                              correlation_ler_1000 %>% 
                                                                                                                                                                                                                                                mutate(is.na(pearson_corr_col)) %>% 
                                                                                                                                                                                                                                                filter(is.na(pearson_corr_col) == "TRUE") %>% 
                                                                                                                                                                                                                                                select(ncRNA, gene, pearson_corr_ler) %>% 
                                                                                                                                                                                                                                                arrange(desc(pearson_corr_ler)) %>% 
                                                                                                                                                                                                                                                head(20)
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Ler_NEW_RNA_R_34309 lncRNA*
                                                                                                                                                                                                                                                ```{r, include= FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Ler_NEW_RNA_R_34309", "AT4G05530")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Correlation positive chez ler, n?gative chez col
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              AT4G05530 = Encodes a peroxisomal member of the short-chain dehydrogenase/reductase (SDR) family of enzymes. Loss of IBR1 function causes increased resistance to indole-3-butyric acid without affecting plant responses to IAA, NAA, and 2,4-D. This enzyme may be responsible for catalyzing a dehydrogenation step in the beta-oxidation-like conversion of IBA to IAA
                                                                                                                                                                                                                                              SDRa is a peroxisomal oxidoreductase-like protein that is required for response to pro-auxins. [SDRa]
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, trop risqu?, pas d'AGI et ler express et anticor?lation pas si marqu?




*XLOC_001547 / At1NC109950 lncRNA*
```{r, include=FALSE}
gene_list <- c("XLOC_001547", "AT1G12560")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 
stat_normcounts_current <- stat_normcounts_tidy %>%
  filter(ID %in% gene_list) 

ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts, group=ecotype)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free")
```
Corr?lation positive chez ler unqiuement


*Araport
AT1G12560 = Member of Alpha-Expansin Gene Family. Naming convention from the Expansin Working Group (Kende et al, 2004. Plant Mol Bio). Containing a conserved root hair-specific cis-element RHE. Expressed specifically in root hair cell and involved in root hair elongation

*GV

*clustering GLM
```{r, include=FALSE}
gene_list <- c() %>% as.data.frame()
``` 


*Gem2Net

**Non**, trop risqu?, pas d'AGI et ler express et anticor?lation pas si marqu?
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                ### Les 20 meilleures corr?lations negatives Ler avec ncRNA non exprim? chez Col (NA ler)
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              correlation_ler_1000_neg <- correlation_Ler %>% arrange(pearson_corr_ler) %>% filter(row_number() <= 1000) %>% as.data.frame() %>% inner_join(correlation_Col)
                                                                                                                                                                                                                                              correlation_ler_1000_neg %>% 
                                                                                                                                                                                                                                                mutate(is.na(pearson_corr_col)) %>% 
                                                                                                                                                                                                                                                filter(is.na(pearson_corr_col) == "TRUE") %>% 
                                                                                                                                                                                                                                                select(ncRNA, gene, pearson_corr_ler) %>% 
                                                                                                                                                                                                                                                arrange(pearson_corr_ler) %>% 
                                                                                                                                                                                                                                                head(20)
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_001547/At1NC109950 lncRNA*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_001547", "AT5G48600")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              corr?lation negative chez ler uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              structural maintenance of chromosome 3
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**, coding pas si sexy
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Ler_NEW_RNA_F_12534/AT2G27440 lncRNA*
                                                                                                                                                                                                                                                ```{r}
                                                                                                                                                                                                                                              gene_list <- c("Ler_NEW_RNA_F_12534", "AT5G05560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation negative chez ler uniquement
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              pseudogene of rac GTPase activating protein
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Non**
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *Ler_NEW_RNA_R_14094/AT2G08790*
                                                                                                                                                                                                                                                ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Ler_NEW_RNA_R_14094", "AT2G33790")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Corr?lation chez ler unqiuement
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              C22 <- correlation_all %>% filter(ncRNA == "Ler_NEW_RNA_R_14094", gene %in% c("AT2G33790")) %>% add_column(type = "Col_count_Neg_1000_spe")
                                                                                                                                                                                                                                              C22
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Araport
                                                                                                                                                                                                                                              Gene
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AT2G33790 = pollen Ole e 1 allergen protein containing 14.6% proline residues, similar to arabinogalactan protein (Daucus carota) GI:11322245, SP:Q03211 Pistil-specific extensin-like protein precursor (PELP) {Nicotiana tabacum}; contains Pfam profile PF01190: Pollen proteins Ole e I family
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *GV
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *clustering GLM
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              gene_list <- c() %>% as.data.frame()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Gem2Net
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              **Bof**, risqu?
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                ``` {r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("Ler_NEW_RNA_R_14094", "AT2G33790")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              GO root dev
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ## 1000 1eres corr?lation positives des comptages Ler
                                                                                                                                                                                                                                              ```{r}
                                                                                                                                                                                                                                              correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 1000) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
                                                                                                                                                                                                                                              corr_Pos_1000_ler <- correlation_Ler %>% 
                                                                                                                                                                                                                                                arrange(desc(pearson_corr_ler)) %>% 
                                                                                                                                                                                                                                                filter(row_number() <= 1000) %>% 
                                                                                                                                                                                                                                                inner_join(correlation_Col) %>% 
                                                                                                                                                                                                                                                mutate(pearson_diff = pearson_corr_col-pearson_corr_ler) 
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Parims ces corr?lations, recherche des 10 couples les plus diff?rentiels entre Col et Ler
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              corr_Pos_1000_ler %>% arrange(desc(pearson_diff)) %>% head(10) %>% filter(ncRNA != "AT5G09710") #c une proteine potentiellement
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              Rien...
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ## 1000 1eres corr?lation negatives des comptages Ler
                                                                                                                                                                                                                                              ```{r}
                                                                                                                                                                                                                                              correlation_Ler %>% arrange(desc(pearson_corr_ler)) %>% filter(row_number() <= 1000) %>% group_by(ncRNA) %>% summarise(n=n()) %>% as.data.frame()
                                                                                                                                                                                                                                              corr_Neg_1000_ler <- correlation_Ler %>% 
                                                                                                                                                                                                                                                arrange(pearson_corr_ler) %>% 
                                                                                                                                                                                                                                                filter(row_number() <= 1000) %>% 
                                                                                                                                                                                                                                                inner_join(correlation_Col) %>% 
                                                                                                                                                                                                                                                mutate(pearson_diff = pearson_corr_col-pearson_corr_ler) 
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              Parims ces corr?lations, recherche des 10 couples les plus diff?rentiels entre Col et Ler
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              corr_Neg_1000_ler %>% arrange(desc(pearson_diff)) %>% head(10) 
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              Rien
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              #*Bilan*
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              candidats <- C1 %>% bind_rows(C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22)
                                                                                                                                                                                                                                              candidats %>% group_by(ncRNA) %>% summarise(n=n())
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              22 couples/cluster s?lectionn?s --> Soit 14 ncRNA
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              14 ncRNA :
                                                                                                                                                                                                                                                *AT1G55525*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_1 <- candidats %>% filter(ncRNA == "AT1G55525") 
                                                                                                                                                                                                                                              candidat_1
                                                                                                                                                                                                                                              candidat_1 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall) %>% filter(ID %in% c("AT1G55525", "AT1G19220", "AT1G50460", "AT1G55580", "AT1G58340", "AT1G78240", "AT2G34650", "AT3G03660", "AT3G14370", "AT3G60630", "AT4G16780", "AT4G17460", "AT4G34410", "AT5G52050", "AT1G19220", "AT1G25220", "AT1G58340", "AT1G70940", "AT1G78240", "AT2G34650", "AT2G42430", "AT3G60630", "AT4G34410", "AT5G26930", "AT5G59030", "AT2G01830", "AT5G35750")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              GO root
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT1G55525", "AT1G19220", "AT1G50460", "AT1G55580", "AT1G58340", "AT1G78240", "AT2G34650", "AT3G03660", "AT3G14370", "AT3G60630", "AT4G16780", "AT4G17460", "AT4G34410", "AT5G52050", "AT1G19220", "AT1G25220", "AT1G58340", "AT1G70940", "AT1G78240", "AT2G34650", "AT2G42430", "AT3G60630", "AT4G34410", "AT5G26930", "AT5G59030", "AT2G01830", "AT5G35750")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              * AT1G19220	*ARF19	auxin response factor 19*	Encodes an auxin response factor that contains the conserved VP1-B3 DNA-binding domain at its N-terminus and the Aux/IAA-like domains III and IV present in most ARFs at its C-terminus. The protein interacts with IAA1 (yeast two hybrid) and other auxin response elements such as ER7 and ER9 (yeast one hybrid). ARF19 protein can complement many aspects of the arf7 mutant phenotype and , together with ARF7, is involved in the response to ethylene.  In the arf7 arf19 double mutant, several auxin-responsive genes (e.g. IAA5, LBD16, LBD29 and LBD33) are no longer upregulated by auxin.
                                                                                                                                                                                                                                              * AT1G25220	ASB1	anthranilate synthase beta subunit 1	*Catalyzes the first step of tryptophan biosynthesis*: Chorismate L-Glutamine = Anthranilate Pyruvate L-Glutamate. Functions as a heterocomplex with anthranilate synthase alpha subunit (ASA1 or ASA2).
                                                                                                                                                                                                                                              * AT1G50460	HKL1	hexokinase-like 1	Involved in glucose-ethylene crosstalk.
                                                                                                                                                                                                                                              AT1G55580	LAS	GRAS family transcription factor	Encodes a member of the GRAS family of putative transcriptional regulators. It is involved in the initiation of axillary meristems during both the vegetative and reproductive growth phases and functions upstream of REV and AXR1 in the regulation of shoot branching.
                                                                                                                                                                                                                                              * AT1G58340	ZF14	MATE efflux family protein	Encodes a plant MATE (multidrug and toxic compound extrusion) transporter that is localized to the Golgi complex and small organelles and is involved in determining the rate of organ initiation.  It is also involved in iron homeostasis when plants are under osmotic stress.
                                                                                                                                                                                                                                              * AT1G70940	*PIN3*	Auxin efflux carrier family protein	A regulator of auxin efflux and involved in differential growth. PIN3 is expressed in gravity-sensing tissues, with PIN3 protein accumulating predominantly at the lateral cell surface. PIN3 localizes to the plasma membrane and to vesicles. In roots, PIN3 is expressed without pronounced polarity in tiers two and three of the columella cells, at the basal side of vascular cells, and to the lateral side of pericycle cells of the elongation zone. PIN3 overexpression inhibits root cell growth. Protein phosphorylation plays a role in regulating PIN3 trafficking to  the plasma membrane.
                                                                                                                                                                                                                                              * AT1G78240	TSD2	S-adenosyl-L-methionine-dependent methyltransferases superfamily protein	Encodes TSD2 (TUMOROUS SHOOT DEVELOPMENT2), a putative methyltransferase with an essential role in cell adhesion and coordinated plant development.
                                                                                                                                                                                                                                              * AT2G01830	WOL	CHASE domain containing histidine kinase protein	*Histidine kinase: cytokinin-binding receptor that transduces cytokinin signals across the plasma membrane*
                                                                                                                                                                                                                                                * AT2G34650	*PID*	Protein kinase superfamily protein	Encodes a protein serine/threonine kinase that may act as a positive regulator of cellular auxin efflux, as a a binary switch for PIN polarity, and as a negative regulator of auxin signaling. Recessive mutants exhibit similar phenotypes as pin-formed mutants in flowers and inflorescence but distinct phenotypes in cotyledons and leaves. Expressed in the vascular tissue proximal to root and shoot meristems, shoot apex, and embryos. Expression is induced by auxin. Overexpression of the gene results in phenotypes in the root and shoot similar to those found in auxin-insensitive mutants. The protein physically interacts with TCH3 (TOUCH3) and PID-BINDING PROTEIN 1 (PBP1), a previously uncharacterized protein containing putative EF-hand calcium-binding motifs.  Acts together with ENP (ENHANCER OF PINOID) to instruct precursor cells to elaborate cotyledons in the transition stage embryo. Interacts with PDK1. PID autophosphorylation is required for the ability of PID to phosphorylate an exogenous substrate. PID activation loop is required for PDK1-dependent PID phosphorylation and requires the PIF domain. Negative regulator of root hair growth. PID kinase activity is critical for the inhibition of root hair growth and for maintaining  the proper subcellular localization of PID.
                                                                                                                                                                                                                                              * AT2G42430	*LBD16*	lateral organ boundaries-domain 16	LOB-domain protein gene LBD16. This gene contains one auxin-responsive element (AuxRE).
                                                                                                                                                                                                                                              * AT3G03660	*WOX11*	WUSCHEL related homeobox 11	Encodes a WUSCHEL-related homeobox gene family member with 65 amino acids in its homeodomain. Proteins in this family contain a sequence of eight residues (TLPLFPMH) downstream of the homeodomain called the WUS box.
                                                                                                                                                                                                                                              * AT3G14370	WAG2	Protein kinase superfamily protein	The WAG2 and its homolog, WAG1 each encodes protein-serine/threonine kinase that are nearly 70% identical to PsPK3 protein. All three together with CsPK3 belong to PsPK3-type kinases. At the N-terminus, all four possess a serine/threonine-rich domain. They are closely *related to Arabidopsis kinases PINOID*. wag1/wag2 double mutants exhibit a pronounced wavy root phenotype when grown vertically on agar plates (while wild-type plants develop wavy roots only on plates inclined to angles less than 90 degrees), indicating an overlapping role for WAG1 and WAG2 as suppressors of root waving. Simultaneous disruption of PID(AT2G34650) and its 3 closest homologs (PID2/AT2G26700, WAG1/AT1G53700, and WAG2/AT3G14370) abolishes the formation of cotyledons.
                                                                                                                                                                                                                                              * AT3G60630	HAM2	GRAS family transcription factor	Belongs to one of the LOM (LOST MERISTEMS) genes: AT2G45160 (LOM1), AT3G60630 (LOM2) and AT4G00150 (LOM3). LOM1 and LOM2 promote cell differentiation at the periphery of shoot meristems and help to maintain their polar organization.
                                                                                                                                                                                                                                              * AT4G16780	HB-2	homeobox protein 2	Encodes a homeodomain-leucine zipper protein that is rapidly and strongly induced by changes in the ratio of red to far-red light.  It is also involved in *cell expansion and cell proliferation and in the response to auxin*.
                                                                                                                                                                                                                                              * AT4G17460	HAT1	Homeobox-leucine zipper protein 4 (HB-4) / HD-ZIP protein	Encodes a class II HD-ZIP protein that regulates meristematic activity in different tissues, and that it is necessary for the correct formation of the gynoecium.
                                                                                                                                                                                                                                              * AT4G34410	RRTF1	redox responsive transcription factor 1	encodes a member of the ERF (ethylene response factor) subfamily B-3 of ERF/AP2 transcription factor family. The protein contains one AP2 domain. There are 18 members in this subfamily including ATERF-1, ATERF-2, AND ATERF-5.
                                                                                                                                                                                                                                              * AT5G26930	GATA23	GATA transcription factor 23	Encodes a member of the GATA factor family of zinc finger transcription factors. Controls lateral root founder cell specification.
                                                                                                                                                                                                                                              * AT5G35750	*HK2*	histidine kinase 2	Encodes histidine kinase AHK2.
                                                                                                                                                                                                                                              * AT5G52050	AT5G52050	MATE efflux family protein	""
                                                                                                                                                                                                                                              * AT5G59030	COPT1	copper transporter 1	encodes a putative copper transport protein that contains copper-binding motif and functionally complements in copper-transport defective yeast strains
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              AIA
                                                                                                                                                                                                                                              les inversement cor?l? (AT5G35750, AT2G01830) --> CK
                                                                                                                                                                                                                                              *TOP1*
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *AT2G31585*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_2 <- candidats %>% filter(ncRNA == "AT2G31585") 
                                                                                                                                                                                                                                              candidat_2 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT2G31585", "AT2G37630")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              Col_ler diff
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT2G31585", "AT2G37630")
                                                                                                                                                                                                                                              gene_list <- c("AT2G31585", "AT2G37630", "AT4G00220", "AT1G65620")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT2G32795*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_3 <- candidats %>% filter(ncRNA == "AT2G32795") 
                                                                                                                                                                                                                                              candidat_3 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT2G32795", "AT2G33790")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT2G32795", "AT2G33790")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G13495*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_4 <- candidats %>% filter(ncRNA == "AT4G13495") 
                                                                                                                                                                                                                                              candidat_4 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT4G13495", "AT1G23010")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              col sm 24
                                                                                                                                                                                                                                              col ler diff
                                                                                                                                                                                                                                              Pi resp
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G13495", "AT1G23010")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *AT4G34881*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_5 <- candidats %>% filter(ncRNA == "AT4G34881") 
                                                                                                                                                                                                                                              candidat_5 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("AT4G34881", "AT5G05560", "AT1G09560")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("AT4G34881", "AT5G05560", "AT1G09560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Ler_NEW_RNA_R_14094*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_6 <- candidats %>% filter(ncRNA == "Ler_NEW_RNA_R_14094") 
                                                                                                                                                                                                                                              candidat_6 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("Ler_NEW_RNA_R_14094", "AT2G33790")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Ler_NEW_RNA_R_14094", "AT2G33790")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *Ler_NEW_RNA_R_34309*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_7 <- candidats %>% filter(ncRNA == "Ler_NEW_RNA_R_34309") 
                                                                                                                                                                                                                                              candidat_7 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("Ler_NEW_RNA_R_34309", "AT4G05530")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("Ler_NEW_RNA_R_34309", "AT4G05530")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              De la merde !!!
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                *XLOC_001547*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_8 <- candidats %>% filter(ncRNA == "XLOC_001547") 
                                                                                                                                                                                                                                              candidat_8 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_001547", "AT1G12560")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_001547", "AT1G12560")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002093/AT1G08925*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_9 <- candidats %>% filter(ncRNA == "XLOC_002093") 
                                                                                                                                                                                                                                              candidat_9 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002093", "AT1G17110", "AT1G19220", "AT1G55580", "AT1G58340", "AT1G70560", "AT1G78240", "AT3G03660", "AT3G07390", "AT3G60630", "AT4G14550", "AT4G34410", "AT5G52050", "AT5G59030", "AT1G17110", "AT1G19220", "AT5G10720", "AT5G52050", "AT1G22530", "AT1G30690", "AT1G72150", "AT1G72160", "AT5G03150", "AT5G35750", "AT5G48000")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002093", "AT1G17110", "AT1G19220", "AT1G55580", "AT1G58340", "AT1G70560", "AT1G78240", "AT3G03660", "AT3G07390", "AT3G60630", "AT4G14550", "AT4G34410", "AT5G52050", "AT5G59030", "AT1G17110", "AT1G19220", "AT5G10720", "AT5G52050", "AT1G22530", "AT1G30690", "AT1G72150", "AT1G72160", "AT5G03150", "AT5G35750", "AT5G48000")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002755/AT1G08173*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_10 <- candidats %>% filter(ncRNA == "XLOC_002755") 
                                                                                                                                                                                                                                              candidat_10 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002755", "AT1G13260", "AT1G67710", "AT1G72150", "AT1G72160", "AT2G23430", "AT2G31090", "AT2G34680", "AT3G04630", "AT3G22400", "AT4G33880", "AT4G34580", "AT5G58010")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002755", "AT1G13260", "AT1G67710", "AT1G72150", "AT1G72160", "AT2G23430", "AT2G31090", "AT2G34680", "AT3G04630", "AT3G22400", "AT4G33880", "AT4G34580", "AT5G58010")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_002900/AT1G08687*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_11 <- candidats %>% filter(ncRNA == "XLOC_002900") 
                                                                                                                                                                                                                                              candidat_11 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_002900", "AT3G04630")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_002900", "AT3G04630")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_005603/AT3G03315*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_12 <- candidats %>% filter(ncRNA == "XLOC_005603") 
                                                                                                                                                                                                                                              candidat_12 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_005603", "AT1G14740", "AT1G16510", "AT1G25220", "AT2G42430", "AT5G12330", "AT1G14740", "AT3G03660", "AT3G07390", "AT3G11260", "AT4G01540", "AT5G52310", "AT1G77740", "AT5G40330")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_005603", "AT1G14740", "AT1G16510", "AT1G25220", "AT2G42430", "AT5G12330", "AT1G14740", "AT3G03660", "AT3G07390", "AT3G11260", "AT4G01540", "AT5G52310", "AT1G77740", "AT5G40330")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_007825*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_13 <- candidats %>% filter(ncRNA == "XLOC_007825") 
                                                                                                                                                                                                                                              candidat_13 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_007825", "AT1G05630")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_007825", "AT1G05630")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              *XLOC_011131*
                                                                                                                                                                                                                                                ```{r, include=FALSE}
                                                                                                                                                                                                                                              candidat_14 <- candidats %>% filter(ncRNA == "XLOC_011131") 
                                                                                                                                                                                                                                              candidat_14 %>% select(gene) %>% as.list()
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              all_info_V201810415 %>% select(ID, Araport11_ID, Araport11_overlap, type, width, Col_Ler,	t0t1,	t1t2,	t0t2,	interaction, longer_ORF_width,	nb_ORFs,	Col_siRNA_detected,	Ler_siRNA_detected,	Col_DicerCall,	Ler_DicerCall,	siRNA_Col_Ler,	GO) %>% filter(ID %in% c("XLOC_011131", "AT3G18780", "AT5G20490")) %>% arrange(desc(type))
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, echo=FALSE}
                                                                                                                                                                                                                                              gene_list <- c("XLOC_011131", "AT3G18780", "AT5G20490")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              normcounts_current <- normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              stat_normcounts_current <- stat_normcounts_tidy %>%
                                                                                                                                                                                                                                                filter(ID %in% gene_list) 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ```
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ESSAYER tema regulation en CIS !
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ```{r, include=FALSE}
                                                                                                                                                                                                                                              stat_normcounts_current_Ler <- stat_normcounts_current %>% filter(ecotype == "ler")
                                                                                                                                                                                                                                              normcounts_current_Ler <- normcounts_current %>% filter(ecotype == "ler")
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              ggplot(data=stat_normcounts_current_Ler, aes(x=time, color=ecotype)) +
                                                                                                                                                                                                                                                geom_point(data=normcounts_current_Ler, aes(y=counts)) +
                                                                                                                                                                                                                                                geom_line(aes(y=mean_counts, group=ecotype)) +
                                                                                                                                                                                                                                                geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
                                                                                                                                                                                                                                                facet_wrap(~ID, scale="free")
                                                                                                                                                                                                                                              ``` 
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              