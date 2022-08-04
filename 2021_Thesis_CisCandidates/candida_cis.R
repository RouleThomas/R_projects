
# Find IAA14 ----------



all_info_V201810415 <- read_excel("data/CIS/all_info_V201810415.xlsx")
GO_root <- read_delim("data/CIS/GO root.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)%>%
  rename(ID="Gene.primaryIdentifier")



#1. isolates the lncRNA precursor of smRNA 24nt non NAT detected in col = 462 unique lncRNA
all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200) %>%
  select(ID)%>%
  unique() # 1671 lincRNA

All_candidates <-all_info_V201810415 %>%
  filter(type=="noncoding",width>200) %>%
  select(ID)%>%
  unique() # 8628 lincRNA
write.table(All_candidates, file="X8628_candidates.csv")



all_info_V201810415 %>%
  filter(type=="coding") %>%
  select(ID)%>%
  unique() # 29,657 lincRNA

neigb_count_up <- all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200) %>%
  select(upstream_ID)%>%
  unique()%>%
  rename(ID=upstream_ID)
neigb_count_down <- all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200) %>%
  select(downstream_ID)%>%
  unique()%>%
  rename(ID=downstream_ID)
neigb_count_up %>%
  bind_rows(neigb_count_down) %>%
  inner_join(GO_root)# 46 genes GO root close to lincRNA

all_info_V201810415_24 <- all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200, Col_siRNA_detected=="TRUE", Col_DicerCall==24, Col_precursor=="siRNA")

X456_candidates <- all_info_V201810415_24 %>% select(ID) %>% unique() # 456 lincRNA col precursor sirNA


write.table(X456_candidates, file="X456_candidates.csv")

all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200, Col_siRNA_detected=="TRUE")%>% select(ID) %>% unique()


#2. isolates the one that are neigb to a gene GO root 

##upstream: 2 lincRNA
all_info_V201810415_24_up <- all_info_V201810415_24 %>%
  select(ID, upstream_ID) %>%
  rename(ID_lincRNA=ID, ID=upstream_ID) %>%
  inner_join(GO_root)


##downstream: 3 lincRNA
all_info_V201810415_24_down <- all_info_V201810415_24 %>%
  select(ID, downstream_ID) %>%
  rename(ID_lincRNA=ID,ID=downstream_ID)%>%
  inner_join(GO_root)


##As so few I stop there and look at the candidates

# Find MARS ----------


#1. isolates the lncRNA non NAT more enriched in Col


all_info_V201810415_lincRNA_col <- all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200, Col_Ler>0)

all_info_V201810415_lincRNA_col %>% select(ID,Araport11_ID)%>%unique
#1599 lincRNA Col enriched
all_info_V201810415 %>%
  filter(type=="coding",Col_Ler>0)%>%
  select(ID)%>%
  unique()
#29,063 coding Col enriched


#2. isolates the neigb also enriched in Col

##upstream: 
all_info_V201810415_lincRNA_up <- all_info_V201810415_lincRNA_col %>%
  select(ID, upstream_ID) %>%
  rename(ID_lincRNA=ID, ID=upstream_ID) %>%
  inner_join(all_info_V201810415) %>%
  filter(NAT=="FALSE",type=="coding",Col_Ler>0, Col_Ler != "NA")



##downstream:

all_info_V201810415_lincRNA_down <- all_info_V201810415_lincRNA_col %>%
  select(ID, downstream_ID) %>%
  rename(ID_lincRNA=ID, ID=downstream_ID) %>%
  inner_join(all_info_V201810415) %>%
  filter(NAT=="FALSE",type=="coding",Col_Ler>0, Col_Ler != "NA")

all_info_V201810415_lincRNA_up %>% 
  bind_rows(all_info_V201810415_lincRNA_down)%>%
  select(ID) %>%
  unique()
# 128 coding genes next to lincRNA both col enriched



# Does neigb have GO root (hope no); only two and they are shit

all_info_V201810415_coding_neigb <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID) %>%
  unique() %>%
  inner_join(GO_root)
#

#3. Lot of candidates so add Pearson correlation analysis between lincRNA and its neigb


normcounts <- read_delim("data/CIS/normcounts.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ID=X1)

lincRNA_coding_pair <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID_lincRNA,ID) %>%
  unique()

#isolates the lincRNA identified above (135 in total)

nc_gene_count <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID_lincRNA)%>%
  rename(ID=ID_lincRNA) %>%
  unique()%>%
  left_join(normcounts)%>%
  column_to_rownames("ID")

#isolates the coding
coding_gene_count <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID)%>%
  unique()%>%
  left_join(normcounts)%>%
  column_to_rownames("ID")



pearson_correlation_df <- # correlation des count ncRNA et count coding gene
  cor(t(nc_gene_count), t(coding_gene_count), method="pearson") %>% # transposition = inversion des labels colonne/ligne
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% # colonne labels portant les ncRNA devient colonne à part label "ncRNA"
  gather(gene, pearson_corr, -ncRNA) %>% # gather toutes les colonnes sauf la colonne "ncRNA" qui reste individuelle et valeur des couples nommer "pearson_corr"
  filter(ncRNA %in% row.names(nc_gene_count), # colonne ncRNA = que les lncRNA
         gene %in% row.names(coding_gene_count))  # colonne gene = que les coding gene


#4. collect the lincRNA_pair pearson score

all_score <- pearson_correlation_df %>% 
  rename(ID_lincRNA=ncRNA,
         ID=gene) %>%
  select(ID_lincRNA, ID, pearson_corr) %>%
  as_tibble() %>%
  inner_join(lincRNA_coding_pair)%>%
  filter(pearson_corr>0.5) %>%
  arrange(desc(pearson_corr))


pearson_correlation_df %>% 
  rename(ID_lincRNA=ncRNA,
         ID=gene) %>%
  select(ID_lincRNA, ID, pearson_corr) %>%
  as_tibble() %>%
  inner_join(lincRNA_coding_pair) %>%
  filter(pearson_corr>0.5) %>%
  arrange(desc(pearson_corr))

write.csv(all_score,file="all_score.csv")


# No CLE14 before !!!!!! Try not Col enriched -----------



#1. isolates the lncRNA non NAT more enriched in Col


all_info_V201810415_lincRNA_col <- all_info_V201810415 %>%
  filter(NAT==FALSE, type=="noncoding",width>200, Col_Ler>0)


#2. isolates the neigb also enriched in Col

##upstream: 
all_info_V201810415_lincRNA_up <- all_info_V201810415_lincRNA_col %>%
  select(ID, upstream_ID) %>%
  rename(ID_lincRNA=ID, ID=upstream_ID) %>%
  inner_join(all_info_V201810415) %>%
  filter(NAT=="FALSE",type=="coding")

##downstream:

all_info_V201810415_lincRNA_down <- all_info_V201810415_lincRNA_col %>%
  select(ID, downstream_ID) %>%
  rename(ID_lincRNA=ID, ID=downstream_ID) %>%
  inner_join(all_info_V201810415) %>%
  filter(NAT=="FALSE",type=="coding")


# Does neigb have GO root (hope no); only two and they are shit

all_info_V201810415_coding_neigb <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID) %>%
  unique() %>%
  inner_join(GO_root)
#

#3. Lot of candidates so add Pearson correlation analysis between lincRNA and its neigb


normcounts <- read_delim("data/CIS/normcounts.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ID=X1)

lincRNA_coding_pair <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID_lincRNA,ID) %>%
  unique()

#isolates the lincRNA identified above (135 in total)

nc_gene_count <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID_lincRNA)%>%
  rename(ID=ID_lincRNA) %>%
  unique()%>%
  left_join(normcounts)%>%
  column_to_rownames("ID")

#isolates the coding
coding_gene_count <- all_info_V201810415_lincRNA_up %>%
  bind_rows(all_info_V201810415_lincRNA_down) %>%
  select(ID)%>%
  unique()%>%
  left_join(normcounts)%>%
  column_to_rownames("ID")



pearson_correlation_df <- # correlation des count ncRNA et count coding gene
  cor(t(nc_gene_count), t(coding_gene_count), method="pearson") %>% # transposition = inversion des labels colonne/ligne
  as.data.frame  %>% # Matrice de correlation
  rownames_to_column("ncRNA") %>% # colonne labels portant les ncRNA devient colonne à part label "ncRNA"
  gather(gene, pearson_corr, -ncRNA) %>% # gather toutes les colonnes sauf la colonne "ncRNA" qui reste individuelle et valeur des couples nommer "pearson_corr"
  filter(ncRNA %in% row.names(nc_gene_count), # colonne ncRNA = que les lncRNA
         gene %in% row.names(coding_gene_count))  # colonne gene = que les coding gene


#4. collect the lincRNA_pair pearson score

all_score <- pearson_correlation_df %>% 
  rename(ID_lincRNA=ncRNA,
         ID=gene) %>%
  select(ID_lincRNA, ID, pearson_corr) %>%
  as_tibble() %>%
  inner_join(lincRNA_coding_pair)


pearson_correlation_df %>% 
  rename(ID_lincRNA=ncRNA,
         ID=gene) %>%
  select(ID_lincRNA, ID, pearson_corr) %>%
  as_tibble() %>%
  inner_join(lincRNA_coding_pair) %>%
  filter(pearson_corr>0.5)


write.csv(pearson_correlation_df %>% 
            rename(ID_lincRNA=ncRNA,
                   ID=gene) %>%
            select(ID_lincRNA, ID, pearson_corr) %>%
            as_tibble() %>%
            inner_join(lincRNA_coding_pair) %>%
            filter(pearson_corr>0.5), file = "mars_candidats_isolation.csv")

### Doesnt work CLE14 is badely corelated with the lincRNA
