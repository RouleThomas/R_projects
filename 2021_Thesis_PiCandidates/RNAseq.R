# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/Pearson/", showWarnings = FALSE, recursive = TRUE)


# Genevestigator perturbation data -----


#AT5G09710 GV ------
correlation_col_100 %>% filter(ncRNA == "AT5G09710") %>% select(gene)

GV <- read_excel("data/Pearson_GV_Araport/AT5G09710_GV.xlsx", col_names=FALSE)


gene_info <-
    GV[1:6, 6:ncol(GV)] %>% # =[remove, take] : remove the 6 first columns, take all columns from the 6
    t %>% # transpose = inverse x and y axis
    as.data.frame %>% # make a table
    as_tibble %>% # make a tibble
    setNames(.[1,] %>% unlist %>% as.character) %>% # use first row as header
    filter(row_number() > 1) %>% # remove values from line 1
    unique # remove duplicate

exp_info <-
    GV[7:nrow(GV), 1:5] %>%
    as.data.frame %>%
    as_tibble %>%
    setNames(.[1,] %>% unlist %>% as.character) %>%
    filter(row_number() > 1) %>%
    group_by(Experiment) %>% # group by Experiment : mutate will focus on Epxperiment column
    mutate(comp=str_c(Experiment, "-", row_number())) %>% # str_c(x,y) = will paste x and y, here for each Experiment, give the number of condition
    ungroup %>%
    unique

GV_formated_data <-
    GV %>%
    dplyr::select(-c(1, 2, 4:6)) %>% # do not select those colums
    filter(!row_number() %in% c(1:3, 5:6)) %>% # do not select those rows
    setNames(str_c(.[1,] %>% unlist %>% as.character %>% replace_na(""), #replace_na by empty
                   .[2,] %>% unlist %>% as.character, sep=":")) %>% # rename header as: row1:row2
    filter(!row_number() %in% 1:2) %>%
    dplyr::rename(Experiment=`:Experiment`) %>%
    group_by(Experiment) %>%
    mutate(comp=str_c(Experiment, "-", row_number())) %>%
    ungroup %>%
    dplyr::select(-Experiment) %>%
    pivot_longer(-comp, names_to="names", values_to="value") %>% # pivot_longer(-all column except comp need to be reshape, names_to=take the header column as value  )
    separate(names, c('Measure', "type"), sep=":")   # in a new column, values_to = title of the corresponding value)


cluster <- data.frame(Gene = c("AT5G09710", "AT3G60630", "AT3G07390", "AT1G70560", "AT1G78240", "AT4G34410", "AT5G52050", "AT3G03660", "AT1G17110", "AT3G14370", "AT1G19220", "AT1G58340", "AT2G36400", "AT5G10720", "AT4G17460", "AT2G34650", "AT5G59030", "AT1G70940", "AT1G55580", "AT4G16780", "AT1G50460", "AT4G09510"),
                               type = c("AT5G09710", "AT3G60630", "AT3G07390", "AT1G70560", "AT1G78240", "AT4G34410", "AT5G52050", "AT3G03660", "AT1G17110", "AT3G14370", "AT1G19220", "AT1G58340", "AT2G36400", "AT5G10720", "AT4G17460", "AT2G34650", "AT5G59030", "AT1G70940", "AT1G55580", "AT4G16780", "AT1G50460", "AT4G09510"))



GV_formated_data <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  left_join(cluster) %>% 
  mutate(type = fct_explicit_na(type, na_level = "None"), type=as.character(type)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
         ) %>%
  ungroup




Correlation_plot <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


res1 <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

pdf("out/Pearson/corrplot_AT5G09710_GV.pdf", width=7, height=7)

corrplot(Correlation_plot, p.mat = res1$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


dev.off()


#AT1G55525 GV ------
correlation_col_100 %>% filter(ncRNA == "AT1G55525") %>% select(gene) 


AT1G55525 <-correlation_col_100 %>% filter(ncRNA == "AT1G55525") %>% select(gene) %>% rbind("AT1G55525")



GV <- read_excel("data/Pearson_GV_Araport/AT1G55525_GV.xlsx", col_names=FALSE)


gene_info <-
  GV[1:6, 6:ncol(GV)] %>% # =[remove, take] : remove the 6 first columns, take all columns from the 6
  t %>% # transpose = inverse x and y axis
  as.data.frame %>% # make a table
  as_tibble %>% # make a tibble
  setNames(.[1,] %>% unlist %>% as.character) %>% # use first row as header
  filter(row_number() > 1) %>% # remove values from line 1
  unique # remove duplicate

exp_info <-
  GV[7:nrow(GV), 1:5] %>%
  as.data.frame %>%
  as_tibble %>%
  setNames(.[1,] %>% unlist %>% as.character) %>%
  filter(row_number() > 1) %>%
  group_by(Experiment) %>% # group by Experiment : mutate will focus on Epxperiment column
  mutate(comp=str_c(Experiment, "-", row_number())) %>% # str_c(x,y) = will paste x and y, here for each Experiment, give the number of condition
  ungroup %>%
  unique

GV_formated_data <-
  GV %>%
  dplyr::select(-c(1, 2, 4:6)) %>% # do not select those colums
  filter(!row_number() %in% c(1:3, 5:6)) %>% # do not select those rows
  setNames(str_c(.[1,] %>% unlist %>% as.character %>% replace_na(""), #replace_na by empty
                 .[2,] %>% unlist %>% as.character, sep=":")) %>% # rename header as: row1:row2
  filter(!row_number() %in% 1:2) %>%
  dplyr::rename(Experiment=`:Experiment`) %>%
  group_by(Experiment) %>%
  mutate(comp=str_c(Experiment, "-", row_number())) %>%
  ungroup %>%
  dplyr::select(-Experiment) %>%
  pivot_longer(-comp, names_to="names", values_to="value") %>% # pivot_longer(-all column except comp need to be reshape, names_to=take the header column as value  )
  separate(names, c('Measure', "type"), sep=":")   # in a new column, values_to = title of the corresponding value)


cluster <- data.frame(Gene = c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220"),
                      type = c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220"))




GV_formated_data <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  left_join(cluster) %>% 
  mutate(type = fct_explicit_na(type, na_level = "None"), type=as.character(type)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup




Correlation_plot <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


res1 <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

pdf("out/Pearson/corrplot_AT1G55525_GV.pdf", width=7, height=7)

corrplot(Correlation_plot, p.mat = res1$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


dev.off()


#XLOC_002093/AT1G08925_ no data for GV ------

#XLOC_002755/AT1G08173 no data for GV ------

#XLOC_002900/AT1G08687 no data for GV ------


#AT5G38005 GV ------

correlation_col_100 %>% filter(ncRNA == "AT5G38005") %>% select(gene)

GV <- read_excel("data/Pearson_GV_Araport/AT5G38005_GV.xlsx", col_names=FALSE)


gene_info <-
  GV[1:6, 6:ncol(GV)] %>% # =[remove, take] : remove the 6 first columns, take all columns from the 6
  t %>% # transpose = inverse x and y axis
  as.data.frame %>% # make a table
  as_tibble %>% # make a tibble
  setNames(.[1,] %>% unlist %>% as.character) %>% # use first row as header
  filter(row_number() > 1) %>% # remove values from line 1
  unique # remove duplicate

exp_info <-
  GV[7:nrow(GV), 1:5] %>%
  as.data.frame %>%
  as_tibble %>%
  setNames(.[1,] %>% unlist %>% as.character) %>%
  filter(row_number() > 1) %>%
  group_by(Experiment) %>% # group by Experiment : mutate will focus on Epxperiment column
  mutate(comp=str_c(Experiment, "-", row_number())) %>% # str_c(x,y) = will paste x and y, here for each Experiment, give the number of condition
  ungroup %>%
  unique

GV_formated_data <-
  GV %>%
  dplyr::select(-c(1, 2, 4:6)) %>% # do not select those colums
  filter(!row_number() %in% c(1:3, 5:6)) %>% # do not select those rows
  setNames(str_c(.[1,] %>% unlist %>% as.character %>% replace_na(""), #replace_na by empty
                 .[2,] %>% unlist %>% as.character, sep=":")) %>% # rename header as: row1:row2
  filter(!row_number() %in% 1:2) %>%
  dplyr::rename(Experiment=`:Experiment`) %>%
  group_by(Experiment) %>%
  mutate(comp=str_c(Experiment, "-", row_number())) %>%
  ungroup %>%
  dplyr::select(-Experiment) %>%
  pivot_longer(-comp, names_to="names", values_to="value") %>% # pivot_longer(-all column except comp need to be reshape, names_to=take the header column as value  )
  separate(names, c('Measure', "type"), sep=":")   # in a new column, values_to = title of the corresponding value)


cluster <- data.frame(Gene = c("AT5G38005","AT5G57620","AT1G56010","AT4G14716","AT3G54870","AT1G72160"),
                      type = c("AT5G38005","AT5G57620","AT1G56010","AT4G14716","AT3G54870","AT1G72160"))




GV_formated_data <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  left_join(cluster) %>% 
  mutate(type = fct_explicit_na(type, na_level = "None"), type=as.character(type)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup




Correlation_plot <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


res1 <- GV_formated_data %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

pdf("out/Pearson/corrplot_AT1G38005_GV.pdf", width=7, height=7)

corrplot(Correlation_plot, p.mat = res1$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


dev.off()


#XLOC_005603/AT3G03315 no data for GV ------

# Araport11 perturbation data -----
Araport11_all_counts <- 
  read_tsv("data/Pearson_GV_Araport/Araport11_all_gene_counts.tsv.gz") %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "normalised_count")
#AT5G09710-----
cluster <- c("AT5G09710", "AT3G60630", "AT3G07390", "AT1G70560", "AT1G78240", "AT4G34410", "AT5G52050", "AT3G03660", "AT1G17110", "AT3G14370", "AT1G19220", "AT1G58340", "AT2G36400", "AT5G10720", "AT4G17460", "AT2G34650", "AT5G59030", "AT1G70940", "AT1G55580", "AT4G16780", "AT1G50460", "AT4G09510")


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

pdf("out/Pearson/AT5G09710_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()

median_AT5G09710 <- 
as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT5G09710")



#AT1G55525-----
cluster <- c("AT1G55525", "AT5G52050", "AT4G34410", "AT2G34650", "AT1G58340", "AT1G55580", "AT3G03660", "AT1G78240", "AT3G60630", "AT4G17460", "AT4G16780", "AT3G14370", "AT1G50460", "AT1G19220")


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

pdf("out/Pearson/AT1G55525_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


median_AT1G55525 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT1G55525")



#XLOC_002093/AT1G08925-----
cluster <- c("AT1G08925", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110", "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550")
XLOC_002093 <- c("XLOC_002093", "AT1G78240", "AT3G07390", "AT3G03660", "AT1G70560", "AT4G34410", "AT1G58340", "AT3G60630", "AT1G17110",
                 "AT1G19220", "AT5G52050", "AT5G59030", "AT1G55580", "AT4G14550") %>% as.tibble %>% rename(gene="value")


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

pdf("out/Pearson/AT1G08925_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()

median_AT1G08925 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT1G08925")



#XLOC_002755/AT1G08173-----
cluster <- c("AT1G08173", "AT2G23430", "AT1G72160", "AT5G58010", "AT1G67710", "AT3G04630", "AT1G13260", "AT2G34680", "AT1G72150", "AT3G22400", "AT4G33880", "AT4G34580", "AT2G31090")


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

pdf("out/Pearson/AT1G08173_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


median_AT1G08173 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT1G08173")



#XLOC_002900/AT1G08687-----
cluster <- c("AT1G08687", "AT4G17500", "AT3G04630", "AT5G61600", "AT1G67710", "AT1G56010", "AT2G26670")
XLOC_002900 <- c("XLOC_002900", "AT4G17500", "AT3G04630", "AT5G61600", "AT1G67710", "AT1G56010", "AT2G26670") %>% as.tibble %>% rename(gene=value)


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

pdf("out/Pearson/AT1G08687_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()



median_AT1G08687 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT1G08687")



#AT5G38005-----
cluster <- c("AT5G38005","AT5G57620","AT1G56010","AT4G14716","AT3G54870","AT1G72160")


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

pdf("out/Pearson/AT5G38005_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()

median_AT5G38005 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT5G38005")

#XLOC_005603/AT3G03315-----
cluster <- c("AT3G03315", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430")
XLOC_005603 <- c("XLOC_005603", "AT1G25220", "AT1G14740", "AT5G12330", "AT1G16510", "AT2G42430") %>% as.tibble %>% rename(gene=value)


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

pdf("out/Pearson/AT3G03315_Araport.pdf", width=9, height=9)
corrplot(cluster_corr, p.mat = cluster_sig$p, method="color",
         col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex = .9, tl.col="black", 
         diag = FALSE)
dev.off()


median_AT3G03315 <- 
  as.tibble(cluster_corr) %>% 
  filter(row_number()==1) %>% 
  select(-(1)) %>% 
  t %>% 
  as.tibble() %>% 
  summarise(median=median(V1), sd=sd(V1), mean=mean(V1), n=n()) %>%
  mutate(se=sd/n) %>%
  add_column("ncRNA"="AT3G03315")


# median correlation all candidats -----

all_median <- median_AT1G08173 %>%
  bind_rows(median_AT1G08687,median_AT1G08925, median_AT1G55525, median_AT3G03315, median_AT5G09710, median_AT5G38005) %>% 
  arrange(median)




my_graph <- 
  all_median  %>%
  ggplot(data=., mapping=aes(x=ncRNA, y=median)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin= (median-se), ymax= (median+se), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "median Pearson correlation")+
  geom_hline(yintercept = 0, linetype=1)
my_graph

my_graph <- 
  all_median  %>%
  ggplot(data=., mapping=aes(x=ncRNA, y=mean)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin= (mean-se), ymax= (mean+se), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_text(aes(label=n), nudge_y = +0.06, size=4)+
  xlab(label = "")+
  ylab(label = "mean Pearson correlation")+
  geom_hline(yintercept = 0, linetype=1)
my_graph



?geom_errorbar
?summarise

