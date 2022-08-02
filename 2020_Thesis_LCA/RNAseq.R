
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Data import ----

GV <- read_excel("data/RNAseq/GV_lca_region.xlsx", col_names=FALSE)


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


GV_formated_data_lca <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup


GV_formated_data_lca_order <- GV_formated_data_lca %>%
  filter(Gene == "AT5G38005") %>% 
  arrange(value) %>%
  dplyr::select(comp) %>%
  mutate(row_idx=as.factor(row_number()))

GV_formated_data_lca_order <-
  GV_formated_data_lca %>%
  left_join(GV_formated_data_lca_order)


# Graphs geom_bind2d

ggplot(GV_formated_data_lca_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=value)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("raw")

ggplot(GV_formated_data_lca_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=scaled)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("scaled")

ggplot(GV_formated_data_lca_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=normalised)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("normalised")



# Correlation plots
#
# Check different option there:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

library(corrplot)

# Final graphs -----

Correlation_plot <- GV_formated_data_lca %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


res1 <- cor.mtest(Correlation_plot, conf.level = .95)

## specialized the insignificant value according to the significant level

pdf("out/RNAseq/corrplot_GV.pdf", width=7, height=7)

corrplot(Correlation_plot, p.mat = res1$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


dev.off()



# RNAseq flg22 ------
GSE63603_FPKM_for_all_samples <- read_delim("C:/Users/roule/OneDrive/Bureau/Covid19 Work/Paper_lca/DAta_ngs/GSE63603_FPKM_for_all_samples.txt", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)



normcounts_tidy <- GSE63603_FPKM_for_all_samples %>% 
  gather(`WT_CK_R1_ FPKM`:`OX9_flg22_R2_FPKM`, key = "condition", value = "counts", na.rm = TRUE) %>%
  separate(condition, c("genotype", "condition", "replicate"), sep="_") %>%
  dplyr::select(-X14,-X15,-X16,-X17)

stat_normcounts_tidy <- normcounts_tidy %>%
  group_by(AGI, genotype, condition) %>%
  summarise(mean_counts=mean(counts), 
            median_counts=median(counts),
            ecart_type=sd(counts), #ecart-type
            n=n(), #nombre d'échantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart


# List gene  -------------

gene_list <- c("AT5G37900","AT5G37910","AT5G05565","AT5G37920","AT5G37930",
               "AT5G37940","AT5G05585","AT5G37950","AT5G37960","AT5G37970",
               "AT5G37980","AT5G37990","AT5G05605","AT5G38000","AT5G38005",
               "AT5G38010","AT5G38020","AT5G38030","AT5G38035","AT5G38037",
               "AT5G38040","AT5G38050","AT5G38060","AT5G38070","AT5G38080",
               "AT5G38090","AT5G38096","AT5G38100")

gene_list <- c("AT5G37960","AT5G37970",
               "AT5G37980","AT5G37990","AT5G05605","AT5G38000","AT5G38005",
               "AT5G38010","AT5G38020","AT5G38030","AT5G38035","AT5G38037",
               "AT5G38040","AT5G38050","AT5G38060")



normcounts_current <- normcounts_tidy %>%
  filter(AGI %in% gene_list) #Selection des genes en etudes dans normcounts


stat_normcounts_current <-
  stat_normcounts_tidy %>%
  filter(AGI %in% gene_list) # Stat que sur les genes en études


#Graphs
ggplot(data=stat_normcounts_current, aes(x=condition, y=mean_counts, fill=genotype)) +
  geom_col(position=position_dodge(.9)) +
  geom_errorbar(position=position_dodge(.9), aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=0.5) +
  facet_wrap(~AGI, scale="free") +
  theme_classic()


# clf -------

GSE61545_genes_fpkm_tracking <- read_excel("data/RNAseq/GSE61545_genes.fpkm_tracking.xlsx") %>%
  gather(clf_inflor_r1:WT_siliques_r3, key = "condition", value = "counts", na.rm = TRUE) %>%
  separate(condition, c("genotype", "organ", "replicate"), sep="_")



stat_normcounts_tidy <- GSE61545_genes_fpkm_tracking %>%
  group_by(ID, genotype, organ) %>%
  summarise(mean_counts=mean(counts), 
            median_counts=median(counts),
            ecart_type=sd(counts), #ecart-type
            n=n(), #nombre d'échantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart


# List gene

gene_list <- c("AT5G37960", "AT5G37970", "AT5G37980", "AT5G37990", "AT5G38000", "AT5G38000", "AT5G38005", "AT5G38010", "AT5G38020", "AT5G38030", "AT5G38035"
               , "AT5G38037", "AT5G38040")

gene_list <- c("AT5G37900","AT5G37910","AT5G05565","AT5G37920","AT5G37930",
               "AT5G37940","AT5G05585","AT5G37950","AT5G37960","AT5G37970",
               "AT5G37980","AT5G37990","AT5G05605","AT5G38000","AT5G38005",
               "AT5G38010","AT5G38020","AT5G38030","AT5G38035","AT5G38037",
               "AT5G38040","AT5G38050","AT5G38060","AT5G38070","AT5G38080",
               "AT5G38090","AT5G38096","AT5G38100")

normcounts_current <- GSE61545_genes_fpkm_tracking %>%
  filter(ID %in% gene_list) #Selection des genes en etudes dans normcounts


stat_normcounts_current <-
  stat_normcounts_tidy %>%
  filter(ID %in% gene_list) # Stat que sur les genes en études


#Graphs
stat_normcounts_current %>% filter(genotype =="WT", organ %in% c("shoots","roots")) %>%
ggplot(data=, aes(x=organ, y=mean_counts, fill=organ)) +
  geom_col(position=position_dodge(.9)) +
  geom_errorbar(position=position_dodge(.9), aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=0.5) +
  facet_wrap(~ID, scale="free") +
  theme_classic()



