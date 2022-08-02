# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/RNAseq/", showWarnings = FALSE, recursive = TRUE)




#data import
normcounts <- read_excel("data/RNAseq/normcounts.xlsx")

normcounts_tidy <- normcounts %>% 
  gather("fresh_1":"SL48_3", key = "condition", value = "counts", na.rm = TRUE) %>%
  separate(condition, c("condition", "replicate"), sep="_")


stat_normcounts_tidy <- normcounts_tidy %>%
  group_by(ID, condition) %>%
  summarise(mean_counts=mean(counts), 
            median_counts=median(counts),
            ecart_type=sd(counts), #ecart-type
            n=n(), #nombre d'échantillons
            erreur_std=ecart_type/sqrt(n)) #erreur standart


# List gene  -------------

gene_list <- c("AT5G42600", "AT5G42580", "AT5G42590", "AT5G00580")
gene_list <- c("AT5G42600", "AT5G42580", "AT5G42590", "AT5G00580","AT5G06335","AT5G06325")

normcounts_current <- normcounts_tidy %>%
  filter(ID %in% gene_list) 

stat_normcounts_current <-
  stat_normcounts_tidy %>%
  filter(ID %in% gene_list) # Stat que sur les genes en études

stat_normcounts_current$condition <- factor(stat_normcounts_current$condition, levels=c("fresh", "dry", "S1", "S12", "S48", "SL1", "SL6", 
                                                                                        "SL12", "SL24", "SL48")) # Choisir ordrer du facte_wrap

normcounts_current$condition <- factor(normcounts_current$condition, levels=c("fresh", "dry", "S1", "S12", "S48", "SL1", "SL6", 
                                                                              "SL12", "SL24", "SL48"))



ggplot(data=stat_normcounts_current, aes(x=condition, color=ID)) +
  geom_point(data=normcounts_current, aes(y=counts)) +
  geom_line(aes(y=mean_counts)) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.1) +
  facet_wrap(~ID, scale="free", nrow=1) +
  theme_bw()

stat_normcounts_current 


stat_normcounts_current_name <- stat_normcounts_current %>%
  mutate(gene=ifelse(ID=="AT5G00580", "MHAL", NA)) %>%
  mutate(gene=ifelse(ID=="AT5G42600", "MRN1", gene)) %>%
  mutate(gene=ifelse(ID=="AT5G42580", "CYP705A12", gene)) %>%
  mutate(gene=ifelse(ID=="AT5G42590", "CYP71A16", gene))

my_graph <- 
stat_normcounts_current_name %>% left_join(selected_gene_qPCR) %>% mutate(gene=as.factor(gene)) %>%
ggplot(data=., aes(x=condition)) +
  geom_line(aes(y=mean_counts, group=gene, color=gene), size = 1) +
  geom_errorbar(aes(ymin=mean_counts - erreur_std, ymax=mean_counts + erreur_std), width=.2) +
  scale_color_manual(breaks=selected_gene_RNAseq$gene,
                    values=selected_gene_RNAseq$gene_color,
                    labels=selected_gene_RNAseq$gene_name,) +
  xlab("") +
  ylab("") +
  theme(axis.text.x  = element_text(angle=45, vjust=1, size=10, hjust=1))
my_graph

#Save
ggsave(filename="out/RNAseq/marneral_cluster_germination.pdf", plot=my_graph, width = 5, height = 3) 


# Genevestigator perturbation data -----
#data import

GV <- read_excel("data/RNAseq/GV_marneral_region.xlsx", col_names=FALSE)

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


marneral_cluster <- data.frame(Gene = c("AT5G42580","AT5G42590","AT5G42600"),
                               type = c("AT5G42580","AT5G42590","AT5G42600"))

GV_formated_data_marneral <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  left_join(marneral_cluster) %>% 
  mutate(type = fct_explicit_na(type, na_level = "None"), type=as.character(type)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
         ) %>%
  ungroup



GV_formated_data_MRN1_order <- GV_formated_data_marneral %>%
  filter(Gene == "AT5G42600") %>% 
  arrange(value) %>%
  dplyr::select(comp) %>%
  mutate(row_idx=as.factor(row_number()))

GV_formated_data_marneral_order <-
  GV_formated_data_marneral %>%
  left_join(GV_formated_data_MRN1_order)
  

# Graphs geom_bind2d
  
ggplot(GV_formated_data_marneral, aes(comp, Gene)) + 
  geom_bin2d(aes(fill=value)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("raw")
  
ggplot(GV_formated_data_marneral, aes(comp, Gene)) + 
  geom_bin2d(aes(fill=scaled)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("scaled")
  
ggplot(GV_formated_data_marneral, aes(comp, Gene)) + 
  geom_bin2d(aes(fill=normalised)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("normalised")

# Graphs geom_line


GV_formated_data_marneral %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=comp, y=value, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green")) +
  ggtitle("raw")

GV_formated_data_marneral %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=comp, y=scaled, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green")) +
  ggtitle("scaled")

GV_formated_data_marneral %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=comp, y=normalised, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green"))
  ggtitle("normalised")


## Ordered according to MRN1 value
  
# Graphs geom_bin2d


ggplot(GV_formated_data_marneral_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=value)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("raw")

  
ggplot(GV_formated_data_marneral_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=scaled)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("scaled")
  
ggplot(GV_formated_data_marneral_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=normalised)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("normalised")

# Graphs geom_line


GV_formated_data_marneral_order %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=row_idx, y=value, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green")) +
  ggtitle("raw")

GV_formated_data_marneral_order %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=row_idx, y=scaled, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green")) +
  ggtitle("scaled")

GV_formated_data_marneral_order %>%
  filter(Gene %in% c("AT5G42550","AT5G42560","AT5G42570","AT5G42580","AT5G42590","AT5G42600","AT5G42610","AT5G42620","AT5G42630")) %>%
  # For a correct ordering. Still did not managed to put grey bellow everybody
  mutate(type=factor(type, levels=c("None", "AT5G42580", "AT5G42590", "AT5G42600"))) %>%
  ggplot(data=., mapping=aes(x=row_idx, y=normalised, group=Gene)) +
  geom_line(aes(color = type), size=0.75)  +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  # grey for None put less impact.
  scale_color_manual(values=c("grey", "blue", "red", "green")) +
  ggtitle("normalised")


# Correlation plots
#
# Check different option there:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

library(corrplot)

# open a pdf file
pdf("out/RNAseq/corrplot.pdf", width=7, height=7)
GV_formated_data_marneral %>%
    # Only columns of interest
    dplyr::select(comp, value, Gene) %>%
    # need a wide format
    pivot_wider(names_from="Gene", values_from="value") %>%
    # format as data.frame with rownames (maybe not needed)
    as.data.frame %>%
    column_to_rownames("comp") %>%
    # Calculate the correlation coefficient between the different genes choose one from the list
    cor(method = "pearson") %>%
    # cor(method = "kendall") %>%
    # cor(method = "spearman") %>%
    # plot it nicely
    corrplot(title="raw")
# close the pdf file
dev.off()

GV_formated_data_marneral %>%
    # Only columns of interest
    dplyr::select(comp, scaled, Gene) %>%
    # need a wide format
    pivot_wider(names_from="Gene", values_from="scaled") %>%
    # format as data.frame with rownames (maybe not needed)
    as.data.frame %>%
    column_to_rownames("comp") %>%
    # Calculate the correlation coefficient between the different genes choose one from the list
    cor(method = "pearson") %>%
    # cor(method = "kendall") %>%
    # cor(method = "spearman") %>%
    # plot it nicely
    corrplot(title="scaled")


GV_formated_data_marneral %>%
    # Only columns of interest
    dplyr::select(comp, normalised, Gene) %>%
    # need a wide format
    pivot_wider(names_from="Gene", values_from="normalised") %>%
    # format as data.frame with rownames (maybe not needed)
    as.data.frame %>%
    column_to_rownames("comp") %>%
    # Calculate the correlation coefficient between the different genes choose one from the list
    cor(method = "pearson") %>%
    # cor(method = "kendall") %>%
    # cor(method = "spearman") %>%
    # plot it nicely
    corrplot(title="normalised", method="color", col = brewer.pal(n = 8, name = "RdYlBu"))




# Final graphs -----

Correlation_plot <- GV_formated_data_marneral %>%
  # Only columns of interest
  dplyr::select(comp, normalised, Gene) %>%
  # need a wide format
  pivot_wider(names_from="Gene", values_from="normalised") %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("comp") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


res1 <- GV_formated_data_marneral %>%
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

pdf("out/RNAseq/corrplot_GV.pdf", width=7, height=7)

corrplot(Correlation_plot, p.mat = res1$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


dev.off()



# Araport11 perturbation data -----
#data import

Araport11_gene_counts <- read_table2("data/RNAseq/Araport11_gene_counts.tsv")
Araport11_sample_info <- read_table2("data/RNAseq/Araport11_sample_info.tsv")
Araport11_SRA_details <- read_table2("data/RNAseq/Araport11_SRA_details.tsv")
Araport11_transcript_counts <- read_table2("data/RNAseq/Araport11_transcript_counts.tsv")

Araport11_formated_data_marneral <- Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, normalised_count) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup


Araport11_formated_data_MRN1_order <- Araport11_formated_data_marneral %>%
  filter(gene_id == "AT5G42600") %>% 
  arrange(scaled) %>%
  dplyr::select(sample_id) %>%
  mutate(row_idx=as.factor(row_number()))

Araport11_formated_data_marneral_order <-
  Araport11_formated_data_MRN1_order %>%
  left_join(Araport11_formated_data_marneral)

Araport11_formated_data_marneral_order$gene_id <- 
  factor(Araport11_formated_data_marneral_order$gene_id, levels=c("AT5G42500","AT5G42510","AT5G42520","AT5G42540","AT5G42560","AT5G42570",
                                                                  "AT5G42580","AT5G42590", "AT5G00580","AT5G06325","AT5G06335", "AT5G42600",
                                                                  "AT5G42620","AT5G42630","AT5G42640",
                                                                  "AT5G42650","AT5G42660","AT5G42670","AT5G42680","AT5G42690","AT5G42700"))



Correlation_plot_Araport <- Araport11_formated_data_marneral_order %>%
  # only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 



res2 <-  Araport11_formated_data_marneral_order %>%
  # only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

pdf("out/RNAseq/corrplot_Araport11.pdf", width=7, height=7)

corrplot(Correlation_plot_Araport, p.mat = res2$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)

corrplot(Correlation_plot_Araport, method="color", col = brewer.pal(n = 10, name = "RdYlBu"), pch.cex = .9, tl.col="black", 
         diag = FALSE)

dev.off()


# Araport11 only coding norm_count ------


Araport11_formated_data_marneral_coding <- Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, normalised_count) %>% 
  filter(gene_id != c("AT5G00580","AT5G06325","AT5G06335")) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup

Araport11_formated_data_MRN1_order_coding <- Araport11_formated_data_marneral %>%
  filter(gene_id == "AT5G42600") %>%
  arrange(scaled) %>%
  dplyr::select(sample_id) %>%
  mutate(row_idx=as.factor(row_number()))

Araport11_formated_data_marneral_order_coding <-
  Araport11_formated_data_MRN1_order %>%
  left_join(Araport11_formated_data_marneral)





Correlation_plot_Araport_coding <- Araport11_formated_data_marneral_order_coding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 



res2 <- Araport11_formated_data_marneral_order_coding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

corrplot(Correlation_plot_Araport_coding, p.mat = res2$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)

corrplot(Correlation_plot_Araport_coding, method="color", col = brewer.pal(n = 10, name = "RdYlBu"), pch.cex = .9, tl.col="black", 
         diag = FALSE)



# tpm without non-coding -----

Araport11_formated_data_marneral_coding <- Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, TPM) %>% 
  filter(! gene_id %in% c("AT5G00580","AT5G06325","AT5G06335")) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(TPM), # mean = 1, sd= 1 = mean / sd
         scaled=TPM / max(abs(TPM)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup



Correlation_plot_Araport_coding <- Araport11_formated_data_marneral_coding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


## specialized the insignificant value according to the significant level

corrplot(Correlation_plot_Araport_coding, method="color", col = brewer.pal(n = 10, name = "RdYlBu"), pch.cex = .9, tl.col="black", 
         diag = FALSE) 


# log(normalized + 1) without non-coding ------

Araport11_formated_data_marneral_coding <- Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, normalised_count) %>% 
  filter(! gene_id %in% c("AT5G00580","AT5G06325","AT5G06335")) %>%
  mutate(normalised_count=log(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup



Correlation_plot_Araport_coding <- Araport11_formated_data_marneral_coding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


## specialized the insignificant value according to the significant level

corrplot(Correlation_plot_Araport_coding, method="color", col = brewer.pal(n = 10, name = "RdYlBu"), pch.cex = .9, tl.col="black", 
         diag = FALSE) 



res2 <- Araport11_formated_data_marneral_coding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level

corrplot(Correlation_plot_Araport_coding, p.mat = res2$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)


# tpm with non-coding -----


Araport11_formated_data_marneral_noncoding <- Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, TPM) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(TPM), # mean = 1, sd= 1 = mean / sd
         scaled=TPM / max(abs(TPM)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup


Araport11_formated_data_marneral_noncoding$gene_id <- 
  factor(Araport11_formated_data_marneral_noncoding$gene_id, levels=c("AT5G42500","AT5G42510","AT5G42520","AT5G42540","AT5G42560","AT5G42570",
                                                                  "AT5G42580","AT5G42590", "AT5G00580","AT5G06325","AT5G06335", "AT5G42600",
                                                                  "AT5G42620","AT5G42630","AT5G42640",
                                                                  "AT5G42650","AT5G42660","AT5G42670","AT5G42680","AT5G42690","AT5G42700"))



Correlation_plot_Araport_noncoding <- Araport11_formated_data_marneral_noncoding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 


## specialized the insignificant value according to the significant level-OK

corrplot(Correlation_plot_Araport_noncoding, method="color", col = brewer.pal(n = 10, name = "RdYlBu"), pch.cex = .9, tl.col="black", 
         diag = FALSE) 



res2 <- Araport11_formated_data_marneral_noncoding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

## specialized the insignificant value according to the significant level



corrplot(Correlation_plot_Araport_noncoding, p.mat = res2$p, method="color", col = brewer.pal(n = 10, name = "RdYlBu"),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, tl.col="black", 
         diag = FALSE)



#log(normalized + 1) with non coding  -----
Araport11_formated_data_marneral_noncoding <-
  Araport11_gene_counts %>%
  dplyr::select(sample_id, gene_id, normalised_count) %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% c("AT5G42520","AT5G42540","AT5G42560","AT5G42570",
                        "AT5G42580","AT5G42590", "AT5G00580","AT5G06325","AT5G06335", "AT5G42600",
                        "AT5G42620","AT5G42630")) %>%
  mutate(gene_id = factor(gene_id,
                          levels=c("AT5G42520","AT5G42540","AT5G42560","AT5G42570",
                                   "AT5G42580","AT5G42590", "AT5G00580","AT5G06325","AT5G06335", "AT5G42600",
                                   "AT5G42620","AT5G42630"))

Correlation_plot_Araport_noncoding <-
  Araport11_formated_data_marneral_noncoding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor(method = "pearson") 

res2 <-
  Araport11_formated_data_marneral_noncoding %>%
  # Only columns of interest
  dplyr::select(sample_id, normalised, gene_id) %>%
  # need a wide format
  spread(gene_id, normalised) %>%
  # format as data.frame with rownames (maybe not needed)
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  # Calculate the correlation coefficient between the different genes choose one from the list
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/RNAseq/corrplot_Araport_log(norm_count).pdf", width=7, height=7)
    corrplot(Correlation_plot_Araport_noncoding, p.mat = res2$p, method="color",
             col = brewer.pal(n = 10, name = "RdYlBu"),
             insig = "label_sig", sig.level = c(.001, .01, .05),
             pch.cex = .9, tl.col="black", 
             diag = FALSE)
dev.off()

thalianol_cluster <-
    c("AT5G47880", "AT5G47890", "AT5G47900", "AT5G47910", "AT5G47920",
      "AT5G47930", "AT5G47940", "AT5G47950", "AT5G47960", "AT5G47970",
      "AT5G47980", "AT5G07035", "AT5G47990", "AT5G48000", "AT5G48010",
      "AT5G48020", "AT5G48030", "AT5G48040", "AT5G48050", "AT5G48060",
      "AT5G48070", "AT5G48080", "AT5G48090", "AT5G48100"
    )

arabidiol_cluster <-
    c( "AT4G15248", "AT4G15250", "AT4G06300", "AT4G06305", "AT4G06310",
      "AT4G06315", "AT4G06320", "AT4G15258", "AT4G15260", "AT4G15270",
      "AT4G15280", "AT4G15290", "AT4G15300", "AT4G15310", "AT4G15320",
      "AT4G06325", "AT4G15330", "AT4G15340", "AT4G15345", "AT4G15350",
      "AT4G15360", "AT4G15370", "AT4G15380", "AT4G15390", "AT4G15393",
      "AT4G15396", "AT4G15398", "AT4G15400", "AT4G15410", "AT4G15415",
      "AT4G15417", "AT4G15420", "AT4G15430", "AT4G15440", "AT4G15450")

tirucadiol_cluster <-
    c("AT5G36223", "AT5G36220", "AT5G36210", "AT5G36200", "AT5G36190",
      "AT5G36185", "AT5G36180", "AT5G36170", "AT5G36160", "AT5G05325",
      "AT5G36150", "AT5G36140", "AT5G36130", "AT5G36125", "AT5G36120",
      "AT5G36110", "AT5G36100", "AT5G36090", "AT5G36080", "AT5G36075",
      "AT5G36070", "AT5G36060", "AT5G36050", "AT5G36040", "AT5G36035")


Araport11_all_counts <- 
    read_tsv("data/RNAseq/Araport11_all_gene_counts.tsv.gz") %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "normalised_count")

# Thalianol
thalianol_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% thalianol_cluster) %>%
  mutate(gene_id = factor(gene_id, levels = thalianol_cluster))

thalianol_corr <-
  thalianol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

thalianol_sig <-
  thalianol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/RNAseq/corrplot_thalianol_Araport_log(norm_count).pdf", width=8, height=8)
    corrplot(thalianol_corr, p.mat = thalianol_sig$p, method="color",
             col = brewer.pal(n = 10, name = "RdYlBu"),
             insig = "label_sig", sig.level = c(.001, .01, .05),
             pch.cex = .9, tl.col="black", 
             diag = FALSE)
dev.off()




# Arabidiol
arabidiol_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% arabidiol_cluster) %>%
  mutate(gene_id = factor(gene_id, levels = arabidiol_cluster))

arabidiol_corr <-
  arabidiol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

arabidiol_sig <-
  arabidiol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/RNAseq/corrplot_arabidiol_Araport_log(norm_count).pdf", width=12, height=12)
    corrplot(arabidiol_corr, p.mat = arabidiol_sig$p, method="color",
             col = brewer.pal(n = 10, name = "RdYlBu"),
             insig = "label_sig", sig.level = c(.001, .01, .05),
             pch.cex = .9, tl.col="black", 
             diag = FALSE)
dev.off()

# Tirucadiol
tirucadiol_expr <-
  Araport11_all_counts %>%
  mutate(normalised_count=log2(normalised_count+1)) %>%
  group_by(gene_id) %>%
  mutate(normalised=scale(normalised_count), # mean = 1, sd= 1 = mean / sd
         scaled=normalised_count / max(abs(normalised_count)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup %>%
  filter(gene_id %in% tirucadiol_cluster) %>%
  mutate(gene_id = factor(gene_id, levels = tirucadiol_cluster))

tirucadiol_corr <-
  tirucadiol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor(method = "pearson") 

tirucadiol_sig <-
  tirucadiol_expr %>%
  dplyr::select(sample_id, normalised, gene_id) %>%
  spread(gene_id, normalised) %>%
  as.data.frame %>%
  column_to_rownames("sample_id") %>%
  cor.mtest(method = "pearson", conf.level = .95)

pdf("out/RNAseq/corrplot_tirucadiol_Araport_log(norm_count).pdf", width=10, height=10)
    corrplot(tirucadiol_corr, p.mat = tirucadiol_sig$p, method="color",
             col = brewer.pal(n = 10, name = "RdYlBu"),
             insig = "label_sig", sig.level = c(.001, .01, .05),
             pch.cex = .9, tl.col="black", 
             diag = FALSE)
dev.off()




























#### BONUS do not delet YET

# lncRNA MHAL is the most induced among lncRNA population?
all_info_noncoding <- read_excel("P:/RNAseq analysis/all_info_V201810415.xlsx") %>%
  select(Araport11_ID, type, width, NAT) %>%
  filter(Araport11_ID != "NA", 
         type == "noncoding",
         width>199, 
         NAT == "FALSE") %>%
  mutate(gene_id=Araport11_ID) %>%
  select(gene_id)


fresh_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/fresh_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

# MHAL 26/166 

dry_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/dry_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#45/116 

S1_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/S1_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#28/91 

S12_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/S12_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))


#8/130


S48_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/S48_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#53/113

SL1_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/SL1_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#23/96

SL6_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/SL6_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#8/79

SL12_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/SL12_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#10/69


SL24_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/SL24_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

#10/45

SL24_vs_SL48_down_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/SL24_vs_SL48_down_expressed.tsv") %>%
  select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange) %>%
  inner_join(all_info_noncoding) %>%
  filter(!duplicated(gene_id))

# Faire lncRNA avec voisins également Up !!



#clustering 

Global_time_diff_expressed <- 
  read.delim("K:/DDEVE (8095)/REGARN (8097)/Thomas/Projects/2019_Germination/TB_190923_01_DESeq2_FDR001/out/Global_time_diff_expressed.tsv") %>%
  select(gene_id) %>%
  transmute(ID=gene_id)


#
diffcounts <- normcounts %>% inner_join(Global_time_diff_expressed) 

diffcounts <- diffcounts %>% 
  mutate(fresh=rowMeans(diffcounts[,c(2:4)])) %>%
  mutate(dry=rowMeans(diffcounts[,c(5:7)])) %>%
  mutate(S1=rowMeans(diffcounts[,c(8:10)])) %>%
  mutate(S12=rowMeans(diffcounts[,c(11:13)])) %>%
  mutate(S48=rowMeans(diffcounts[,c(14:16)])) %>%
  mutate(SL1=rowMeans(diffcounts[,c(17:19)])) %>%
  mutate(SL6=rowMeans(diffcounts[,c(20:22)])) %>%
  mutate(SL12=rowMeans(diffcounts[,c(23:25)])) %>%
  mutate(SL24=rowMeans(diffcounts[,c(26:28)])) %>%
  mutate(SL48=rowMeans(diffcounts[,c(29:31)]))


diffcounts <- diffcounts %>%
  select(fresh, dry, S1, S12, S48, SL1, SL6, SL12, SL24, SL48, ID) %>%
  filter(!duplicated(ID)) 

#diffcounts <- diffcounts %>%
filter(col0 > 0)


diffcounts <- as.data.frame(diffcounts)



row.names(diffcounts)=diffcounts$ID
diffcounts$ID=NULL


library(pheatmap)

heatmap=pheatmap(diffcounts,scale = "row",cluster_cols = F )


cluster_gene=cbind(diffcounts,
                   cluster = cutree(heatmap$tree_row,
                                    k = 15))%>%add_rownames(var="ID")%>%arrange(desc(cluster))


hm=pheatmap::pheatmap(diffcounts[heatmap$tree_row[["order"]],],
                      cluster_cols = F,cluster_rows = T,
                      scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      fontsize_row = 2,cutree_rows = 15,
                      annotation_row =data.frame(
                        row.names = cluster_gene$ID,
                        cluster=paste0("cluster_",cluster_gene$cluster)))


#cluster selection ------
cl_1 <- cluster_gene %>% 
  filter(cluster == 1) %>%
  rename(ID = "code") %>%
  inner_join(clustering_DE)



write.table(cl_13,"cl_13.txt",sep="\t",row.names=FALSE)



