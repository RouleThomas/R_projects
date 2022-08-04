
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)




# Data import ----

GV <- read_excel("data/RNAseq/GV_solinc_region.xlsx", col_names=FALSE)


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


solinc_cluster <- data.frame(Gene = c("AT4G14548","AT4G14550","AT4G14560"),
                               type = c("AT4G14548","AT4G14550","AT4G14560"))

GV_formated_data_solinc <- GV_formated_data %>%
  filter(type == "Log2-ratio") %>%
  left_join(gene_info) %>%
  dplyr::select(comp, value, Gene) %>%
  mutate(value=as.numeric(value)) %>%
  left_join(solinc_cluster) %>% 
  mutate(type = fct_explicit_na(type, na_level = "None"), type=as.character(type)) %>%
  group_by(Gene) %>%
  mutate(normalised=scale(value), # mean = 1, sd= 1 = mean / sd
         scaled=value / max(abs(value)) # between -1 and +1, keeping the symetry
  ) %>%
  ungroup


GV_formated_data_solinc_order <- GV_formated_data_solinc %>%
  filter(Gene == "AT4G14548") %>% 
  arrange(value) %>%
  dplyr::select(comp) %>%
  mutate(row_idx=as.factor(row_number()))

GV_formated_data_solinc_order <-
  GV_formated_data_solinc %>%
  left_join(GV_formated_data_solinc_order)


# Graphs geom_bind2d

ggplot(GV_formated_data_solinc_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=value)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("raw")

ggplot(GV_formated_data_solinc_order, aes(row_idx, Gene)) + 
  geom_bin2d(aes(fill=scaled)) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_gradient2(low="#ca0020", mid="#f7f7f7", high="#0571b0") +
  ggtitle("scaled")

ggplot(GV_formated_data_solinc_order, aes(row_idx, Gene)) + 
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

Correlation_plot <- GV_formated_data_solinc %>%
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

