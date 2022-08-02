#theme graph
theme_set(theme_bw())

# Genotype metadata----------------------------------------------------------------
genotype_metadata <-
  tibble(genotype=c("col", "rnai92", "rnai38","rnai18", "pro", "lhp1", "clf", "35smrn1", "mrn1", "ubq","mro"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2","RNAi MARS 3", "SALK_133089", "lhp1", "clf", "35S:MRN1", "mrn1", "UBQ:MARS","mro"),
         genotype_color=c("#F8766D", "#C77CFF", "#00BFC4","#7CAE00", "#a6bddb", "#DDCC77", "#DDCC77", "#3f2d54", "#332288", "#AD1457","#003333")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))

genotype_metadata_root <-
  tibble(genotype=c("col", "rnai92", "rnai38","rnai18", "pro","35smrn1", "mrn1", "mro"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2","RNAi MARS 3", "SALK_133089", "35S:MRN1", "mrn1", "mro"),
         genotype_color=c("#F8766D", "#C77CFF", "#00BFC4","#7CAE00", "#a6bddb", "#004529", "#238443", "#78c679")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))

genotype_metadata_germination <-
  tibble(genotype=c("col", "rnai92", "rnai38", "pro", "35smrn1", "mrn1"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2", "SALK_133089", "35S:MRN1", "mrn1"),
         genotype_color=c("#F8766D", "#C77CFF", "#00BFC4","#7CAE00", "#3f2d54", "#336666")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))


#VERSION 2 COLOR -----
genotype_metadata <-
  tibble(genotype=c("col", "rnai92", "rnai38","rnai18", "pro", "lhp1", "clf", "35smrn1", "mrn1", "ubq","mro"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2","RNAi MARS 3", "SALK_133089", "lhp1", "clf", "35S:MRN1", "mrn1", "UBQ:MARS","mro"),
         genotype_color=c("#636363", "#023858", "#0570c0", "#3690c0", "#a6bddb", "#DDCC77", "#DDCC77", "#3f2d54", "#332288", "#AD1457","#003333")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))

genotype_metadata_root <-
  tibble(genotype=c("col", "rnai92", "rnai38","rnai18", "pro","35smrn1", "mrn1", "mro"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2","RNAi MARS 3", "SALK_133089", "35S:MRN1", "mrn1", "mro"),
         genotype_color=c("#636363", "#023858", "#0570c0","#3690c0", "#a6bddb", "#004529", "#238443", "#78c679")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))


genotype_metadata_ubq <-
  tibble(genotype=c("col", "ubq"),
         genotype_name=c("Col-0", "UBQ:MARS"),
         genotype_color=c("#636363","#f03b20")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))

genotype_metadata_germination <-
  tibble(genotype=c("col", "rnai92", "rnai38", "pro", "35smrn1", "mrn1"),
         genotype_name=c("Col-0", "RNAi MARS 1", "RNAi MARS 2", "SALK_133089", "35S:MRN1", "mrn1"),
         genotype_color=c("#636363", "#023858", "#3690c0","#a6bddb", "#004529", "#238443")) %>%
  # Transform genotype name as a factor keeping the order
  mutate(genotype_name=factor(genotype_name, levels=genotype_name))


# Gene metadata---------------------------------------------------------------------
selected_gene_RIP <-
  tibble(gene=c("ref1", "MRN1", "APOLO","MHAL","MHAL.3'"),
         gene_color=c("#636363", "#2CA02C", "#87CEEB", "#D62728","#D62728"))

selected_gene <-
  tibble(gene=c("CYP705A12", "CYP71A16", "MHAL", "MRN1", "MHAL_1","MHAL_2","MHAL_3","MHAL_4","ASCO","APOLO","COLDAIR"),
         gene_name=c("CYP705A12", "CYP71A16", "MHAL", "MRN1", "MHAL.1","MHAL.2","MHAL.3","MHAL.4","ASCO","APOLO","COLDAIR"),
         gene_color=c("#1F77B4", "#9467BD", "#D62728","#2CA02C","#D62728","#D62728","#D62728","#D62728","#87CEEB","#87CEEB","#87CEEB"))

selected_gene_qPCR <-
  tibble(gene=c("CYP705A12", "CYP71A16", "MHAL", "MRN1","RAB18","RD29B"),
         gene_name=c("CYP705A12", "CYP71A16", "MHAL", "MRN1","RAB18","RD29B"),
         gene_color=c("#1F77B4", "#9467BD", "#2CA02C","#D62728", "#999933","#DDCC77"))


selected_gene2_qPCR <-
  tibble(gene=c("CYP705A12", "CYP71A16", "MHAL", "MRN1","RAB18","RD29B"),
         gene_name=c("CYP705A12", "CYP71A16", "MHAL", "MRN1","RAB18","RD29B"),
         gene_color=c("#1F77B4", "#9467BD", "#D62728","#2CA02C", "#999933","#DDCC77"))

selected_gene_RNAseq <-
  tibble(gene=c("CYP705A12", "CYP71A16", "MHAL", "MRN1"),
         gene_name=c("CYP705A12", "CYP71A16", "MHAL", "MRN1"),
         gene_color=c("#1F77B4", "#9467BD", "#D62728","#2CA02C"))

selected_gene_CPC <-
  tibble(gene=c("CYP705A12", "CYP71A16","MRN1", "AT5G00580.1","AT5G00580.2","AT5G00580.3","AT5G00580.4","AT5G06335","AT5G06325", "ASCO","APOLO","COLDAIR"),
         gene_color=c("#1F77B4", "#9467BD", "#2CA02C","#D62728","#D62728","#D62728","#D62728","#D62728","#D62728", "#87CEEB","#87CEEB","#87CEEB"))

selected_gene_nucleus <-
  tibble(gene=c("ref1", "APOLO","ASCO", "MHAL","U6"),
         gene_name=c("ref1", "APOLO", "ASCO", "MHAL","U6"),
         gene_color=c("#636363", "#87CEEB","#87CEEB","#D62728","#87CEEB"))
