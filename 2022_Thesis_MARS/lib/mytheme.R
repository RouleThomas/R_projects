#Packages
library("tidyverse")
library("readxl")
library("readr")
library("ggpubr")
library("psych")
library("multcompView")
library("RColorBrewer")
library("lubridate")

#theme graph
theme_set(theme_bw())

# Genotype metadata----------------------------------------------------------------
selected_genotype_root <-
  tibble(genotype=c("col", "rnai38", "rnai92", "pro","35smrn1"),
         genotype_name=c("Col-0", "RNAi 3.8", "RNAi 9.2", "SALK_133089", "35S:MRN1"),
         genotype_color=c("#F8766D", "#00BFC4", "#C77CFF","#7CAE00", "#3f2d54")) 


selected_genotype <-
  tibble(genotype=c("col", "RNAi3.8", "RNAi9.2", "pro"),
         genotype_name=c("Col-0", "RNAi 3.8", "RNAi 9.2", "SALK"),
         genotype_color=c("#F8766D", "#7CAE00", "#00BFC4","#C77CFF"))

selected_genotype_3C <-
  tibble(genotype=c("col", "rnai38", "rnai92"),
         genotype_name=c("Col-0", "RNAi 3.8", "RNAi 9.2"),
         genotype_color=c("#636363", "#023858", "#0570c0"))

selected_genotype_CHIP <-
  tibble(genotype=c("Col-0", "RNAi 9.2"),
         genotype_color=c("#F8766D", "#C77CFF"))

selected_genotype_lhp1 <-
  tibble(genotype=c("Col-0", "lhp1"),
         genotype_color=c("#F8766D", "#DDCC77"))

selected_genotype_clf <-
  tibble(genotype=c("Col-0", "clf"),
         genotype_color=c("#DDCC77", "#F8766D"))

selected_genotype_clf2 <-
  tibble(genotype=c("Col-0", "clf"),
         genotype_color=c("#F8766D", "#DDCC77"))

# Gene metadata---------------------------------------------------------------------
selected_gene_RIP <-
  tibble(gene=c("ref1", "MRN1", "APOLO","MHAL.5'","MHAL.3'"),
         gene_color=c("#636363", "#2CA02C", "#87CEEB", "#D62728","#D62728"))

selected_gene_meDIP <-
  tibble(gene=c("APOLO.5'","MARS.5'","MARS.3'"),
         gene_color=c("#87CEEB", "#D62728","#D62728"))

selected_gene <-
  tibble(gene=c("CYP705A12", "CYP71A16", "MHAL", "MRN1", "MHAL_1","MHAL_2","MHAL_3","MHAL_4","ASCO","APOLO","COLDAIR"),
         gene_name=c("CYP705A12", "CYP71A16", "MHAL", "MRN1", "MHAL.1","MHAL.2","MHAL.3","MHAL.4","ASCO","APOLO","COLDAIR"),
         gene_color=c("#1F77B4", "#9467BD", "#D62728","#2CA02C","#D62728","#D62728","#D62728","#D62728","#87CEEB","#87CEEB","#87CEEB"))


selected_gene_CPC <-
  tibble(gene=c("CYP705A12", "CYP71A16","MRN1", "MHAL.1","MHAL.2","MHAL.3","MHAL.4","ASCO","APOLO","COLDAIR"),
         gene_name=c("CYP705A12", "CYP71A16", "MRN1", "MHAL.1","MHAL.2","MHAL.3","MHAL.4","ASCO","APOLO","COLDAIR"),
         gene_color=c("#1F77B4", "#9467BD", "#2CA02C","#D62728","#D62728","#D62728","#D62728","#87CEEB","#87CEEB","#87CEEB"))

selected_gene_nucleus <-
  tibble(gene=c("ref1", "APOLO","ASCO", "MHAL","U6"),
         gene_name=c("ref1", "APOLO", "ASCO", "MHAL.5'","U6"),
         gene_color=c("#636363", "#87CEEB","#87CEEB","#D62728","#87CEEB"))

selected_gene_invitro <-
  tibble(RNA=c("GFP", "MARS"),
         gene_color=c("#636363", "#D62728"))

