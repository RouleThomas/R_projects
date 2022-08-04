# MHAL and colateral genes expression in Col0/RNAi MARS 1 and 2/T-DNA/CRISPR lines during ABA treatment_0/4H
# Aur?lie
# 06/08/2021
#


library(tidyverse)
library(readxl)
library(multcompView)
library(PMCMRplus)



load_sample_sheet <- function(excel_file){
  sheets <- excel_sheets(excel_file)
  sample_sheet <-
    read_xlsx(excel_file) %>%
    dplyr::rename(row=`...1`) %>%
    pivot_longer(-row, names_to="column", values_to="values") %>%
    unite(Pos, c(row, column), sep="") %>%
    dplyr::select(-values)
  for (sheet_name in sheets) {
      print(str_c("processing sheet ", sheet_name))
    sheet_info <-
      read_xlsx(excel_file, sheet = sheet_name) %>%
      dplyr::rename(row=`...1`) %>%
      mutate_all(as.character) %>%
      pivot_longer(-row, names_to="column", values_to=sheet_name) %>%
      unite(Pos, c(row, column), sep="")
    sample_sheet <-
      left_join(sample_sheet, sheet_info)
  }
  return(sample_sheet)
}

#chargement du/des plan(s) de plaque 
plan_plaque_1 <-
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/plan_plaque_RNA_MARS_CRISPR_kinetic_1.xlsx")



# Chargement des Cq 
qPCR_1 <-
  read_xlsx("data/qPCR_cDNA/MARS CRISPR/MARS crispr1.xlsx", na="N/A") %>%
  dplyr::rename(Pos="Well")

donnees_1 <-
  qPCR_1 %>%
  left_join(plan_plaque_1, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)



ref_gene <- "ref2"
ref_data <- donnees_1 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- donnees_1 %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "WC",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

ddCt_data_crispr <- ddCt_data


stat_crispr <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_crispr$genotype <- factor(stat_crispr$genotype, levels=c("WC","CRISPR"))






#time 0 ----

stat_qPCR_cDNA <- ddCt_data %>%
  filter(time =="0") %>%
  group_by(gene, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            sem=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(time =="0", gene !="MHAL") %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("WC","CRISPR"))


my_graph <- 
  stat_qPCR_cDNA %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=genotype, y=-mean)) +
  geom_bar(stat="identity", position=position_dodge(0.9), aes(fill=gene)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



my_graph <- 
  stat_qPCR_cDNA %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  left_join(selected_gene2_qPCR)%>%
  ggplot(data=., mapping=aes(x=genotype, y=mean)) +
  geom_bar(stat="identity", position=position_dodge(0.9), aes(fill=gene)) +
  geom_errorbar(mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/CRISPR/CRISPR_t0.pdf", plot=my_graph, width = 7, height = 3)













my_graph <- 
stat_crispr %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_wrap(~genotype)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )  + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/CRISPR/CRISPR_kinetics_not_norm.pdf", plot=my_graph, width = 8, height = 3.5)


dCt_data_control <- dCt_data %>%
  filter(genotype == "WC",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_WC <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))%>%
  filter(genotype=="WC")

dCt_data_control <- dCt_data %>%
  filter(genotype == "CRISPR",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data_CRISPR <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))%>%
  filter(genotype=="CRISPR")

ddCt_data_crispr <- ddCt_data_WC %>%
  bind_rows(ddCt_data_CRISPR)


stat_crispr <- ddCt_data_crispr%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_crispr$genotype <- factor(stat_crispr$genotype, levels=c("WC","CRISPR"))


my_graph <- 
  stat_crispr %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_wrap(~genotype)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )  + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/CRISPR/CRISPR_kinetics_norm.pdf", plot=my_graph, width = 5, height = 3.5)




# ANOVA ------

anova_MRN1 <- ddCt_data_crispr %>% 
  filter(gene == "MRN1") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="MRN1") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP705A12 <- ddCt_data_crispr %>% 
  filter(gene == "CYP705A12") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP705A12") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_CYP71A16 <- ddCt_data_crispr %>% 
  filter(gene == "CYP71A16") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="CYP71A16") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)


anova_all_gene <- anova_MRN1 %>%
  bind_rows(anova_CYP705A12,anova_CYP71A16) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble() 

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene


my_graph <- 
  anova_all_gene %>% filter(!gene %in% c("MHAL")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=-diff + sign(-diff + 0.0001))) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph



# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_CRISPR_marneral.pdf", plot=my_graph, width = 3.5, height = 3.5)







# compil with the RNAi and T-DNA; before run qPCR_cDNA2 to get "stat_rnai" --------------

rnai_crispr <- stat_rnai %>%	
  filter(gene%in% c("CYP705A12","CYP71A16","MRN1","MHAL")) %>%
  bind_rows(stat_crispr%>%filter(genotype!="WC", gene!="MARS")) %>%
  filter(genotype !="rnai18")

rnai_crispr$genotype <-   factor(rnai_crispr$genotype, c("col","rnai92","rnai38","pro","CRISPR"))

my_graph <- 
rnai_crispr %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_crispr.pdf", plot=my_graph, width = 9, height = 4)




my_graph <- 
  rnai_crispr %>% filter(time<1.5) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/re_all_marneral_crispr_zoom.pdf", plot=my_graph, width = 9, height = 5)








# control condition

ddCt_data_rnai_process <- ddCt_data_rnai %>% 
  filter(genotype !="rnai18",
         gene%in%c("CYP705A12","CYP71A16","MHAL","MRN1")) 

ddCt_data_crispr_process <- ddCt_data_crispr %>% 
  filter(genotype =="CRISPR",
         gene%in%c("CYP705A12","CYP71A16","MRN1"))

ddCt_data_crispr_rnai <- ddCt_data_rnai_process %>%
  bind_rows(ddCt_data_crispr_process)


stat_qPCR_cDNA <-ddCt_data_crispr_rnai %>%
  filter(time =="0") %>%
  group_by(gene, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data_crispr_rnai %>%
  filter(time =="0") %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <-   factor(stat_qPCR_cDNA$genotype, c("col","rnai92","rnai38","pro","CRISPR"))


my_graph <- 
  stat_qPCR_cDNA  %>%
  ggplot(data=., mapping=aes(x=genotype, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/marneral cluster control condition_crispr_gene.pdf", plot=my_graph, width = 7, height = 3.5)












# Aurelie root shoot and TRANS target-------------

#chargement du/des plan(s) de plaque 



plan_plaque_1 <-
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/root shoot/plan_plaque_20210805.xlsx")

plan_plaque_2 <-
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/root shoot/plan_plaque_root_shoot_1.xlsx")

plan_plaque_3 <-
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/root shoot/plan_plaque_root_shoot_2.xlsx")

plan_plaque_1_trans <- 
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/target trans/plan_plaque_root_shoot_TRANS1.xlsx")
plan_plaque_2_trans <- 
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/target trans/plan_plaque_root_shoot_TRANS2.xlsx")
plan_plaque_3_trans <- 
  load_sample_sheet("data/qPCR_cDNA/MARS CRISPR/target trans/plan_plaque_root_shoot_TRANS3.xlsx")

# Chargement des Cq 
qPCR_1 <-
  read_xlsx("data/qPCR_cDNA/MARS CRISPR/root shoot/20210805_qPCR1.xlsx", na="N/A") %>%
  dplyr::rename(Pos="Well")

qPCR_2 <-
  read_xlsx("data/qPCR_cDNA/MARS CRISPR/root shoot/root_shoot_1.xlsx", na="N/A") %>%
  dplyr::rename(Pos="Well")

qPCR_3 <-
  read_xlsx("data/qPCR_cDNA/MARS CRISPR/root shoot/root_shoot_2.xlsx", na="N/A") %>%
  dplyr::rename(Pos="Well")

qPCR_1_trans <-
  read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/root_shoot_trans_qPCR_output.xlsx", na="N/A",sheet=1)
qPCR_2_trans <-
  read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/root_shoot_trans_qPCR_output.xlsx", na="N/A",sheet=2)
qPCR_3_trans <-
  read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/root_shoot_trans_qPCR_output.xlsx", na="N/A",sheet=3)



donnees_1 <-
  qPCR_1 %>%
  left_join(plan_plaque_1, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)

donnees_2 <-
  qPCR_2 %>%
  left_join(plan_plaque_2, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)

donnees_3 <-
  qPCR_3 %>%
  left_join(plan_plaque_3, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)

donnees_trans1 <-
  qPCR_1_trans %>%
  left_join(plan_plaque_1_trans, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)
donnees_trans2 <-
  qPCR_2_trans %>%
  left_join(plan_plaque_2_trans, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)
donnees_trans3 <-
  qPCR_3_trans %>%
  left_join(plan_plaque_3_trans, by = c("Pos")) %>%
  filter(!is.na(gene),
         !is.na(Cp),
         !is.na(genotype),
         !is.na(rep_bio),
         !is.na(time)) %>%
  dplyr::select(-Pos)


donnees <- donnees_1 %>%
  bind_rows(donnees_2) %>%
  bind_rows(donnees_3)
donnees <- donnees_1 %>%
  bind_rows(donnees_trans1)%>%
  bind_rows(donnees_trans2) %>%
  bind_rows(donnees_trans3)

donnees <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/trans_root_shoot_qPCR_output.xlsx")



ref_gene <- "ref2"
ref_data <- donnees %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- donnees %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "WC",
         time==0,
         organ =="Root") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_root_shoot <- ddCt_data%>%			
  group_by(gene, time, genotype,organ) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            sem=SD/sqrt(n)) 		



stat_root_shoot$genotype <- factor(stat_root_shoot$genotype, levels=c("WT","RNAi1","RNAi2","TDNA","WC","CRISPR"))

my_graph <- 
stat_root_shoot %>%
  filter(!(gene %in% c("MHAL","MHAL isof1","MHAL isof2", "MHAL isof3","MHAL isof4")))%>%
  filter(genotype %in% c("WT","TDNA","RNAi1","RNAi2"),
         organ=="Shoot")%>%
  ggplot(data=.,	 mapping=aes(x=as.character(time),	 y=mean,	 group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +		
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_wrap(~gene, scale="free")

#my_graph



#ggsave(filename="out/qPCR_cDNA_2/CRISPR/trans/RNAi_tdna_trans.pdf", plot=my_graph, width = 6, height = 6)


my_graph <- 
  stat_root_shoot %>%
  filter(!(gene %in% c("MHAL","MHAL isof1","MHAL isof2", "MHAL isof3","MHAL isof4")))%>%
  filter(genotype %in% c("WC","CRISPR"),
         organ=="Root")%>%
  ggplot(data=.,	 mapping=aes(x=as.character(time),	 y=mean,	 group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +		
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_wrap(~gene, scale="free")

my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/trans/CRISPR_trans.pdf", plot=my_graph, width = 6, height = 6)











stat_root_shoot %>%
  filter(gene %in% c("CYP705A12","CYP71A16","MRN1","MHAL"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(organ~genotype)


stat_root_shoot %>%
  filter(gene %in% c("MHAL","MHAL isof1","MHAL isof2", "MHAL isof3","MHAL isof4"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(organ~genotype)


stat_root_shoot %>%
  filter(!(gene %in% c("MHAL","MHAL isof1","MHAL isof2", "MHAL isof3","MHAL isof4")))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - sem,	 ymax=mean + sem),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(organ~genotype)

stat_root_shoot %>%
  filter(!(gene %in% c("MHAL","MHAL isof1","MHAL isof2", "MHAL isof3","MHAL isof4")))%>%
  filter(genotype %in% c("WT","TDNA","RNAi1","RNAi2"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=genotype)) +
  geom_line(aes(color = genotype),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - sem,	 ymax=mean + sem),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(organ~gene)






# PLOT Root vs shoot in WT -----


stat_root_shoot$organ <- factor(stat_root_shoot$organ, levels=c("Shoot","Root"))


my_graph <- 
  stat_root_shoot %>%
  filter(genotype =="WT", time =="0",
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1"))%>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=organ, y=mean, group=organ)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=gene)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  facet_wrap(~gene,nrow=1)+
  scale_fill_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 
my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/WT_ROOTvsSHOOT.pdf", plot=my_graph, width = 6, height = 3)





## STAT -----
# Root or Shoot ----

##STAT entre time je prefere
stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(!gene %in% c("RAB18","RD29B"),
                          gene%in% c("MRN1","CYP705A12","CYP71A16","MHAL")) %>%
  mutate(time=as.character(time))%>%
  filter(organ =="Shoot") %>%
  dplyr::select(-organ)%>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "time")) %>%
  dplyr::select(-data) %>%
  unnest(stat)


##STAT entre genotype
stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(!gene %in% c("RAB18","RD29B"),
         genotype %in% c("WT","RNAi1","RNAi2","TDNA")) %>%
  mutate(time=as.character(time))%>%
  filter(organ =="Root") %>%
  select(-organ)%>%
  group_by(gene,time) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)



#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("WT","RNAi1","RNAi2","TDNA","WC","CRISPR"))


my_graph <- 
  stat_qPCR_cDNA %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=time, y=-mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



genotype_corr <- data.frame (genotype  = c("WT", "RNAi1", "RNAi2","TDNA"),
                  genotype_name = c("Col-0", "RNAi MARS 1", "RNAi MARS 2","mrs1-1")
)




my_graph <- 
  stat_qPCR_cDNA %>% filter(!gene %in% c("RAB18","RD29B"),
                            genotype%in% c("WT","RNAi1","RNAi2","TDNA")) %>%
  mutate(time=as.character(time))%>%
  left_join(genotype_corr)%>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean, group=time)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=time)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_grey()
my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/Shoot_cluster.pdf", plot=my_graph, width = 6, height = 3)
ggsave(filename="out/qPCR_cDNA_2/CRISPR/Root_cluster.pdf", plot=my_graph, width = 6, height = 3)






stat_root_shoot %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(organ~genotype)



# voisin1 CRISPR ------


voisin1 <- read_excel("data/qPCR_cDNA/MARS CRISPR/voisin_mars/voisin1.xlsx") %>%
  dplyr::select(-Position) %>%
  mutate(Cp=as.numeric(Cp))



ref_gene <- "ref2"
ref_data <- voisin1 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- voisin1 %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "WC",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_voisin <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_voisin %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(~genotype)



## STAT -----

##STAT entre time je prefere
stat_qPCR_cDNA <- 
  ddCt_data%>%
  mutate(time=as.character(time))%>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "time")) %>%
  dplyr::select(-data) %>%
  unnest(stat)


##STAT entre genotype
stat_qPCR_cDNA <- 
  ddCt_data %>%
  mutate(time=as.character(time))%>%
  group_by(gene,time) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)



#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("WC","CRISPR"))


my_graph <- 
  stat_qPCR_cDNA %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=time, y=-mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/voisin_CRISPR.pdf", plot=my_graph, width = 6, height = 3)



# voisin2 RNAi ------


voisin2 <- read_excel("data/qPCR_cDNA/MARS CRISPR/voisin_mars/voisin2.xlsx") %>%
  dplyr::select(-Pos) %>%
  mutate(Cp=as.numeric(Cp))



ref_gene <- "ref2"
ref_data <- voisin2 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- voisin2 %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "Col0",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_voisin <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_voisin %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(~genotype)



## STAT -----

##STAT entre time 
stat_qPCR_cDNA <- 
  ddCt_data%>%
  mutate(time=as.character(time))%>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "time")) %>%
  dplyr::select(-data) %>%
  unnest(stat)


##STAT entre genotype je prefere
stat_qPCR_cDNA <- 
  ddCt_data %>%
  mutate(time=as.character(time))%>%
  group_by(gene,time) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)



#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("Col0","RNAi 1","RNAi 2", "mrs11"))


my_graph <- 
  stat_qPCR_cDNA %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=time, y=-mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/voisin_RNAi.pdf", plot=my_graph, width = 6, height = 3)



# trans candidats ------


trans_1 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/trans-1.xlsx") %>%
  select(-Position,-ID)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))
trans_3 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/trans-3.xlsx") %>%
  select(-Position)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))
trans <- trans_1 %>%
  bind_rows(trans_3)

ref_gene <- "ref"
ref_data <- trans %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- trans %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "col",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))


stat_trans1 <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_trans1$gene <- factor(stat_trans1$gene, levels=c("GGT2","ABC1K3","NIR1","PPC2","QQS"))


stat_trans1 %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+
  facet_wrap(~gene,nrow=1)+
  theme_bw()


# RNAi


trans_2 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/trans-2.xlsx") %>%
  select(-Pos,-ID)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))
trans_4 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/trans-4.xlsx") %>%
  select(-Pos)%>%
  drop_na()%>%
  mutate(Cp=as.numeric(Cp))

trans_RNAi <- trans_2 %>%
  bind_rows(trans_4)


ref_gene <- "ref"
ref_data <- trans_RNAi %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- trans_RNAi %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "col",
         time==0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))


stat_trans2 <- ddCt_data%>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 	

stat_trans2$gene <- factor(stat_trans2$gene, levels=c("GGT2","ABC1K3","NIR1","PPC2","QQS"))


stat_trans2 %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+
  facet_wrap(~gene,nrow=1)



# crispr kinetics isoforms and ref genes ------------



qPCR_384_crispr <- read_excel("data/qPCR_cDNA/MARS CRISPR/qPCR_384_crispr.xlsx")%>%
  drop_na() %>%
  dplyr::select(-Pos)


ref_gene <- "ref"
ref_data <- qPCR_384_crispr %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_384_crispr %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "WC") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_crispr <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


stat_crispr$genotype <- factor(stat_crispr$genotype, c("WC","CRISPR"))

my_graph <- 
stat_crispr %>%	filter(gene%in%c("iso1","iso2","iso3"),genotype=="WC")%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw()+			
  facet_wrap(~genotype,	 nrow=1)
my_graph

# Save  
ggsave(filename="out/qPCR_cDNA_2/CRISPR/isoforms_kin.pdf", plot=my_graph, width = 5, height = 3.5)




my_graph <- 
stat_crispr %>%	filter(gene%in%c("RAB18","RD29B"))%>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw()+			
  facet_wrap(~genotype,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) + theme(legend.position = "none")
my_graph

# Save  
ggsave(filename="out/qPCR_cDNA_2/CRISPR/kin_ABA_genes.pdf", plot=my_graph, width = 5, height = 3.5)


# ANOVA ------

anova_RAB18 <- ddCt_data %>% 
  filter(gene == "RAB18") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RAB18") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)
anova_RD29B <- ddCt_data %>% 
  filter(gene == "RD29B") %>% 
  aov(-ddCt ~ genotype + factor(time), data=.) %>% 
  TukeyHSD("genotype") %>%
  .$genotype %>%
  as.data.frame() %>%
  dplyr::select(diff, "p adj") %>%
  rownames_to_column() %>%
  mutate(gene ="RD29B") %>%
  separate(rowname, sep = "-", c("genotype","comp")) %>%
  dplyr::select(-comp)

anova_all_gene <- anova_RAB18 %>%
  bind_rows(anova_RD29B) %>%
  rename(pval="p adj") %>%
  mutate(pval=as.numeric(pval)) %>%
  as.tibble() 

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

anova_all_gene$pval <- makeStars(anova_all_gene$pval)
anova_all_gene



my_graph <- 
  anova_all_gene %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=gene, y=-diff, fill=genotype)) +
  geom_bar(stat="identity", position=position) +
  geom_text(position=position, size = 2, aes(label=pval%>%formatC(format = "e", digits = 2), y=-diff - 0.1)) +
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  xlab(label = "") +
  ylab(label = "") +
  theme(strip.text = element_text(size = rel(1)))+
  theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/anova_CRISPR_ABA_genes.pdf", plot=my_graph, width = 3.5, height = 3.5)




# ROOT shoot CRISPR repet -------------

root_shoot_crispr_qPCRoutput_repet <- read_excel("data/qPCR_cDNA/MARS CRISPR/root shoot/root_shoot_crispr_qPCRoutput_repet.xlsx",sheet=1)%>%
  drop_na()%>%
  dplyr::select(-Pos)


ref_gene <- "ref"
ref_data <- root_shoot_crispr_qPCRoutput_repet %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- root_shoot_crispr_qPCRoutput_repet %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)%>%
  drop_na()

dCt_data_control <- dCt_data %>%
  filter(genotype == "WT",
         time==0,
         organ =="Root") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


stat_root_shoot <- ddCt_data%>%			
  group_by(gene, time, genotype,organ) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            sem=SD/sqrt(n)) 			


stat_root_shoot$genotype <- factor(stat_root_shoot$genotype, levels=c("WT","CRISPR"))




my_graph <- 
  stat_root_shoot %>% filter(!gene %in% c("RAB18","RD29B"),
                             genotype%in% c("WT","CRISPR"),
                             organ=="Root") %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=genotype, y=mean, group=time)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=time)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - sem), ymax=(mean + sem), width=.5))+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_grey()
my_graph


ggsave(filename="out/qPCR_cDNA_2/CRISPR/Shoot_CRISPR_cluster.pdf", plot=my_graph, width = 4, height =2.5)
ggsave(filename="out/qPCR_cDNA_2/CRISPR/Root_CRISPR_cluster.pdf", plot=my_graph, width = 4, height = 2.5)






































## STAT -----
# Root or Shoot ----

##STAT entre time je prefere
stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(!gene %in% c("RAB18","RD29B"),
         gene%in% c("MRN1","CYP705A12","CYP71A16")) %>%
  mutate(time=as.character(time))%>%
  filter(organ =="Root") %>%
  dplyr::select(-organ)%>%
  group_by(gene,genotype) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "time")) %>%
  dplyr::select(-data) %>%
  unnest(stat)


##STAT entre genotype
stat_qPCR_cDNA <- 
  ddCt_data %>%
  filter(!gene %in% c("RAB18","RD29B","MARS"),
         genotype %in% c("WT","CRISPR")) %>%
  mutate(time=as.character(time))%>%
  filter(organ =="Root") %>%
  dplyr::select(-organ)%>%
  group_by(gene,time) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)



#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, levels=c("WT","RNAi1","RNAi2","TDNA","WC","CRISPR"))


my_graph <- 
  stat_qPCR_cDNA %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=time, y=-mean, group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph




my_graph <- 
  stat_qPCR_cDNA %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  mutate(time=as.character(time))%>%
  ggplot(data=., mapping=aes(x=genotype, y=-mean, group=time)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=time)) +
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter, y=-mean - sign(mean - 0.0001) * (sem + .6)), size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  theme(strip.text = element_text(size = rel(1)))
my_graph



ggsave(filename="out/qPCR_cDNA_2/CRISPR/Shoot_cluster.pdf", plot=my_graph, width = 6, height = 3)






# MRNO kin TRANS target -----

MRO_output_1 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/MRO_output_1.xlsx")%>%
  select(-Pos)


ref_gene <- "ref"
ref_data <- MRO_output_1 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- MRO_output_1 %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_mro <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			





my_graph <- 
  stat_mro %>%
  filter(time %in% c(0,4), genotype !="RNAi")%>%
  ggplot(data=.,	 mapping=aes(x=as.character(time),	 y=mean,	 group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +		
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(~gene,scale="free")

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/CRISPR/trans/MRO_trans.pdf", plot=my_graph, width = 7, height =2)


# MRN1 kin TRANS target -----

MRN_output_1 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/MRN1_output.xlsx",sheet=1)%>%
  select(-Pos)
MRN_output_2 <- read_excel("data/qPCR_cDNA/MARS CRISPR/target trans/MRN1_output.xlsx",sheet=2)%>%
  select(-Pos)

MRN_output <- MRN_output_1 %>%
  bind_rows(MRN_output_2) %>%
  drop_na()

ref_gene <- "ref"
ref_data <- MRN_output %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- MRN_output %>%
  filter(gene != ref_gene) 

dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_mrn <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


stat_mrn$genotype <- factor(stat_mrn$genotype, levels=c("col","35S:MRN","mrn1"))



my_graph <-  
  stat_mrn %>% filter(time %in% c(0,8), genotype !="RNAi9.2") %>%
  ggplot(data=.,	 mapping=aes(x=as.character(time),	 y=mean,	 group=genotype)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), aes(fill=genotype)) +		
  geom_errorbar(position=position_dodge(width=0.9),mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std), width=.5))+	
  geom_hline(yintercept=0) +
  geom_hline(yintercept=1, linetype=2) +	
  geom_hline(yintercept=-1, linetype=2) +	
  theme_bw() +
  facet_grid(~gene)

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA_2/CRISPR/trans/MRN_trans.pdf", plot=my_graph, width = 7, height = 2)













