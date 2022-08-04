# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/qPCR_cDNA/", showWarnings = FALSE, recursive = TRUE)

#qPCR analysis of Col and MHAL deregulated lines

#ABA treatment----
##Normalization time0 col-0 -----

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1"),
         time == 0) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V2.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92","rnai18"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1"),
         time == 0) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V3.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1"),
         time == 0) %>%
  clean_genotype_code

position=position_dodge(.10)

stat_qPCR_cDNA <- qPCR_cDNA_ABA %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  qPCR_cDNA_ABA %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster control condition_gene.pdf", plot=my_graph, width = 7, height = 3.5)

##Normalization time0 each genotype ----

#
qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA.xlsx",sheet=1) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V2.xlsx",sheet=1) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92","rnai18"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V3.xlsx",sheet=1) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","pro","rnai38","rnai92"),
         gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>%
  clean_genotype_code




qPCR_cDNA_ABA$time <- as.POSIXct(strptime(qPCR_cDNA_ABA$time, format="%H/%M/%S"))

qPCR_cDNA_ABA_RealTime <- qPCR_cDNA_ABA %>% 
  mutate(hour = format(qPCR_cDNA_ABA$time,'%H:%M')) %>%
  mutate(time2=as.numeric(time2))

col_ABA <- qPCR_cDNA_ABA_RealTime %>% 
  filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% 
  left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype == "col") %>%
  add_column(condition ="ABA stimulus") %>%
  dplyr::select(-time, -hour, -...5) %>%
  rename(time=time2) %>%
  mutate(time=as.numeric(time)) 

stat_qPCR_cDNA <- qPCR_cDNA_ABA_RealTime %>%
  group_by(gene, time, genotype,hour,time2) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA %>% filter(time2==8) 

#Analysis 
#Short

my_graph <- 
stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(hour %in% c("00:00","00:15","00:30","01:00")) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster short time.pdf", plot=my_graph, width = 10, height = 4)

#Medium
my_graph <- 
  stat_qPCR_cDNA %>% filter(hour %in% c("00:00","02:00","04:00","08:00")) %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) 

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster medium time.pdf", plot=my_graph, width = 7, height = 4)

#Short and medium
my_graph <- 
  stat_qPCR_cDNA %>% filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  )

my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster short&medium time.pdf", plot=my_graph, width = 10, height = 4)


#Col-0 only
my_graph <- 
  stat_qPCR_cDNA %>% filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype == "col") %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph


#Col-0 and SALK ----
my_graph <- 
  stat_qPCR_cDNA %>% filter(hour %in% c("00:00","02:00","04:00","08:00")) %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype %in% c("col","pro")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) +
  facet_wrap(~genotype)

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster SALK.pdf", plot=my_graph, width = 10, height = 4)


#ABA marker genes ----
##Normalization time0 col-0 -----


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_MarkerGenes.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(
         time == 0) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_MarkerGenes_V2.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(
    time == 0) %>%
  clean_genotype_code

position=position_dodge(.10)

stat_qPCR_cDNA <- qPCR_cDNA_ABA %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  qPCR_cDNA_ABA %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster control condition_ABA genes.pdf", plot=my_graph, width = 7, height = 3.5)

##Normalization time0 each genotype ----

#


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_MarkerGenes.xlsx",sheet=1) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_MarkerGenes_V2.xlsx",sheet=1) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code

qPCR_cDNA_ABA$time <- as.POSIXct(strptime(qPCR_cDNA_ABA$time, format="%H/%M/%S"))

qPCR_cDNA_ABA_RealTime <- qPCR_cDNA_ABA %>% 
  mutate(hour = format(qPCR_cDNA_ABA$time,'%H:%M')) %>%
  mutate(time2=as.numeric(time2))

col_ABA <- qPCR_cDNA_ABA_RealTime %>% 
  filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% 
  left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype == "col") %>%
  add_column(condition ="ABA stimulus") %>%
  dplyr::select(-time, -hour, -...5) %>%
  rename(time=time2) %>%
  mutate(time=as.numeric(time)) 

stat_qPCR_cDNA <- qPCR_cDNA_ABA_RealTime %>%
  group_by(gene, time, genotype,hour,time2) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA %>% filter(time2==8) 

#Analysis 
#Short

my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(hour %in% c("00:00","00:15","00:30","01:00")) %>%
  filter(gene %in% c("RD29B","RAB18")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster short time_ABA.pdf", plot=my_graph, width = 10, height = 4)

#Medium
my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("RD29B","RAB18")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster all time_ABA Marker Gene.pdf", plot=my_graph, width = 10, height = 4)


# marneral cluster gene in MRN mutants ----
##Normalization time0 col-0 -----


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V2.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","35smrn1","mrn1"), 
    time == 0) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V3.xlsx",sheet=2) %>%
  mutate(time=as.numeric(time),ddCt=as.numeric(ddCt)) %>%
  filter(genotype %in% c("col","35smrn1","mrn1"), 
         time == 0) %>%
  clean_genotype_code

position=position_dodge(.10)

stat_qPCR_cDNA <- qPCR_cDNA_ABA %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  qPCR_cDNA_ABA %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster control condition_MRN1 mutants marneral genes.pdf", plot=my_graph, width = 7, height = 3.5)

##Normalization time0 each genotype ----

#


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V2.xlsx",sheet=1) %>%
  filter(genotype %in% c("col","35smrn1","mrn1")) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V3.xlsx",sheet=1) %>%
  filter(genotype %in% c("col","35smrn1","mrn1")) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code

qPCR_cDNA_ABA$time <- as.POSIXct(strptime(qPCR_cDNA_ABA$time, format="%H/%M/%S"))

qPCR_cDNA_ABA_RealTime <- qPCR_cDNA_ABA %>% 
  mutate(hour = format(qPCR_cDNA_ABA$time,'%H:%M')) %>%
  mutate(time2=as.numeric(time2))

col_ABA <- qPCR_cDNA_ABA_RealTime %>% 
  filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% 
  left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype == "col") %>%
  add_column(condition ="ABA stimulus") %>%
  dplyr::select(-time, -hour, -...5) %>%
  rename(time=time2) %>%
  mutate(time=as.numeric(time)) 

stat_qPCR_cDNA <- qPCR_cDNA_ABA_RealTime %>%
  group_by(gene, time, genotype,hour,time2) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA %>% filter(time2==8) 

#Analysis 
#Short

my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(hour %in% c("00:00","00:15","00:30","01:00")) %>%
  filter(gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster short time_marneral genes marneral mutants.pdf", plot=my_graph, width = 10, height = 4)

#Medium
my_graph <- 
  stat_qPCR_cDNA %>% filter(hour %in% c("00:00","02:00","04:00","08:00")) %>% left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16", "MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  facet_wrap(~genotype_name, nrow=1) +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster all time_MRN mutants medium marneral genes.pdf", plot=my_graph, width = 10, height = 4)


# MARS expression in mars mutants --------

##Normalization time0 each genotype ----

#


qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V2.xlsx",sheet=3) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code

qPCR_cDNA_ABA <- read_excel("data/qPCR_cDNA/qPCR_cDNA_ABA_V3.xlsx",sheet=3) %>%
  mutate(ddCt=as.numeric(ddCt)) %>%
  clean_genotype_code



qPCR_cDNA_ABA$time <- as.POSIXct(strptime(qPCR_cDNA_ABA$time, format="%H/%M/%S"))

qPCR_cDNA_ABA_RealTime <- qPCR_cDNA_ABA %>% 
  mutate(hour = format(qPCR_cDNA_ABA$time,'%H:%M')) %>%
  mutate(time2=as.numeric(time2))

col_ABA <- qPCR_cDNA_ABA_RealTime %>% 
  filter(hour %in% c("00:00","00:15","00:30","01:00", "02:00","04:00","08:00")) %>% 
  left_join(genotype_metadata) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL"),
         genotype == "col") %>%
  add_column(condition ="ABA stimulus") %>%
  dplyr::select(-time, -hour, -...5) %>%
  rename(time=time2) %>%
  mutate(time=as.numeric(time)) 

stat_qPCR_cDNA <- qPCR_cDNA_ABA_RealTime %>%
  group_by(gene, time, genotype,hour,time2) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA %>% filter(time2==8) 

#Analysis 
#Short

my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(hour %in% c("00:00","00:15","00:30","01:00")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster short time_MARS in mars mutants.pdf", plot=my_graph, width = 10, height = 4)

#Medium
my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>% filter(hour %in% c("00:00","02:00","04:00","08:00")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position = position_dodge(width = 0.5))+
  geom_hline(yintercept=0, linetype=1) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  scale_x_datetime() +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  ) 

my_graph



# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster all time_MARS in mars mutants.pdf", plot=my_graph, width = 10, height = 4)








# Nitrogen starvation ----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=1)


##Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_nitrate <- ddCt_data %>%
  add_column(condition="NO3- starvation")

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time %in% c("0","12","24")) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral nitrogen starvation.pdf", plot=my_graph, width = 5, height = 3)


# Pi starvation -----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=2) %>%
    clean_genotype_code

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_Pi <- ddCt_data %>%
  add_column(condition="Pi starvation") %>%
  filter(time %in% c("0","2","1", "12","24"))

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time %in% c("0","2","1", "12","24")) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral Pi starvation.pdf", plot=my_graph, width = 5, height = 3)


# NAA treatment ----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=3)

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) %>%
  filter(gene != "ref1")
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_naa <- ddCt_data %>%
  add_column(condition="NAA stimulus") %>%
  filter(time != 24)

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>% filter(time != 24) %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral NAA treatment.pdf", plot=my_graph, width = 5, height = 3)

#Heat Stress -----

qPCR_cDNA_OtherStress <- read_excel("data/qPCR_cDNA/qPCR_cDNA_OtherStress.xlsx", sheet=4)

#Data processing
ref_gene <- "ref2"
ref_data <- qPCR_cDNA_OtherStress %>%
  filter(gene == ref_gene) %>%
  rename(Cp_ref=Cp) %>%
  dplyr::select(-gene)
gene_data <- qPCR_cDNA_OtherStress %>%
  filter(gene != ref_gene) 
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)
dCt_data_control <- dCt_data %>%
  filter(time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))
ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt))

col_hs <- ddCt_data %>%
  add_column(condition="Heat stress")

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 


# Graph

my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/marneral Heat Stress treatment.pdf", plot=my_graph, width = 5, height = 3)

# Col-0 all stress combined ----
Col_all_stress <- col_nitrate %>% bind_rows(col_ABA, col_Pi,col_naa, col_hs) %>%
  dplyr::select(time,gene,ddCt,condition)


stat_qPCR_cDNA <- Col_all_stress %>%
  group_by(gene, time,condition) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA$condition <- factor(stat_qPCR_cDNA$condition, c("Pi starvation","NO3- starvation","NAA stimulus", "Heat stress", "ABA stimulus"))


my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) +
  facet_wrap(~condition, nrow=1, scales = "free_x") +
  xlab("time")

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/Col-0 all treatment.pdf", plot=my_graph, width = 10, height = 3)

# Epigentic mutant ----

# LHP1 mutant -------

# Data import
lhp1_qPCR_output <- read_excel("data/qPCR_cDNA/lhp1_qPCR_output.xlsx") %>%
    clean_genotype_code


ref_gene <- "ref2"
ref_data <- lhp1_qPCR_output %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  select(-gene, -Cp)
gene_data <- lhp1_qPCR_output %>%
  filter(gene != ref_gene) 
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 

stat_ddCt_LHP1 <- ddCt_data %>%
  group_by(genotype, gene) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

ddCt_data$gene <- factor(ddCt_data$gene, c("CYP705A12","CYP71A16","MHAL","MRN1"))


my_graph <- 
ddCt_data %>% filter(genotype %in% c("col","lhp1"), gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16")) %>% mutate(ddCt=-ddCt) %>%
  ggbarplot(., x = "gene", y = "ddCt", add = "mean_se",
            fill = "genotype",
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="t.test", label = "p.signif")  +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10)) +
  ylab("log2(FC)") +
  scale_fill_manual(
    labels=genotype_metadata$genotype_name,
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,)


my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/col vs lhp1 t0.pdf", plot=my_graph, width = 5, height = 3)



# CLF mutant----
#normalization time0 Col-0 ----
clf_output_qPCR <- read_excel("data/qPCR_cDNA/clf_output_qPCR.xlsx") %>%
    clean_genotype_code

ref_gene <- "ref2"
ref_data <- clf_output_qPCR %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)
gene_data <- clf_output_qPCR %>%
  filter(gene != ref_gene) 
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)

dCt_data_control <- dCt_data %>%
  filter(genotype == "col", time =="0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 


ddCt_data$gene <- factor(ddCt_data$gene, c("CYP705A12","CYP71A16","MHAL","MRN1"))
ddCt_data$genotype <- factor(ddCt_data$genotype, c("col","clf"))


my_graph <- 
  ddCt_data %>% filter(gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16"), time ==0) %>% 
  mutate(ddCt=-ddCt) %>% left_join(genotype_metadata) %>% 
  ggbarplot(., x = "gene", y = "ddCt", add = "mean_se",
            fill = "genotype_name",
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype_name), method="t.test", label = "p.format")  +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1, size=10))+
  ylab("log2(FC)") +
  scale_fill_manual(
    labels=genotype_metadata$genotype_name,
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,)

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/col vs clf t0.pdf", plot=my_graph, width = 7, height = 4)


#normalization time0 each genotype ----
dCt_data_control <- dCt_data %>%
  filter(genotype == "col", time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_col <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>% 
  filter(genotype =="col")

dCt_data_control <- dCt_data %>%
  filter(genotype == "clf", time =="0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data_clf <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>% 
  filter(genotype =="clf")

ddCt_data <- ddCt_data_col %>% 
  bind_rows(ddCt_data_clf)

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, c("col","clf"))

#Graph

my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype_name) +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/col vs clf kinetic_t0 genotype.pdf", plot=my_graph, width = 5, height = 3)


#normalization time0 col kinetic ----
dCt_data_control <- dCt_data %>%
  filter(genotype == "col", time == "0") %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))

ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) 


stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA$genotype <- factor(stat_qPCR_cDNA$genotype, c("col","clf"))

#Graph

my_graph <- 
  stat_qPCR_cDNA %>%
  filter(gene %in% c("MRN1", "CYP705A12", "CYP71A16","MHAL")) %>% 
  left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=gene)) +
  geom_line(aes(color = gene), size=0.75)+
  geom_point(aes(y=mean), size = 0.75, shape=15)+
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), width = 0.5) +
  geom_hline(yintercept=0, linetype=1) +
  xlab(label = "")+
  ylab(label = "log2(FC)") +
  facet_wrap(~genotype_name) +
  theme(strip.text = element_text(size = rel(1))) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/col vs clf kinetic_t0 col.pdf", plot=my_graph, width = 5, height = 3)



# lhp1 ABA -----

lhp1_output_qPCR <- read_excel("data/qPCR_cDNA/lhp1_output_qPCR.xlsx")

ref_gene <- "ref2"
ref_data <- lhp1_output_qPCR %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- lhp1_output_qPCR %>%
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


ddCt_data_col <- ddCt_data %>% filter(genotype == "col")
ddCt_data_lhp1 <- ddCt_data %>% filter(genotype == "lhp1")

ddCt_data_0 <- ddCt_data_col %>%
  bind_rows(ddCt_data_lhp1)

stat_lhp1 <- ddCt_data_0 %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_lhp1 <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 



stat_lhp1  %>% 	
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=genotype)) +
  geom_line(aes(color = genotype),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)	



my_graph <- 
stat_lhp1 %>%	filter(time %in% c(0,4)) %>% 	
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

my_graph <- 
ddCt_data %>% mutate(ddCt=-ddCt) %>% filter(gene =="MRN1") %>%
  ggbarplot(., x = "time", y = "ddCt", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 6) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=0)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/col vs lhp1 t0t4_t0Norm_MRN1.pdf", plot=my_graph, width = 2, height = 3)


my_graph <- 
  ddCt_data %>% mutate(ddCt=-ddCt) %>% filter(gene =="CYP705A12") %>%
  ggbarplot(., x = "time", y = "ddCt", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 6) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=0)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/col vs lhp1 t0t4_t0Norm_cyp705.pdf", plot=my_graph, width = 2, height = 3)



my_graph <- 
  ddCt_data %>% mutate(ddCt=-ddCt) %>% filter(gene =="CYP71A16") %>%
  ggbarplot(., x = "time", y = "ddCt", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 6) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=0)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/col vs lhp1 t0t4_t0Norm_cyp71A.pdf", plot=my_graph, width = 2, height = 3)


my_graph <- 
  ddCt_data %>% mutate(ddCt=-ddCt) %>% filter(gene =="MHAL") %>%
  ggbarplot(., x = "time", y = "ddCt", add = "mean_se",
            fill = "genotype", 
            position = position_dodge(0.8))+
  stat_compare_means(aes(group = genotype), method="anova", label = "p.format", label.y = 8) +
  theme_bw()+
  theme(axis.text.x  = element_text(angle=0, vjust=0, size=10, hjust=0.5))+
  xlab("")+
  ylab("")+
  geom_hline(yintercept=0)+
  scale_fill_manual(breaks=genotype_metadata$genotype,
                    values=genotype_metadata$genotype_color,
                    labels=genotype_metadata$genotype_name,) + theme(legend.position = "none") 
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/col vs lhp1 t0t4_t0Norm_mhal.pdf", plot=my_graph, width = 2, height = 3)





# High Treatment ABA 100uM -------

qPCR_output_100uM <- read_excel("data/qPCR_cDNA/qPCR_High_ABA_100uM.xlsx",sheet=1) %>%
  mutate(Cp=as.numeric(Cp))

ref_gene <- "ref"
ref_data <- qPCR_output_100uM %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_output_100uM %>%
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

#not normalized


stat_High <- ddCt_data %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


stat_High %>% filter(gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1 )+
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 



#normalized t0

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

ddCt_data_col <- ddCt_data %>% filter(genotype == "col")


dCt_data_control <- dCt_data %>%
  filter(genotype == "18") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

ddCt_data_18 <- ddCt_data %>% filter(genotype == "18")

dCt_data_control <- dCt_data %>%
  filter(genotype == "38") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



ddCt_data_38 <- ddCt_data %>% filter(genotype == "38")


dCt_data_control <- dCt_data %>%
  filter(genotype == "92") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


ddCt_data_92 <- ddCt_data %>% filter(genotype == "92")

ddCt_data_0 <- ddCt_data_col %>%
  bind_rows(ddCt_data_18) %>%
  bind_rows(ddCt_data_38) %>%
  bind_rows(ddCt_data_92) 



stat_High <- ddCt_data_0 %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) %>%
  mutate(genotype = as.factor(genotype))	

stat_High$genotype <- factor(stat_High$genotype, c("col","92","38","18"))

my_graph <- 
stat_High %>% filter(gene %in% c("CYP705A12","CYP71A16","MHAL","MRN1")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1 )+
  scale_color_manual(
    breaks=selected_gene$gene,
    values=selected_gene$gene_color,
    labels=selected_gene$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/100uM ABA marneral cluster.pdf", plot=my_graph, width = 5, height = 2)


my_graph <- 
  stat_High %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1 )+
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph

# Save 
ggsave(filename="out/qPCR_cDNA/100uM ABA resp genes.pdf", plot=my_graph, width = 5, height = 2)


#UBQ:MARS -----

qPCR_output_UBQ <- read_excel("data/qPCR_cDNA/qPCR_UBQ.xlsx") 

ref_gene <- "ref"
ref_data <- qPCR_output_UBQ %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)

gene_data <- qPCR_output_UBQ %>%
  filter(gene != ref_gene) 

# Fusion dataframe (sans le gene de reference) avec Cp du gene de reference pour chaque echantillon (=ref_data)
#& ajout colonne dCp (=variation en nb de cycle/gene par rapport au gene de reference)
dCt_data <- gene_data %>%
  left_join(ref_data) %>%
  mutate(dCp=Cp - Cp_ref)



# Ajout dCp/gene par rapport Ã  la condition controle (ici Col NoABA) = dCt_data_control
dCt_data_control <- dCt_data %>%
  filter(type == "WT") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp, na.rm=TRUE))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))



stat_ubq <- ddCt_data %>%
  group_by(gene, time, type,genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_ubq$type <- factor(stat_ubq$type, c("WT","WT_1","UBQ_1","WT_2","UBQ_2","WT_3","UBQ_3"))

my_graph <- 
stat_ubq %>% left_join(genotype_metadata_ubq) %>%
  ggplot(data=., mapping=aes(x=type, y=mean)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill=genotype)) +
  geom_errorbar(width = .2, position = position_dodge(0.5), mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std))) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  facet_grid(time~gene, scale="free") +
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))  +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata_ubq$genotype,
    values=genotype_metadata_ubq$genotype_color,
    labels=genotype_metadata_ubq$genotype_name,
  )


my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/UBQ_MARS.pdf", plot=my_graph, width = 8, height = 3.5)


stat_ubq %>% left_join(genotype_metadata) %>% filter(!gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=., mapping=aes(x=type, y=mean)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, aes(fill=genotype_name)) +
  geom_errorbar(width = .2, position = position_dodge(0.5), mapping=aes(ymin=(mean - erreur_std), ymax=(mean + erreur_std))) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1, linetype=2) +
  geom_hline(yintercept = -1, linetype=2) +
  facet_grid(time~gene, scale="free") +
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))  +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_color_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )




# RNAi1.8, line 3 kinetic ------

qPCR_cDNA_RNAi_line3 <- read_excel("data/qPCR_cDNA/qPCR_cDNA_RNAi_line3.xlsx")


ref_gene <- "ref2"
ref_data <- qPCR_cDNA_RNAi_line3 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_cDNA_RNAi_line3 %>%
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


ddCt_data_col <- ddCt_data %>% filter(genotype == "col")


dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai18") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))

ddCt_data_rnai<- ddCt_data %>% filter(genotype == "rnai18")

ddCt_data_0 <- ddCt_data_col %>%
  bind_rows(ddCt_data_rnai)

stat_rnai <- ddCt_data_0 %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			


my_graph <- 
stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1","MHAL")) %>% left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1)+
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) 

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/kin RNAi 3 marneral genes.pdf", plot=my_graph, width = 5, height = 2)


my_graph <- 
  stat_rnai %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1")) %>% filter(time <1.5) %>% 
  left_join(genotype_metadata) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype_name,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene2_qPCR$gene,
    values=selected_gene2_qPCR$gene_color,
    labels=selected_gene2_qPCR$gene_name,
  ) 

my_graph
# Save 
ggsave(filename="out/qPCR_cDNA/kin RNAi 3 marneral genes_short.pdf", plot=my_graph, width = 5, height = 2)


#time 0 Line3 ----

qPCR_cDNA_RNAi_line3 <- read_excel("data/qPCR_cDNA/qPCR_cDNA_RNAi_line3.xlsx")


ref_gene <- "ref2"
ref_data <- qPCR_cDNA_RNAi_line3 %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_cDNA_RNAi_line3 %>%
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
  mutate(time=as.numeric(time)) %>%
  filter(time==0)


position=position_dodge(.10)

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/time0 RNAi MARS 3.pdf", plot=my_graph, width = 7, height = 3.5)






# mro kin ------
#kin
qPCR_cDNA_mro_kin <- read_excel("data/qPCR_cDNA/qPCR_cDNA_mro_kin.xlsx")



ref_gene <- "ref"
ref_data <- qPCR_cDNA_mro_kin %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_cDNA_mro_kin %>%
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


ddCt_data_col <- ddCt_data %>% filter(genotype == "col")


dCt_data_control <- dCt_data %>%
  filter(genotype == "mro") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))
ddCt_data_mro <- ddCt_data %>% filter(genotype == "mro")

ddCt_data_0 <- ddCt_data_col %>%
  bind_rows(ddCt_data_mro)

stat_mro <- ddCt_data_0 %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			





stat_mro %>%	
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=genotype)) +
  geom_line(aes(color = genotype),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~gene,	 nrow=1)



my_graph <- 
  stat_mro %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1", "MHAL")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/mro marneral genes.pdf", plot=my_graph, width = 5, height = 2)



my_graph <- 
  stat_mro %>%	filter(gene%in% c("CYP705A12","CYP71A16","MRN1","MHAL")) %>% filter(time <1.5) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) +
  scale_color_manual(
    breaks=selected_gene_qPCR$gene,
    values=selected_gene_qPCR$gene_color,
    labels=selected_gene_qPCR$gene_name,
  ) 

my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/mro marneral genes_short.pdf", plot=my_graph, width = 5, height = 2)

#time 0



dCt_data_control <- dCt_data %>%
  filter(genotype == "col") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time)) %>%
  filter(time=="0")



position=position_dodge(.10)

stat_qPCR_cDNA <- ddCt_data %>%
  group_by(gene, time, genotype) %>%
  summarise(mean=mean(-ddCt), 
            median= median(-ddCt),
            SD=sd(-ddCt), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=SD/sqrt(n)) 

stat_qPCR_cDNA <- 
  ddCt_data %>%
  group_by(gene) %>%
  nest %>%
  mutate(stat = map(data, one_way_anova, "ddCt", "genotype")) %>%
  dplyr::select(-data) %>%
  unnest(stat)

#ANOVa GOOD ------
label_y_shift <-
  max(stat_qPCR_cDNA$mean + stat_qPCR_cDNA$sem) * 0.05



my_graph <- 
  stat_qPCR_cDNA %>% left_join(genotype_metadata) %>%
  ggplot(data=., mapping=aes(x=genotype_name, y=-mean)) +
  geom_bar(stat="identity", position=position, aes(fill=genotype)) +
  geom_errorbar(mapping=aes(ymin=(-mean - sem), ymax=(-mean + sem), width=.5))+
  geom_text(aes(label=letter), nudge_y = -0.6, size=4)+
  theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10))+
  xlab(label = "")+
  ylab(label = "log2(FC)")+
  facet_wrap(~gene, nrow=1) +
  geom_hline(yintercept = 0, linetype=1) +
  theme(strip.text = element_text(size = rel(1)))+
  scale_fill_manual(
    breaks=genotype_metadata$genotype,
    values=genotype_metadata$genotype_color,
    labels=genotype_metadata$genotype_name,
  )
my_graph


# Save 
ggsave(filename="out/qPCR_cDNA/marneral cluster control condition mro mutant_gene.pdf", plot=my_graph, width = 7, height = 3.5)


# water kin ------
#kin

qPCR_output_water <- read_excel("data/qPCR_cDNA/qPCR_water_control.xlsx", sheet=1) 


ref_gene <- "ref"
ref_data <- qPCR_output_water %>%
  filter(gene == ref_gene) %>%
  mutate(Cp_ref=Cp) %>%
  dplyr::select(-gene, -Cp)


gene_data <- qPCR_output_water %>%
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


ddCt_data_col <- ddCt_data %>% filter(genotype == "col")


dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai18") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


ddCt_data_rnai18 <- ddCt_data %>% filter(genotype == "rnai18")



dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai38") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


ddCt_data_rnai38 <- ddCt_data %>% filter(genotype == "rnai38")



dCt_data_control <- dCt_data %>%
  filter(genotype == "rnai92") %>%
  filter(time == 0) %>%
  group_by(gene) %>%
  summarise(dCt_control=mean(dCp))


ddCt_data <- dCt_data %>%
  left_join(dCt_data_control) %>%
  mutate(ddCt = dCp - dCt_control,
         FC = 2^(-ddCt)) %>%
  mutate(time=as.numeric(time))


ddCt_data_rnai92 <- ddCt_data %>% filter(genotype == "rnai92")

ddCt_data_0 <- ddCt_data_col %>%
  bind_rows(ddCt_data_rnai18) %>%
  bind_rows(ddCt_data_rnai38) %>%
  bind_rows(ddCt_data_rnai92)

stat_water <- ddCt_data_0 %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_water %>% filter(gene %in% c("RAB18","RD29B")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 

stat_water %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MARS")) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 


stat_water %>% filter(gene %in% c("CYP705A12","CYP71A16","MRN1", "MARS"), time <2) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +			
  facet_wrap(~genotype,	 nrow=1) 




stat_water %>%	
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=genotype)) +
  geom_line(aes(color = genotype),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +			
  theme_bw()  +
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1)+
  facet_wrap(~gene,	 nrow=1)






























#TEST for norm ------------------


# normalization time 0 col for ABA kin ----------------------

ddCt_data_NotNorm <- read_excel("data/qPCR_cDNA/ddCt_data_NotNorm.xlsx",sheet=4) %>%
  dplyr::select(-exp) %>%
  mutate(ddCt=as.numeric(ddCt),
         time=as.numeric(time))



stat_test <- ddCt_data_NotNorm %>%			
  group_by(gene, time, genotype) %>%	
  summarise(mean=mean(-ddCt),	 		
            median= median(-ddCt),			
            SD=sd(-ddCt),	 #ecart-type		
            n=n(),	 #nombre d'?chantillons		
            erreur_std=SD/sqrt(n)) 			

stat_test %>% filter(gene %in% c("MHAL","MRN1","CYP705A12","CYP71A16"), genotype %in% c("col","pro","rnai38","rnai92"), 
                     time %in% c(0,0.25,0.5,1,2,4,8)) %>%
  ggplot(data=.,	 mapping=aes(x=time,	 y=mean,	 group=gene)) +
  geom_line(aes(color = gene),	 size=0.75) +		
  geom_errorbar(mapping=aes(ymin=mean - erreur_std,	 ymax=mean + erreur_std),	 width=.2)+	
  geom_point(aes(y=mean),	 size = 0.75,	 shape=15)+	
  geom_hline(yintercept=0) +	
  geom_hline(yintercept=1) +	
  geom_hline(yintercept=-1) +	
  theme_bw() +
  facet_wrap(~genotype,	 nrow=1) 




