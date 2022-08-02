# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/CPC/", showWarnings = FALSE, recursive = TRUE)

# Coding capacity calculation
Coding_capacity_MHAL <- read_excel("data/CPC/Coding_capacity_MHAL.xlsx") %>%
  mutate(score=as.numeric(score))


Coding_capacity_MHAL$gene <- factor(Coding_capacity_MHAL$gene, levels=c("CYP705A12", "CYP71A16", "MRN1", "AT5G00580.1", "AT5G00580.2", "AT5G00580.3", 
                                                                        "AT5G00580.4","AT5G06335","AT5G06325", "COLDAIR", "APOLO", "ASCO")) # Choisir ordrer du facte_wrap

my_graph <- 
  Coding_capacity_MHAL %>% filter(methods %in% c("CPC1")) %>%
  ggplot(., aes(gene, score)) +
  geom_col(aes(fill=gene), position="dodge") +
  theme (axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(
    breaks=selected_gene_CPC$gene,
    values=selected_gene_CPC$gene_color,
  ) +
  ggtitle("CPC1")


my_graph

# Save 
ggsave(filename="out/CPC/CPC1.pdf", plot=my_graph, width = 5, height = 3)

my_graph <- 
  Coding_capacity_MHAL %>% filter(methods %in% c("CPC2")) %>%
  ggplot(., aes(gene, score)) +
  geom_col(aes(fill=gene), position="dodge") +
  theme (axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10)) +
  geom_hline(yintercept=0.5) +
  scale_fill_manual(
    breaks=selected_gene_CPC$gene,
    values=selected_gene_CPC$gene_color,
  )+
  ggtitle("CPC2")


my_graph

# Save 
ggsave(filename="out/CPC/CPC2.pdf", plot=my_graph, width = 5, height = 3)

 
## ARES ------

Coding_capacity_MHAL <- 
  read_excel("C:/Users/roule/Box/Thesis data add R script and genome wide browser and fulbright ucr if size ok/mhal/data/CPC/Coding_capacity_MHAL.xlsx")%>%
  mutate(score=as.numeric(score))


Coding_capacity_MHAL$gene <- factor(Coding_capacity_MHAL$gene, levels=c("AT5G00580.1","ARES","IAA14.1","IAA14.2","IAA14.3","IAA14.4")) # Choisir ordrer du facte_wrap

my_graph <- 
  Coding_capacity_MHAL %>% filter(methods %in% c("CPC1"),
                                  gene%in%c("AT5G00580.1","ARES","IAA14.1","IAA14.2","IAA14.3","IAA14.4")) %>%
  ggplot(., aes(gene, score)) +
  geom_col(aes(fill=gene), position="dodge") +
  theme (axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10)) +
  geom_hline(yintercept=0) +
  ggtitle("CPC1")+
  theme_bw()+
  theme(legend.position = "none")

my_graph

# Save 
ggsave(filename="CPC1.pdf", plot=my_graph, width = 2.5, height = 2)

my_graph <- 
  Coding_capacity_MHAL %>% filter(methods %in% c("CPC2"),
                                  gene%in%c("AT5G00580.1","ARES","IAA14.1","IAA14.2","IAA14.3","IAA14.4")) %>%
  ggplot(., aes(gene, score)) +
  geom_col(aes(fill=gene), position="dodge") +
  theme (axis.text.x  = element_text(angle=45, vjust=1, hjust=1, size=10)) +
  geom_hline(yintercept=0.5) +
  ggtitle("CPC2")+
  theme_bw()+
  theme(legend.position = "none")


my_graph

# Save 
ggsave(filename="CPC2.pdf", plot=my_graph, width = 2.5, height = 2)







