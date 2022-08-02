# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/heat stress/", showWarnings = FALSE, recursive = TRUE)


# Data preparation

observation_HS_ABA <- read_excel("./data/germination/observation_HS_ABA.xlsx")

stat_hs <- observation_HS_ABA %>%
  group_by(genotype,exp) %>%
  summarise(mean=mean(percent), 
            (median=median(percent)),
            ecart_type=sd(percent), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n))



#graph 
position=position_dodge(.9)

ggplot(data=stat_hs, mapping=aes(x=genotype, y=mean)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10)) +
  facet_grid(~exp)


my_graph <- 
  observation_HS_ABA %>% 
  ggbarplot(., x = "genotype", y = "percent",add = "mean_se") +
  stat_compare_means(method="wilcox.test", 
                     comparisons = list(c("col", "mro"), c("col", "pro"), c("col", "RNAi.3.8"), c("col", "RNAi.9.2"), c("col","RNAi3") ), 
                     label = "p.label", bracket.size= 0.5) +
  theme_bw() 

my_graph








stat_hs <- observation_HS_ABA %>%
  group_by(genotype) %>%
  summarise(mean=mean(percent), 
            (median=median(percent)),
            ecart_type=sd(percent), #ecart-type
            n=n(), #nombre d'?chantillons
            erreur_std=ecart_type/sqrt(n))



ggplot(data=stat_hs, mapping=aes(x=genotype, y=mean)) +
  geom_bar(stat="identity", position=position) +
  geom_errorbar(mapping=aes(ymin=mean - erreur_std, ymax=mean + erreur_std), position=position, width=.5)+
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))
