# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/germination/", showWarnings = FALSE, recursive = TRUE)

# CRISPR perfect -----------------------


crispr_no_ABA <- read_tsv("data/germination/20210913_germination_noABA.tsv", col_type=cols())
crispr_ABA <- read_tsv("data/germination/20210913_germination_ABA.tsv", col_type=cols())

raw_data <-
  bind_rows(crispr_no_ABA, crispr_ABA) %>%
  mutate(total_nb=50,
         max_germination=`112h`) %>%
  pivot_longer(-c(concentration, plate_nb, genotype, total_nb, max_germination),
               values_to="CumGerm",
               names_to="time") %>%
  mutate(time=as.numeric(str_replace(time, "h", "")),
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor,
         radicule_out=CumGerm/total_nb)

genotype_info <-
  tibble(genotype=c("Col-0", "CRIPSR")) %>%
  mutate(genotype=str_replace_all(genotype, "[:. -]", "") %>%
           tolower,
         genotype=factor(genotype, levels=genotype))


# the measurements by themself
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm)) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, time, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  arrange(concentration, genotype, plate_nb, time) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Set the lower bounds of the time intervals: time from before point
    timeBef=lag(time),
    # if time 0, it give NA. et it to 0
    timeBef=if_else(is.na(timeBef), 0, timeBef),
    # correct time 16 to have time before to 1 (we need time range, and we set it from 0-1 for time 1)
    timeBef=if_else(time==16, 1, timeBef),
    # Set the upper bounds of the time intervals: time of current point
    timeAf=time,
    # for time 0 set it to 1h.
    timeAf=if_else(time==0, 1, timeAf)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric
  )

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(time == 112) %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 112,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-time, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)

formated_data %>%
  tail(10)



# 0.5 ABA
## Model fitting


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL3, log = "", legendPos = c(100, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(100, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(100, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(100, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)

## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames('col-crispr') %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames('col-crispr') %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>% left_join(genotype_metadata_germination) %>%
  ggplot(., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
ggsave(filename="out/germination/CRISPR_Speed of germination at T50 0.5uM ABA.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames('col-crispr') %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_germination) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/CRISPR_Time for T50 0.5uM.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA

mod_00_W22 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = W2.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_00_W22, log = "", legendPos = c(100, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_00_W22)

mod_00_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_00_LL3, log = "", legendPos = c(100, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_00_LL3)


## Classical time estimation

anova(mod_00_LL3, mod_00_W22, test = "Chisq")


ED(mod_00_LL3, c(5,50,95), interval="delta")

ED(mod_00_LL3,c(50),display=FALSE)[,2]*100/ED(mod_00_W22,c(50),display=FALSE)[,1]
ED(mod_00_LL3,c(5),display=FALSE)[,2]*100/ED(mod_00_W22,c(5),display=FALSE)[,1]
ED(mod_00_LL3,c(95),display=FALSE)[,2]*100/ED(mod_00_W22,c(95),display=FALSE)[,1]


slope_comp_00 <-
  compParm(mod_00_LL3, "Slope")



slope_comp_letters_00 <- 
  slope_comp_00[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames('col-crispr') %>%
  multcompLetters
slope_comp_letters_00$Letters



slope_stats_00 <-
  coef(summary(mod_00_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_00$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)



# Plot Speed of germination at T50 0uM ABA ----

shift <- max(slope_stats_00$se + slope_stats_00$slope) * .05

my_graph <- 
  slope_stats_00 %>% left_join(genotype_metadata_germination) %>%
  ggplot(data = ., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  ) + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

# Save 
ggsave(filename="out/germination/CRISPR_Speed of germination at T50 0uM ABA.pdf", plot=my_graph, width = 3, height = 5) 


### T50



T50_comp_00 <-
  compParm(mod_00_LL3, "T50")



T50_comp_letters_00 <- 
  T50_comp_00[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames('col-crispr') %>%
  multcompLetters
T50_comp_letters_00




T50_stats_00 <-
  coef(summary(mod_00_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_00$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)





# Plot Time for T50 0uM ABA ----
shift <- max(T50_stats_00$se + T50_stats_00$T50) * .05
my_graph <- 
  T50_stats_00 %>% left_join(genotype_metadata_germination) %>%
  ggplot(data = , aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  ) + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

# Save 
ggsave(filename="out/germination/CRISPR_Time for T50 0uM ABA.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----


stat_count_ABA_radicule_out_12 <-
  raw_data %>% 
  group_by(genotype, time, concentration) %>%
  mutate(radicule_out=radicule_out*100)%>%
  summarise(mean=mean(radicule_out), 
            median=median(radicule_out),
            SD=sd(radicule_out), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12  %>%
  left_join(genotype_metadata_germination) %>%
  filter(concentration == 0) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  xlab(label = "Time after sowing (hours)") +
  ylim(0,100)+
  scale_color_manual(labels=genotype_metadata_germination$genotype_name,
                     breaks=genotype_metadata_germination$genotype_name,
                     values=genotype_metadata_germination$genotype_color,)+ theme(legend.position = "none")

my_graph



# Save 
ggsave(filename="out/germination/CRISPR_Time course germination 0ABA.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12  %>%
  left_join(genotype_metadata_germination) %>%
  filter(concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  xlab(label = "Time after sowing (hours)") +
  ylim(0,100)+
  scale_color_manual(labels=genotype_metadata_germination$genotype_name,
                     breaks=genotype_metadata_germination$genotype_name,
                     values=genotype_metadata_germination$genotype_color,)+ theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/germination/CRISPR_Time course germination 0.5ABA.pdf", plot=my_graph, width = 5, height = 3) 






















































# Data preparation exp1 ------

# time between plant out and midnight
time_shift <- 24 - 10.5 # out at 10h30

raw_data <-
  read_tsv("data/germination/germination.txt", col_type=cols()) %>%
  mutate(CumGerm=total_nb * radicule_out / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor,
         time=case_when(
           day == "0" ~ 0,
           day == "1AM" ~ 10.5 + time_shift,
           day == "1PM" ~ 18.5 + time_shift,
           day == "2AM" ~ 10.5 + 24 + time_shift,
           day == "2PM" ~ 18.5 + 24 + time_shift,
           day == "3AM" ~ 10.5 + 48 + time_shift,
           day == "3PM" ~ 18.5 + 48 + time_shift
         ),
  )

genotype_info <-
  tibble(genotype=c("Col-0", "RNAi 3.8", "RNAi 9.2", "pro",
                         "35S:MRN1", "mrn1")) %>%
  mutate(genotype=str_replace_all(genotype, "[:. -]", "") %>%
           tolower %>%
           str_replace("tdna", ""),
         genotype=factor(genotype, levels=genotype))




# the measurements by themself
formated_measurements <-
  raw_data %>%
  select(-time) %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 10.5 + time_shift,
      day == "2AM" ~ 18.5 + time_shift,
      day == "2PM" ~ 10.5 + 24 + time_shift,
      day == "3AM" ~ 18.5 + 24 + time_shift,
      day == "3PM" ~ 10.5 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 10.5 + time_shift,
      day == "1PM" ~ 18.5 + time_shift,
      day == "2AM" ~ 10.5 + 24 + time_shift,
      day == "2PM" ~ 18.5 + 24 + time_shift,
      day == "3AM" ~ 10.5 + 48 + time_shift,
      day == "3PM" ~ 18.5 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-day, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>% left_join(genotype_metadata_germination) %>%
ggplot(., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
  

my_graph

# Save 
ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
T50_stats_05 %>% left_join(genotype_metadata_germination) %>%
ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA

mod_00_W22 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = W2.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_00_W22, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_00_W22)


## Classical time estimation



ED(mod_00_W22, c(5,50,95), interval="delta")



ED(mod_00_W22,c(50),display=FALSE)[,2]*100/ED(mod_00_W22,c(50),display=FALSE)[,1]
ED(mod_00_W22,c(5),display=FALSE)[,2]*100/ED(mod_00_W22,c(5),display=FALSE)[,1]
ED(mod_00_W22,c(95),display=FALSE)[,2]*100/ED(mod_00_W22,c(95),display=FALSE)[,1]



slope_comp_00 <-
  compParm(mod_00_W22, "Slope")



slope_comp_letters_00 <- 
  slope_comp_00[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_00$Letters



slope_stats_00 <-
  coef(summary(mod_00_W22))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_00$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)



# Plot Speed of germination at T50 0uM ABA ----

shift <- max(slope_stats_00$se + slope_stats_00$slope) * .05

my_graph <- 
slope_stats_00 %>% left_join(genotype_metadata_germination) %>%
ggplot(data = ., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  ) + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

# Save 
ggsave(filename="out/germination/Speed of germination at T50 0uM ABA.pdf", plot=my_graph, width = 3, height = 5) 


### T50



T50_comp_00 <-
  compParm(mod_00_W22, "T50")



T50_comp_letters_00 <- 
  T50_comp_00[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_00




T50_stats_00 <-
  coef(summary(mod_00_W22))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_00$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)





# Plot Time for T50 0uM ABA ----
shift <- max(T50_stats_00$se + T50_stats_00$T50) * .05
my_graph <- 
T50_stats_00 %>% left_join(genotype_metadata_germination) %>%
ggplot(data = , aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_germination$genotype,
    values=genotype_metadata_germination$genotype_color,
    labels=genotype_metadata_germination$genotype_name,
  ) + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

# Save 
ggsave(filename="out/germination/Time for T50 0uM ABA.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----


stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, time, concentration) %>%
  summarise(mean=mean(radicule_out), 
            (median =median(radicule_out)),
            SD=sd(radicule_out), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12  %>% left_join(genotype_metadata_germination) %>%
  filter(concentration == 0) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  scale_color_manual(labels=genotype_metadata_germination$genotype_name,
                    breaks=genotype_metadata_germination$genotype_name,
                    values=genotype_metadata_germination$genotype_color,)+ theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12  %>% left_join(genotype_metadata_germination) %>%
  filter(concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  scale_color_manual(labels=genotype_metadata_germination$genotype_name,
                     breaks=genotype_metadata_germination$genotype_name,
                     values=genotype_metadata_germination$genotype_color,)+ theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA.pdf", plot=my_graph, width = 5, height = 3) 







# Data preparation exp2, plante dès le debut ------

# time between plant out and midnight
time_shift <- 24 - 10.5 # out at 10h30

raw_data <-
  read_excel("data/germination/germination_exp2_corr.xlsx")  %>%
  dplyr::select(-cotyledon_open)%>%
  mutate( total_nb = 50,
          radicule_out= CumGerm/ total_nb * 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor,
         time=case_when(
           day == "0" ~ 0,
           day == "1AM" ~ 10.5 + time_shift,
           day == "1PM" ~ 18.5 + time_shift,
           day == "2AM" ~ 10.5 + 24 + time_shift,
           day == "2PM" ~ 18.5 + 24 + time_shift,
           day == "3AM" ~ 10.5 + 48 + time_shift,
           day == "3PM" ~ 18.5 + 48 + time_shift,
           day == "4AM" ~ 10.5 + 72 + time_shift
         ),
  )%>%
  drop_na()%>%
  dplyr::select(-day)
  



raw_data <-
  read_excel("data/germination/germination_2.xlsx",sheet=3) %>%
  mutate(CumGerm=total_nb * radicule_out / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)


genotype_info <-
  tibble(genotype=c("Col", "rnai18","mro")) %>%
  mutate(genotype=str_replace_all(genotype, "[:. -]", "") %>%
           tolower %>%
           str_replace("tdna", ""),
         genotype=factor(genotype, levels=genotype))



# time between plant out and midnight
time_shift <- 24 - 10.5 # out at 10h30




############ TEST
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm)) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, time, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  arrange(concentration, genotype, plate_nb, time) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Set the lower bounds of the time intervals: time from before point
    timeBef=lag(time),
    # if time 0, it give NA. et it to 0
    timeBef=if_else(is.na(timeBef), 0, timeBef),
    # correct time 24 to have time before to 1 (we need time range, and we set it from 0-1 for time 1)
    timeBef=if_else(time==24, 1, timeBef),
    # Set the upper bounds of the time intervals: time of current point
    timeAf=time,
    # for time 0 set it to 1h.
    timeAf=if_else(time==0, 1, timeAf)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric
  )

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(time == 96) %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 96,
    timeAf = Inf
  )


############TEST 




















# the measurements by themself
formated_measurements <-
  raw_data %>%
  dplyr::select(-time)%>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 10.5 + time_shift,
      day == "2AM" ~ 18.5 + time_shift,
      day == "2PM" ~ 10.5 + 24 + time_shift,
      day == "3AM" ~ 18.5 + 24 + time_shift,
      day == "3PM" ~ 10.5 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 10.5 + time_shift,
      day == "1PM" ~ 18.5 + time_shift,
      day == "2AM" ~ 10.5 + 24 + time_shift,
      day == "2PM" ~ 18.5 + 24 + time_shift,
      day == "3AM" ~ 10.5 + 48 + time_shift,
      day == "3PM" ~ 18.5 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

















# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-time, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM_exp2.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(., aes(x = genotype_name, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
ggsave(filename="out/germination/Speed of germination at T50 0uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


# Save 
ggsave(filename="out/germination/Time for T50 0uM_exp2.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----

raw_data <-
  read_excel("data/germination/germination_2.xlsx",sheet=2) %>%
  mutate(CumGerm=total_nb * radicule_out / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)



stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, time, concentration) %>%
  summarise(mean=mean(CumGerm), 
            (median =median(CumGerm)),
            SD=sd(CumGerm), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12  %>% left_join(genotype_metadata_root) %>%
  filter(concentration == 0) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype_name)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  scale_color_manual(labels=genotype_metadata_root$genotype_name,
                     breaks=genotype_metadata_root$genotype_name,
                     values=genotype_metadata_root$genotype_color,)+ theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA_exp2.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12  %>% left_join(genotype_metadata_root) %>%
  filter(concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=time, y=mean, group=genotype)) +
  geom_line(aes(color = genotype_name), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") +
  scale_color_manual(labels=genotype_metadata_root$genotype_name,
                     breaks=genotype_metadata_root$genotype_name,
                     values=genotype_metadata_root$genotype_color,)+ theme(legend.position = "none")
my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA_exp2.pdf", plot=my_graph, width = 5, height = 3) 



# CRISPR  ------

raw_data <-
  read_excel("data/germination/20210621_germination_ABA.xlsx") %>%
  mutate(CumGerm=germ/total_nb * 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)


genotype_info <-
  tibble(genotype=c("col", "crispr")) %>%
  mutate(genotype=as.factor(genotype))



# time between plant out and midnight
time_shift <- 24 - 10.5 # out at 10h30

# the measurements by themself
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 10.5 + time_shift,
      day == "2AM" ~ 18.5 + time_shift,
      day == "2PM" ~ 10.5 + 24 + time_shift,
      day == "3AM" ~ 18.5 + 24 + time_shift,
      day == "3PM" ~ 10.5 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 10.5 + time_shift,
      day == "1PM" ~ 18.5 + time_shift,
      day == "2AM" ~ 10.5 + 24 + time_shift,
      day == "2PM" ~ 18.5 + 24 + time_shift,
      day == "3AM" ~ 10.5 + 48 + time_shift,
      day == "3PM" ~ 18.5 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-day, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA useless ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

# run this because not statistic so can't run the one before
T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) 


# Plot Time for T50 0.5uM ABA original not statistic so can't apply ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters))+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


# Save 
ggsave(filename="out/germination/Time for T50 0uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----

raw_data <-
  read_excel("data/germination/20210621_germination_ABA.xlsx",sheet=1) %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)



stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, day, concentration) %>%
  summarise(mean=mean(CumGerm), 
            (median =median(CumGerm)),
            SD=sd(CumGerm), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 




# CRISPR respecting time every 12 hrs ------

raw_data <-
  read_excel("data/germination/20210621_germination_ABA.xlsx") %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)


genotype_info <-
  tibble(genotype=c("col", "crispr")) %>%
  mutate(genotype=as.factor(genotype))



# time between plant out and midnight
time_shift <- 24 - 9 # out at 9h

# the measurements by themself
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 9 + time_shift,
      day == "2AM" ~ 21 + time_shift,
      day == "2PM" ~ 9 + 24 + time_shift,
      day == "3AM" ~ 21 + 24 + time_shift,
      day == "3PM" ~ 9 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 9 + time_shift,
      day == "1PM" ~ 21 + time_shift,
      day == "2AM" ~ 9 + 24 + time_shift,
      day == "2PM" ~ 21 + 24 + time_shift,
      day == "3AM" ~ 9 + 48 + time_shift,
      day == "3PM" ~ 21 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-day, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA useless ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

# run this because not statistic so can't run the one before
T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) 


# Plot Time for T50 0.5uM ABA original not statistic so can't apply ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters))+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


# Save 
ggsave(filename="out/germination/Time for T50 0uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----

raw_data <-
  read_excel("data/germination/20210621_germination_ABA.xlsx",sheet=1) %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)



stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, day, concentration) %>%
  summarise(mean=mean(CumGerm), 
            (median =median(CumGerm)),
            SD=sd(CumGerm), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 




## CRISPR second exp --------------


raw_data <-
  read_excel("data/germination/20210913_germination_ABA.xlsx") %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)


genotype_info <-
  tibble(genotype=c("col", "crispr")) %>%
  mutate(genotype=as.factor(genotype))



# time between plant out and midnight
time_shift <- 24 - 10.5 # out at 10h30

# the measurements by themself
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 10.5 + time_shift,
      day == "2AM" ~ 18.5 + time_shift,
      day == "2PM" ~ 10.5 + 24 + time_shift,
      day == "3AM" ~ 18.5 + 24 + time_shift,
      day == "3PM" ~ 10.5 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 10.5 + time_shift,
      day == "1PM" ~ 18.5 + time_shift,
      day == "2AM" ~ 10.5 + 24 + time_shift,
      day == "2PM" ~ 18.5 + 24 + time_shift,
      day == "3AM" ~ 10.5 + 48 + time_shift,
      day == "3PM" ~ 18.5 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-day, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA useless ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

# run this because not statistic so can't run the one before
T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) 


# Plot Time for T50 0.5uM ABA original not statistic so can't apply ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters))+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


# Save 
ggsave(filename="out/germination/Time for T50 0uM_crispr.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----

raw_data <-
  read_excel("data/germination/20210913_germination_ABA.xlsx",sheet=1) %>%
  mutate(CumGerm=germ / total_nb * 100)



stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, day, concentration) %>%
  summarise(mean=mean(CumGerm), 
            (median =median(CumGerm)),
            SD=sd(CumGerm), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA_crispr2.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA_crispr2.pdf", plot=my_graph, width = 5, height = 3) 




# CRISPR respecting time every 12 hrs ------

raw_data <-
  read_excel("data/germination/20210913_germination_ABA.xlsx") %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)


genotype_info <-
  tibble(genotype=c("col", "crispr")) %>%
  mutate(genotype=as.factor(genotype))



# time between plant out and midnight
time_shift <- 24 - 9 # out at 9h

# the measurements by themself
formated_measurements <-
  raw_data %>%
  # filter the time
  filter(!is.na(plate_nb),
         !is.na(CumGerm),
         day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM")) %>%
  # Arrange the data so the different times are on following lines per plate
  # to be able to calculate daily germination
  arrange(concentration, genotype, day, plate_nb) %>%
  group_by(concentration, genotype, plate_nb) %>%
  mutate(
    # Correct cumulative germination
    # more germinated seeds than total number of seeds ;)
    CumGerm=if_else(CumGerm > total_nb, total_nb, CumGerm),
    # more seeds on previous day.
    CumGerm=if_else(CumGerm < lag(CumGerm), lag(CumGerm), CumGerm),
    CumGerm=if_else(is.na(CumGerm), 0, CumGerm),
    # Calculate the daily germination
    DailyGerm=CumGerm - lag(CumGerm),
    # Set the daily germination of time point 0 to 0
    DailyGerm=if_else(is.na(DailyGerm), 0, DailyGerm)) %>%
  ungroup %>%
  mutate(
    # Clean genotype name
    genotype=str_replace_all(genotype, "[:.]", ""),
    # Set the lower bounds of the time intervals
    timeBef=case_when(
      day == "0" ~ 0,
      day == "1AM" ~ 1,
      day == "1PM" ~ 9 + time_shift,
      day == "2AM" ~ 21 + time_shift,
      day == "2PM" ~ 9 + 24 + time_shift,
      day == "3AM" ~ 21 + 24 + time_shift,
      day == "3PM" ~ 9 + 48 + time_shift
    ),
    # Set the upper bounds of the time intervals
    timeAf=case_when(
      day == "0" ~ 1,
      day == "1AM" ~ 9 + time_shift,
      day == "1PM" ~ 21 + time_shift,
      day == "2AM" ~ 9 + 24 + time_shift,
      day == "2PM" ~ 21 + 24 + time_shift,
      day == "3AM" ~ 9 + 48 + time_shift,
      day == "3PM" ~ 21 + 48 + time_shift
    ),
    # compute a unique plate identification
    plate_id=str_c(concentration, genotype, plate_nb) %>% as.factor %>% as.numeric)

# The plant that do not germinate
ungerminated_data <-
  formated_measurements %>%
  # Keep only last time point
  filter(day == "3PM") %>%
  # remove unneeded columns that will be generated
  dplyr::select(-(DailyGerm:timeAf)) %>%
  unique %>%
  mutate(
    # The daily germination seeds (the seeds that did not germinated) is
    # the total number of seeds minus the seeds that already germinated
    DailyGerm=total_nb - CumGerm,
    # The cumulative number of seed is all the seeds
    CumGerm=total_nb,
    # Set the bound of the time interval
    timeBef = 18.5 + 48 + time_shift,
    timeAf = Inf
  )

# fusion of the data
formated_data <-
  bind_rows(formated_measurements, ungerminated_data) %>%
  arrange(plate_id, timeBef) %>%
  group_by(plate_id) %>%
  # Generate a line id that should look like:
  # 1    for first plate, zero time point
  # 1.1  for first plate, time point 1
  # 1.2  ...
  # ...
  # 2    for second plate, zero time point
  # 2.1  for second plate, time point 1
  # 1.2  ...
  # ...
  mutate(rownames=str_c(plate_id, ".", row_number() - 1) %>%
           str_replace("\\.0", "")
  ) %>%
  # Remove unneeded columns
  dplyr::select(-day, -plate_nb) %>%
  # transform it as data.frame since `dcr` package do not work with tibble
  as.data.frame %>%
  # set the row names with the constructed ID
  column_to_rownames("rownames")

formated_data %>%
  head(10)



# 0.5 ABA
## Model fitting



mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0.5),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0.5))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA useless ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters)) +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0.5uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

# run this because not statistic so can't run the one before
T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) 


# Plot Time for T50 0.5uM ABA original not statistic so can't apply ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)")+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph
# Save 
ggsave(filename="out/germination/Time for T50 0.5uM_crispr2.pdf", plot=my_graph, width = 3, height = 5) 


# 0 ABA


mod_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0))

plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0 ABA")
summary(mod_05_LL3)

mod_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           curveid = genotype,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
# Add this model to the plot in red
plot(mod_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(2), add=TRUE)
summary(mod_05_LL2)

anova(mod_05_LL3, mod_05_LL2, test = "Chisq")



## Global genotype effect


modNull_05_LL3 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.3(names=c("Slope", "Max", "T50")),
           type = "event", subset=c(concentration == 0),
           upperl=c(NA,1,NA))
plot(mod_05_LL3, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL3, log = "", type = "none", legendPos = c(80, .5), col=c(3), add=TRUE)

anova(mod_05_LL3, modNull_05_LL3, test = "Chisq")


modNull_05_LL2 <- 
  drc::drm(DailyGerm ~ timeBef + timeAf,
           data = formated_data,
           fct = LL.2(names=c("Slope", "T50")),
           type = "event", subset=c(concentration == 0))
plot(mod_05_LL2, log = "", legendPos = c(80, .5),
     ylab="Proportion of germinated seed", xlab="Time (h)", main="0.5 ABA")
plot(modNull_05_LL2, log = "", type = "none", legendPos = c(80, .5), col=c(4), add=TRUE)

anova(mod_05_LL2, modNull_05_LL2, test = "Chisq")



## Model comparison

AIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)
BIC(mod_05_LL3, mod_05_LL2, modNull_05_LL3, modNull_05_LL2)


## Classical time estimation


ED(mod_05_LL3, c(5,50,95), interval="delta")



ED(mod_05_LL2,c(50),display=FALSE)[,2]*100/ED(mod_05_LL2,c(50),display=FALSE)[,1]
ED(mod_05_LL2,c(5),display=FALSE)[,2]*100/ED(mod_05_LL2,c(5),display=FALSE)[,1]
ED(mod_05_LL2,c(95),display=FALSE)[,2]*100/ED(mod_05_LL2,c(95),display=FALSE)[,1]



## Comparison of the parameters

### Maximum germination rate



compParm(mod_05_LL3, "Max")



### Slope


slope_comp_05 <-
  compParm(mod_05_LL3, "Slope")


slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters

slope_comp_letters_05 <- 
  slope_comp_05[,"p-value"] %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
slope_comp_letters_05$Letters


slope_stats_05 <-
  coef(summary(mod_05_LL2))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "Slope") %>%
  dplyr::select(genotype, slope=Estimate, se=`Std. Error`) %>%
  left_join(slope_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  mutate(slope=-slope) %>%
  left_join(genotype_info)


# Plot Speed of germination at T50 0.5uM ABA ----
shift <- max(slope_stats_05$se + slope_stats_05$slope) * .05
my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  geom_text(aes(y=slope + se + shift, label=letters))+
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

my_graph <- 
  slope_stats_05 %>%
  ggplot(., aes(x = genotype, y = slope, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = slope - se, ymax = slope + se), width=.2) +
  xlab("") + ylab("Speed of germination at T50 (proportion/hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))


my_graph

# Save 
#ggsave(filename="out/germination/Speed of germination at T50 0uM ABA_exp2.pdf", plot=my_graph, width = 3, height = 5) 



T50_comp_05 <-
  compParm(mod_05_LL3, "T50")



T50_comp_letters_05 <- 
  T50_comp_05[,"p-value"] %>%
  p.adjust(method="fdr") %>%
  setNames(str_replace(names(.), "/", "-")) %>%
  multcompLetters
T50_comp_letters_05

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`) %>%
  left_join(T50_comp_letters_05$Letters %>%
              as.data.frame %>%
              setNames("letters") %>%
              rownames_to_column("genotype")) %>%
  left_join(genotype_info)

T50_stats_05 <-
  coef(summary(mod_05_LL3))[,1:2] %>%
  as.data.frame %>%
  rownames_to_column("parameters") %>%
  separate(parameters, c("parameter", "genotype"), sep=':') %>%
  filter(parameter == "T50") %>%
  dplyr::select(genotype, T50=Estimate, se=`Std. Error`)


# Plot Time for T50 0.5uM ABA ----

shift <- max(T50_stats_05$se + T50_stats_05$T50) * .05

my_graph <- 
  T50_stats_05 %>% left_join(genotype_metadata_root) %>%
  ggplot(data = ., aes(x = genotype_name, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  geom_text(aes(y=T50 + se + shift, label=letters)) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  scale_fill_manual(
    breaks=genotype_metadata_root$genotype,
    values=genotype_metadata_root$genotype_color,
    labels=genotype_metadata_root$genotype_name,
  )+
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph

my_graph <- 
  T50_stats_05 %>% 
  ggplot(data = ., aes(x = genotype, y = T50, fill=genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = T50 - se, ymax = T50 + se), width=.2) +
  xlab("") + ylab("Time to 50% of germination (hour)") +
  xlab("") + ylab("Time to 50% of germination (hour)") + theme(legend.position = "none")+
  theme (axis.text.x  = element_text(angle=90, vjust=.5, size=10))
my_graph


# Save 
ggsave(filename="out/germination/Time for T50 0uM_crispr2.pdf", plot=my_graph, width = 3, height = 5) 




# Time course of seed germination -----

raw_data <-
  read_excel("data/germination/20210621_germination_ABA.xlsx",sheet=1) %>%
  mutate(CumGerm=total_nb * germ / 100,
         concentration=factor(concentration),
         plate_nb=factor(plate_nb),
         genotype=genotype %>% tolower %>% factor)



stat_count_ABA_radicule_out_12 <- raw_data %>% 
  group_by(genotype, day, concentration) %>%
  summarise(mean=mean(CumGerm), 
            (median =median(CumGerm)),
            SD=sd(CumGerm), #ecart-type
            n=n(), #nombre d'?chantillons
            SE=SD/sqrt(n)) 

my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph




# Save 
ggsave(filename="out/germination/Time course germination 0ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 


my_graph <- 
  stat_count_ABA_radicule_out_12 %>%
  filter(day %in% c("0", "1AM", "1PM", "2AM", "2PM", "3AM", "3PM"), concentration == 0.5) %>%
  ggplot(data=., mapping=aes(x=day, y=mean, group=genotype)) +
  geom_line(aes(color = genotype), size=0.75) +
  geom_errorbar(mapping=aes(ymin=mean - SE, ymax=mean + SE), width=.1)+
  geom_point(aes(y=mean), size = 1, shape=15) +
  theme_bw() +
  theme (axis.text.x  = element_text(angle=0, vjust=.5, size=10))+
  ylab(label = "percent") + theme(legend.position = "none")

my_graph


# Save 
ggsave(filename="out/germination/Time course germination 0.5ABA_crispr.pdf", plot=my_graph, width = 5, height = 3) 







