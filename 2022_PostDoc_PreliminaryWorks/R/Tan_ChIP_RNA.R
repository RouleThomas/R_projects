# Coverage plot to compare peak size EMF2 Vs H3K27me3

## Import/Tidy Greenscreen bed peaks
EMF2_Greenscreen <- read_delim("2022_PostDoc_PreliminaryWorks/R/Tan_ChIP_RNA/ChIP/EMF2_pool_peaks.narrowPeak.bed", 
                                         delim = "\t", escape_double = FALSE, 
                                         col_names = FALSE, trim_ws = TRUE) %>%
  mutate(length=X3-X2) %>%
  add_column(peak="EMF2",analyses="greenscreen") %>%
  dplyr::select(peak, length,analyses) 

H3K27me3_Greenscreen <- read_delim("2022_PostDoc_PreliminaryWorks/R/Tan_ChIP_RNA/ChIP/H3K27me3_pool_peaks.broadPeak.bed", 
                                         delim = "\t", escape_double = FALSE, 
                                         col_names = FALSE, trim_ws = TRUE)%>%
  mutate(length=X3-X2) %>%
  add_column(peak="H3K27me3",analyses="greenscreen") %>%
  dplyr::select(peak, length,analyses)


### Plot

dummy <- EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  group_by(peak) %>%
  summarize(median = median(length))

EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  ggplot(., aes(x=length, fill=peak)) +
  geom_density(alpha=0.4) +
  xlim(0,10000) +
  geom_vline(data=dummy,aes(xintercept = median, color = peak))


EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  ggplot(., aes(x=peak, y=length, fill=peak)) +
  geom_boxplot()+
  ylim(0,10000)


## Import/Tidy Tan bed peaks

EMF2_Tan <- read_delim("2022_PostDoc_PreliminaryWorks/R/Tan_ChIP_RNA/ChIP/EMF2_peaks_Tan.bed", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE) %>%
  mutate(length=X3-X2) %>%
  add_column(peak="EMF2",analyses="Tan") %>%
  dplyr::select(peak, length,analyses)


H3K27me3_peaks_Tan <- read_delim("2022_PostDoc_PreliminaryWorks/R/Tan_ChIP_RNA/ChIP/H3K27me3_peaks_Tan.bed", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE) %>%
  mutate(length=X3-X2) %>%
  add_column(peak="H3K27me3",analyses="Tan") %>%
  dplyr::select(peak, length,analyses)



### Plot

dummy <- EMF2_Tan %>% bind_rows(H3K27me3_peaks_Tan) %>%
  group_by(peak) %>%
  summarize(median = median(length))

EMF2_Tan %>% bind_rows(H3K27me3_peaks_Tan) %>%
  ggplot(., aes(x=length, fill=peak)) +
  geom_density(alpha=0.4) +
  xlim(0,10000) +
  geom_vline(data=dummy,aes(xintercept = median, color = peak))



### Plot both


## median
dummy <- EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  bind_rows(EMF2_Tan) %>% bind_rows(H3K27me3_peaks_Tan) %>%
  group_by(peak, analyses) %>%
  summarize(median = median(length))

EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  bind_rows(EMF2_Tan) %>% bind_rows(H3K27me3_peaks_Tan) %>%
  ggplot(., aes(x=peak, y=length, fill=peak)) +
  geom_boxplot()+
  facet_wrap(~analyses)+
  ylim(0,10000)+ 
  geom_text(data = dummy, aes(x = peak, y = median, label = median), 
            size = 5, vjust = -0.5)

## mean
dummy <- EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  bind_rows(EMF2_Tan) %>% bind_rows(H3K27me3_peaks_Tan) %>%
  group_by(peak, analyses) %>%
  summarize(mean = mean(length))

EMF2_Greenscreen %>% bind_rows(H3K27me3_Greenscreen) %>%
  bind_rows(EMF2_Tan) %>% bind_rows(H3K27me3_peaks_Tan) %>%
  ggplot(., aes(x=peak, y=length, fill=peak)) +
  geom_boxplot()+
  facet_wrap(~analyses)+
  ylim(0,10000)+ 
  geom_text(data = dummy, aes(x = peak, y = mean, label = mean), 
            size = 5, vjust = -0.5)



## 








