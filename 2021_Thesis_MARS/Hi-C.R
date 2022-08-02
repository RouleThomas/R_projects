# Theme loading
# chrdir = TRUE to allow the included source to work
source("lib/preambule.R", chdir = TRUE)

# Out directory if it does not exist
dir.create("out/Hi-C/", showWarnings = FALSE, recursive = TRUE)


# PNAS root Shoot -----
marneral_interaction_shoot_root <- read_excel("data/Hi-C/marneral_interaction_shoot_root.xlsx")


ggplot(data=marneral_interaction_shoot_root, aes(x=start, y=end, fill=organ)) + 
  geom_raster() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlim(17023080, 17059473) + 
  ylim(17023080, 17059473) + 
  coord_fixed()

ggplot(data=marneral_interaction_shoot_root, aes(x=start, y=end, fill=organ)) + 
  geom_raster() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed()



setwd("~/mhal")




setwd("C:/Bureau/Covid19 Work/Copie Covid/Analyse Hi-C/R")

sig_HiC <- read_excel("data/Hi-C/Hi-C - modified.XLSX", sheet=1)
sig_HiC_lhp1 <- read_excel("data/Hi-C/Hi-C - modified.XLSX", sheet=2)


marneral <-
  sig_HiC %>%
  dplyr::select(chr1=`Chr(1)`, start_1=`start(1)`,
         chr2=`chr(2)`, start_2=`start(2)`,
         z_score=`Z-score`,
         LogP,
         Circos_Thickness=`Circos Thickness`,
         Distance) %>%
  filter(chr1=="Chr5",
         chr2=="Chr5",
         start_1 > 16994451,
         start_2 > 16994451,
         start_1 < 17087850,
         start_2 < 17087850
  )


ggplot(data=marneral, aes(x=start_1, y=start_2, fill=z_score)) + geom_raster() + 
  geom_abline(intercept = 0, slope = 1) + xlim(17023080, 17059473) + ylim(17023080, 17059473) + coord_fixed()


#TFL1 

5:1015493-1048007

TFL1 <-
  sig_HiC_lhp1 %>%
  dplyr::select(chr1=`Chr(1)`, start_1=`start(1)`,
                chr2=`chr(2)`, start_2=`start(2)`,
                z_score=`Z-score`,
                LogP,
                Circos_Thickness=`Circos Thickness`,
                Distance) %>%
  filter(chr1=="Chr5",
         chr2=="Chr5",
         start_1 > 1015493,
         start_2 > 1015493,
         start_1 < 1048007,
         start_2 < 1048007
  )


ggplot(data=TFL1, aes(x=start_1, y=start_2, fill=z_score)) + geom_raster() + 
  geom_abline(intercept = 0, slope = 1) + xlim(1015493, 1048007) + ylim(1015493, 1048007) + coord_fixed()






# Liu

table_S2_chromatinLoop <- read_delim("C:/Users/roule/OneDrive/Bureau/Covid19 Work/Paper_mhal/qPCR_tables/Figures/Hi-C/table_S2_chromatinLoop.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)



MHAL <-
  table_S2_chromatinLoop %>%
  dplyr::select(start_1=`Region_a_to`,
         start_2=`Region_b_to`,
         p_value,
         q_value,
         Distance, 
         Chromosome) %>%
  filter(Chromosome == 5,
         start_1 > 16994451,
         start_2 > 16994451,
         start_1 < 17087850,
         start_2 < 17087850
  )


ggplot(data=MHAL, aes(x=start_1, y=start_2, fill=p_value)) + geom_raster() + geom_abline(intercept = 0, slope = 1) + xlim(17023080, 17059473) + ylim(17023080, 17059473) + coord_fixed()

write.csv(MHAL, file="MHAL")




# LATERALINC ----





marneral <-
  sig_HiC_lhp1 %>%
  dplyr::select(chr1=`Chr(1)`, start_1=`start(1)`,
                chr2=`chr(2)`, start_2=`start(2)`,
                z_score=`Z-score`,
                LogP,
                Circos_Thickness=`Circos Thickness`,
                Distance) %>%
  filter(chr1=="Chr4",
         chr2=="Chr4",
         start_1 > 8325768,
         start_2 > 8325768,
         start_1 < 8370570,
         start_2 < 8370570
  )


ggplot(data=marneral, aes(x=start_1, y=start_2, fill=z_score)) + geom_raster() + 
  geom_abline(intercept = 0, slope = 1) + xlim(8325768, 8370570) + ylim(8325768, 8370570) + coord_fixed()



table_S2_chromatinLoop <- read_delim("data/Hi-C/table_S2_chromatinLoop.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)



MHAL <-
  table_S2_chromatinLoop %>%
  select(start_1=`Region_a_from`,
         start_2=`Region_b_from`,
         p_value,
         q_value,
         Distance, 
         Chromosome) %>%
  filter(Chromosome == 4,
         start_1 > 8325768,
         start_2 > 8325768,
         start_1 < 8370570,
         start_2 < 8370570
  )


ggplot(data=MHAL, aes(x=start_1, y=start_2, fill=p_value)) + geom_raster() + geom_abline(intercept = 0, slope = 1) + 
  xlim(8325768, 8370570) + ylim(8325768, 8370570) + coord_fixed()


# TFL1 -----

5:1015493-1048007

TFL1 <-
  table_S2_chromatinLoop %>%
  dplyr::select(start_1=`Region_a_to`,
                start_2=`Region_b_to`,
                p_value,
                q_value,
                Distance, 
                Chromosome) %>%
  filter(Chromosome == 5,
         start_1 > 1015493,
         start_2 > 1015493,
         start_1 < 1048007,
         start_2 < 1048007
  )


ggplot(data=TFL1, aes(x=start_1, y=start_2, fill=p_value)) + geom_raster() + geom_abline(intercept = 0, slope = 1) + xlim(1015493, 1048007) + ylim(1015493, 1048007) + coord_fixed()

write.csv(TFL1, file="TFL1")


