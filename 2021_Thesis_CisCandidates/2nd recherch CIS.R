

X2nd_research_CIS <- read_excel("data/CIS/2nd research CIS.xlsx")
all_info_V201810415 <- read_excel("data/CIS/all_info_V201810415.xlsx") %>% 
  dplyr::select(ID, Araport11_ID, Col_siRNA_detected, Col_DicerCall,downstream_ID,upstream_ID) 

smRNA <- all_info_V201810415 %>%
  inner_join(X2nd_research_CIS) %>%
  filter(Col_siRNA_detected=="TRUE",Col_DicerCall==24)
