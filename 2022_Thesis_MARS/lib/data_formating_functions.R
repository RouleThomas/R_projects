#' Correct the genotype code to have something homogeneous.
#' 
#' @param data A data frame.
#'   The colunm genotype is expected and will be processed
#' @return A tibble containg with the corrected genotype code in the genotype column
clean_genotype_code <- function(data){
  data %>%
    mutate(
      genotype = tolower(genotype),
      genotype = case_when(
        genotype == "col-0" ~ "col",
        genotype == "col 0" ~ "col",
        genotype == "rnai9.2" ~ "rnai92",
        genotype == "rnai 9.2" ~ "rnai92",
        genotype == "rnai3.8" ~ "rnai38",
        genotype == "rnai 3.8" ~ "rnai38",
        TRUE ~ genotype
      )
    )
}
