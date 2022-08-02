#' Compute descriptive statistics from a grouped tibble
#' 
#' @param data A data frame.
#' @param measure A string. Contain the column name used to calculate the statistics
#' @return A tibble conating for each grouping the average (mean), the median,
#'   the standard deviation (sd), the number of samples (n) and the stander
#'   deviation of the mean (sem)
#' @examples
#' library(dplyr)
#' root_length <- tibble(
#'   genotype  = rep(c("A", "B"), each=10),
#'   condition = rep(c("T1", "T2"), each=5, 2),
#'   length    = c(rnorm(5, 6, 3), rnorm(5, 7, 3), rnorm(5, 6.5, 3), rnorm(5, 7.4, 3))
#'   )
#' root_length %>% group_by(genotype) %>% compute_desc_stat("length")
#' root_length %>% group_by(condition) %>% compute_desc_stat("length")
#' root_length %>% group_by(genotype, condition) %>% compute_desc_stat("length")
compute_desc_stat <- function(data, measure){
  data %>%
  summarise(
    mean   = mean(!!sym(measure)),
    median = median(!!sym(measure)),
    n      = n(), # sample number
    sd     = sd(!!sym(measure)), # standard deviation
    sem    = sd / sqrt(n) # standard error
  )
}

#' Mark differencial groups with letters according to a One way ANOVA
#' 
#' @param data A data frame.
#' @param measure A string. Contain the column name used to calculate the statistics
#' @param group A string. Contain the column name used as grouping factor
#' @return 
#' @examples
#' library(dplyr)
#' root_length <- tibble(
#'   genotype  = rep(c("A", "B"), each=10),
#'   condition = rep(c("T1", "T2"), each=5, 2),
#'   length    = c(rnorm(5, 6, 3), rnorm(5, 7, 3), rnorm(5, 6.5, 3), rnorm(5, 7.4, 3))
#'   )
#' root_length %>% one_way_anova("length", "genotype")
#' root_length %>% one_way_anova("length", "condition")
#' root_length %>%
#'   group_by(condition) %>%
#'   nest() %>%
#'   mutate(stat=map(data, one_way_anova, "length", "genotype")) %>%
#'   unnest(stat)
one_way_anova <- function(df, measure, group){
  # Extract descriptive stats
  desc_stat <-
    df %>%
    group_by(!!sym(group)) %>%
    compute_desc_stat(measure)

  # Construct the model measure ~ group.
  # Use get to extract the real name of the variables
  aov_model <-
    aov(get(measure) ~ get(group), data=df)

  # Extract the letters according to TukeyHSD for each group
  extract_letters <-
    TukeyHSD(aov_model) %>%
    # extract the result for the grouping variable
    .$`get(group)` %>%
    # extract the p value
    extract_p %>%
    # get the grouping letters
    multcompLetters %>%
    .$Letters %>%
    # format theme for one colmun letter and one for the group
    tibble(letter=., group=names(.))%>%
    setNames(c("letter", group))

  # join the two dataframe and return it
  left_join(desc_stat, extract_letters)
}

