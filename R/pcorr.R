
library(tidyverse)
library(DT)

#' @title Partial correlation test with covariates and x y variables
#'
#' @export
#' @param anth string vector variable
#' @param clean a dataframe variable
#' @param data a dataframe variable
#' @param covs a string vector variable for all covariates
#'
#' @return list of two dataframes (partial correlation table and a p-value count table)

get_pcor <- function(anth = anthro, clean = clean_dat, data = dat, covs){
  m_tab <- data.frame(matrix(, nrow=length(clean), ncol=0))
  anthro_n <- c()
  for (ant in anth){
    m0 <- data[,c(covs, ant)] %>% na.omit()
    anthro_n <- c(anthro_n, length(m0[,1]))
    m0r <- c();m0p <- c()
    for (met in colnames(clean)){
      m1 <- data[,c(covs, ant, met)] %>% na.omit()
      x_data <- m1[, c(covs, ant)]; y_data <- m1[, c(covs, met)]
      x_f <- as.formula(paste(ant, paste(covs, collapse = " +"), sep = " ~ "))
      y_f <- as.formula(paste(met, paste(covs, collapse = " +"), sep = " ~ "))
      x_fit <- lm(x_f, data = x_data); y_fit <- lm(y_f, data = y_data)
      ex <- resid(x_fit); ey <- resid(y_fit)
      xy_cor <- cor.test(ex, ey); rvalue <- xy_cor$estimate; pvalue <- xy_cor$p.value
      m0r <- c(m0r, rvalue); m0p <- c(m0p, pvalue)
    }
    m_tab <- cbind(m_tab, m0r, m0p)
  }
  colnames(m_tab) <- sapply(anthro, function(x){c(paste0("corr_", x), paste0("p_", x))}) %>% c()
  rownames(m_tab) <- colnames(clean)
  m_tab %>% dplyr::select(starts_with("p_")) -> m_p
  c(apply(m_p, 2, FUN = function(x) {sum(x < 0.05/(1013*5))}),
    sum(apply(m_p, 1, FUN = function(x) {any(x < (0.05/(1013*5)))})),
    sum(apply(m_p, 1, FUN = function(x) {all(x < (0.05/(1013*5)))}))) -> m_pcount
  m_pcount <- cbind(m_pcount, c(anthro_n, " ", " ")) %>% as.data.frame()
  m_pcount %>% magrittr::set_colnames(c("# of metabolite P < 0.05/(1013*5)", "Sample_size")) %>%
    magrittr::set_rownames(c(anth, "at least one anthro_var", "all anthro_vars"))-> m_pcount
  return(list(m_tab, m_pcount))
}

#' @title Display a table in rmarkdown file only
#' @export
#' @param df a dataframe (table) variable
#' @param title a string variable
#' @param row_name a boolean variable
#'
#' @return nothing returns, just to display a table in rmarkdown

table_display <- function(df, title, row_name = TRUE){
  df %>% datatable(rownames = row_name, caption = title,filter="top",
                   extensions = 'Buttons',options = list(pageLength = 15,scrollX=T,
                                                         autoWidth = TRUE,
                                                         dom = 'Blfrtip',
                                                         buttons = c('copy', 'csv', 'excel','pdf', 'print')))
}


