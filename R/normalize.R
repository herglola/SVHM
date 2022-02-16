
#' Normalize
#'
#' normalizes a vector
#'
#' @param df dataframe
#' @param col columns to be normalized
#'
#' @return normalized columns of dataframe
#'
#' @examples {
#' Example with the preloaded mtcars dataset
#' normalize(mtcars, c('disp',  'hp'))
#' }
#'
normalize <- function(df, covariates){
  return(apply(df[covariates], 2, function(x) {(x - min(x)) / (max(x) - min(x))}))
}
