library(distances)
#' Gaussian Kernel
#'
#' calculates the Gaussian Kernel value of two inputs
#'
#'
#' @param x first input vector
#' @param y second input vector
#' @param gamma_squared width of the kernel
#'
#' @return gaussian kernel value
#'
#' @examples {
#' x <- runif(n=10)
#' y <- runif(n=10)
#' SVHM:::radial_kernel(x,y,.5)
#' }
#'
radial_kernel <- function(x, y, gamma_squared) {
  return(exp(-(as.numeric(crossprod(x-y)))/(2*gamma_squared)))
}

#' Gaussian Kernel Matrix
#'
#' calculates the gaussian kernel value of the covariates with each other. calculated matrix will be symmetric
#'
#' @param covariates dataset of covariates of the subjects in a dataset
#' @param gamma_squared width of the kernel
#'
#' @return gaussian kernel matrix
#'
#' @import distances
#'
#' @examples {
#' # Example with the preloaded mtcars dataset
#' covariates <- subset( mtcars, select = c('drat', 'wt') )
#' SVHM:::radial_kernel_mat(covariates,.5)
#' }
#'
radial_kernel_mat <- function(covariates, gamma_squared) {
  norm_matrix <- distances(covariates)
  kernel_matrix <- apply(norm_matrix, 1:2, function (x) exp(-(x**2)/(2*gamma_squared)) )
  return(kernel_matrix)
}
