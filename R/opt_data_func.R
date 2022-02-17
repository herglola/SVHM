library(Matrix)
library(matchingR)

#' Optimization Data
#'
#' calculates all needed values to execute quadratic optimization in SVHM
#'
#'
#' @param covariates dataset of covariates of the subjects in a dataset
#' @param training_dataset data frame representing the trainings dataset
#' @param ordered_event_times data frame of all event times ordered in ascending order
#' @param gamma_squared width of the kernel
#'
#' @return {List
#' \code{$r_vec}            vector representing at which event times the subjects are under risk
#' \code{$adap_kernel_mat}  matrix on which quadratic optimization will be performed
#' \code{$c_mat}            matrix representing the constraints of the optimization problem
#' \code{$w_vec}            vector of weights at any event time for all subjects
#' \code{$kernel_mat}       Gram matrix of covariates
#' \code{$e_vec}            vector indicating vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m,
#' }
#'
#'
#' @import Matrix
#' @import matchingR
#'
#'

optimization_data <- function(covariates, training_dataset, ordered_event_times, gamma_squared) {

  kernel_matrix <- radial_kernel_mat(covariates=covariates, gamma_squared=gamma_squared)

  risk_and_event_matrix <- create_risk_and_event_matrix(training_dataset = training_dataset, ordered_event_times = ordered_event_times)
  risk_vector <- as(t(risk_and_event_matrix$r_mat), 'sparseVector')
  event_matrix <- risk_and_event_matrix$e_mat
  event_vector <- as(t(event_matrix), 'sparseVector')

  rows_in_list <- apply(kernel_matrix, 2, function(x){
                                          Matrix(rep(as(t(event_matrix*x), "sparseVector"),nrow(ordered_event_times)),nrow=length(event_vector), ncol=nrow(ordered_event_times), sparse=TRUE)
                                          }
                        )

  #rows_in_list <- lapply(rows_in_list, as, "sparseMatrix")
  adapted_kernel_matrix <- Matrix(do.call(cbind, rows_in_list),sparse = TRUE)
  #adapted_kernel_matrix<- Matrix(matchingR:::repcol(adapted_kernel_matrix,nrow(ordered_event_times)), sparse=TRUE)
  adapted_kernel_matrix <- t(t(adapted_kernel_matrix)*event_vector)

  'quad_matrix <-Matrix(0, length(event_vector), length(event_vector), sparse=TRUE)
  n <- nrow(ordered_event_times)
  for (i in 1:nrow(kernel_matrix)) {
    quad_matrix[,((i-1)*n+1):(i*n)] <- as(t(event_matrix*kernel_matrix[,i]), "sparseVector")
  }
  adapted_kernel_matrix <-  t(t(quad_matrix)*event_vector)'

  'adapted_kernel_matrix<-matchingR:::repcol(kernel_matrix,nrow(ordered_event_times))
  adapted_kernel_matrix<-matchingR:::reprow(adapted_kernel_matrix,nrow(ordered_event_times))

  Matrix(adapted_kernel_matrix, sparse = TRUE)
  adapted_kernel_matrix <- Matrix(t(t(event_vector*adapted_kernel_matrix)*event_vector), sparse = TRUE)'


  cond_mat <- condition_mat(event_vector, nrow(ordered_event_times))
  weight_vec <- as.vector(t(weight_mat(training_dataset, ordered_event_times)))

  return(list('r_vec'=risk_vector,
              'adap_k_mat'=adapted_kernel_matrix,
              'c_mat'=cond_mat,
              'w_vec'=weight_vec,
              'k_mat'=kernel_matrix,
              'e_vec'=event_vector))
}
