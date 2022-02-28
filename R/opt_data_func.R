library(Matrix)
#' Optimization Data
#'
#' calculates all needed values to execute quadratic optimization in SVHM
#'
#'
#' @param covariates dataset of covariates of the subjects in a dataset
#' @param training_dataset data frame representing the trainings dataset
#' @param ordered_event_times data frame of all event times ordered in ascending order
#' @param type Type of kernel, either 'gauss' or 'poly' for gaussian or polynomial kernel
#' @param gamma_squared width of gaussian kernel
#' @param d degree of polynomial kernel
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
#'
#'

optimization_data <- function(covariates, training_dataset, ordered_event_times, gamma_squared=.5, d=1) {

  kernel_matrix <- radial_kernel_mat(covariates=covariates, gamma_squared=gamma_squared)

  risk_and_event_matrix <- create_risk_and_event_matrix(training_dataset = training_dataset, ordered_event_times = ordered_event_times)
  risk_vector <- as(t(risk_and_event_matrix$r_mat), 'sparseVector')
  event_matrix <- risk_and_event_matrix$e_mat
  event_vector <- as(t(event_matrix), 'sparseVector')

  rows_in_list <- apply(kernel_matrix, 2, function(x){
                                          Matrix(rep(as(t(event_matrix*x), "sparseVector"),nrow(ordered_event_times)),nrow=length(event_vector), ncol=nrow(ordered_event_times), sparse=TRUE)
                                          }
                        )

  adapted_kernel_matrix <- Matrix(do.call(cbind, rows_in_list),sparse = TRUE)
  adapted_kernel_matrix <- t(t(adapted_kernel_matrix)*event_vector)

  cond_mat <- condition_mat(event_vector, nrow(ordered_event_times))
  weight_vec <- as.vector(t(weight_mat(training_dataset, ordered_event_times)))

  return(list('r_vec'=risk_vector,
              'adap_k_mat'=adapted_kernel_matrix,
              'c_mat'=cond_mat,
              'w_vec'=weight_vec,
              'k_mat'=kernel_matrix,
              'e_vec'=event_vector))
}


#' Time Dependent Optimization Data
#'
#' calculates all needed values to execute quadratic optimization for the time dependent SVHM at the given event time.
#'
#'
#' @param covariates dataset of covariates of the subjects in a dataset
#' @param mat_train matrix of all individuals under risk at the event time
#' @param event_time event time for which data is calculated
#' @param gamma_squared width of gaussian kernel
#'
#' @return {List
#' \code{$adap_k_mat}  matrix on which quadratic optimization will be performed
#' \code{$w_vec}       vector of weights at any event time for all subjects
#' \code{$k_mat}       Gram matrix of covariates
#' \code{$e_vec}       vector indicating vector containing information if a subject experiences an event.
#' }
#'
#'
#' @import Matrix
#'
#'

optimization_time_data <- function(covariates, mat_train, event_time, gamma_squared=.5) {

  event_vector <- rep(-1, nrow(mat_train))
  event_vector[which(mat_train[, 'futime'] == event_time)] <- 1

  kernel_matrix <- radial_kernel_mat(covariates=covariates, gamma_squared=gamma_squared)

  adapted_kernel_matrix <- t(t(event_vector*kernel_matrix)*event_vector)

  weight_vec <- rep(1/nrow(mat_train), nrow(mat_train))
  weight_vec[which(mat_train[, 'futime'] == event_time)] <- 1-1/nrow(mat_train)

  return(list('adap_k_mat'=adapted_kernel_matrix,
              'w_vec'=weight_vec,
              'k_mat'=kernel_matrix,
              'e_vec'=event_vector))
}

