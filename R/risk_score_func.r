
#' Training risk scores
#'
#' calculates the risk scores for all individuals in the training dataset.
#'
#'
#' @param gamma_sol optimal solution of the SVHM
#' @param kernel_mat Gram matrix of the covariates
#' @param num_event_times number of event times that occour in the training data set
#'
#' @return vector of risk scores for all training subjects
#'
#' @note The calculated risk scores are not the actual risk scores defined by the Risk function but the induce an ordering of the risk scores. For detailed information see reference
#'
#'
#' @references Wang, Y., Chen, T., and Zeng, D. Support vector hazards machine: A counting process framework for learning risk scores for censored outcomes. Journal of Machine Learning Research, 17(167):1-37, 2016
#'
#'

risk_score_training <- function(gamma_sol, kernel_mat, event_vec, num_event_times, training_set_size) {
  risk_training <- rep(Inf, nrow(kernel_mat))
  v <- gamma_sol*event_vec
  m <- matrix(v, nrow = num_event_times)
  jsum_vec <- colSums(m)
  for (i in 1:length(risk_training)) {
    # kernel_mat is symmetric
    risk_training[i] <- sum(jsum_vec*kernel_mat[,i])
  }
  return(risk_training)
}

#' risk scores
#'
#' calculates the risk scores for one individual with the help of the calculated optimal solution to the quadratic programming problem of SVHM and the kernel matrix of the covariates of the test dataset.
#'
#'
#' @param gamma_sol optimal solution of the SVHM
#' @param event_vec vector containing information of the training if a subject in the training dataset is at risk or if an event happens.
#'                  If n are the number of subjects in the training dataset and m the number of event times in the training dataset, then event_vec has length n*m
#' @param v covariates of the individual for which the risk is to be calculated
#' @param covariates_train dataset of covariates of the subjects in the training dataset
#' @param num_event_times number of event times that occour in the training data set
#' @param gamma_squared width of the kernel
#'
#' @return risk score of the individual
#'
#' @note The calculated risk score is not the actual risk scores defined by the Risk function but it induce an ordering of the risk scores. For detailed information see reference
#'
#' @references Wang, Y., Chen, T., and Zeng, D. Support vector hazards machine: A counting process framework for learning risk scores for censored outcomes. Journal of Machine Learning Research, 17(167):1-37, 2016
#'
#'
risk_score <- function(gamma_sol, event_vec, v, covariates_train, num_event_times, gamma_squared){
  kernel_vec <- sapply(1:nrow(covariates_train), function(x) {
                      radial_kernel(v, as.numeric(covariates_train[x,]),gamma_squared)
                      })
  kernel_vec <- rep(as.matrix(kernel_vec), each=num_event_times)
  return(sum(gamma_sol*event_vec*as.vector(kernel_vec)))
}
