library(osqp)

#' Optimal solution of SVHM
#'
#' Uses the osqp package to solve the quadratic optimization problem defined by SVHM.
#'
#'
#' @param optimization_data all values needed for optimization in a list with order (risk_vector, adapted_kernel_matrix, cond_mat, weight_vec)
#' @param num_event_times number of event times in the training dataset
#' @param cost cost parameter of the support vector machine of type numeric
#'
#' @return optimal solution for the SVHM
#'
#'
#' @import osqp
#'


opt_sol_osqp <- function(optimizazion_data, num_event_times, cost) {
  risk_vec <- optimizazion_data$r_vec
  adap_kernel_mat <- optimizazion_data$adap_k_mat
  cond_mat <- optimizazion_data$c_mat
  weight_vec <- optimizazion_data$w_vec

  lower_bound <- rep(0, num_event_times + length(risk_vec))
  upper_bound <- c(cost*weight_vec, rep(0, num_event_times))


  # Formulierung des Optimierungsproblems
  invisible(capture.output(model <- osqp(P=adap_kernel_mat, q=-risk_vec, A=cond_mat, l=lower_bound, u=upper_bound, pars=osqpSettings())))
  res_osqp <- model$Solve()
  return(-res_osqp$x)
}

#' Optimal solution of time dependent SVHM
#'
#' Uses the osqp package to solve the quadratic optimization problem defined by the time dependent SVHM.
#'
#'
#' @param e_vec vector indicating if a subject experienced an event at an event time
#' @param k_mat matrix
#' @param w_vec weight vector
#' @param cost cost parameter of the support vector machine of type numeric
#'
#' @return optimal solution for the time dependent SVHM
#'
#'
#' @import osqp
#'


opt_time_sol_osqp <- function(e_vec, k_mat, w_vec, cost) {

  lower_bound <- rep(0, length(e_vec))
  upper_bound <- cost*w_vec


  # Formulierung des Optimierungsproblems
  invisible(capture.output(model <- osqp(P=k_mat, q=-rep(1, length(e_vec)), A=Matrix(diag(length(e_vec)), sparse = TRUE), l=lower_bound, u=upper_bound, pars=osqpSettings())))
  res_osqp <- model$Solve()
  return(-res_osqp$x)
}
