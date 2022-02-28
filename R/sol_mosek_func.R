library(Matrix)
library(Rmosek)

#' Optimal solution of SVHM
#'
#' Uses the Rmosek package to solve the quadratic optimization problem defined by SVHM.
#'
#'
#' @param optimization_data all values needed for optimization in a list with order (risk_vector, adapted_kernel_matrix, cond_mat, weight_vec)
#' @param num_event_times number of event times in the training dataset
#' @param cost cost parameter of the support vector machine of type numeric
#'
#' @return optimal solution for the SVHM
#'
#' @note Rmosek requires a license to use!
#'
#' @import Matrix
#' @import Rmosek
#'


opt_sol_mosek <- function(optimizazion_data, num_event_times, cost) {
  risk_vec <- optimizazion_data$r_vec

  adap_kernel_mat <- optimizazion_data$adap_k_mat
  stopifnot(isSymmetric(adap_kernel_mat) == TRUE)
  adap_kernel_mat <- as(adap_kernel_mat, "dgTMatrix")


  adap_kernel_mat_as_triplet <- data.frame(
    i = adap_kernel_mat@i + 1,  # m@i fängt bei 0 an mit der Nummerrierung
    j = adap_kernel_mat@j + 1,  # m@j fängt bei 0 an mit der Nummerrierung
    x = adap_kernel_mat@x
  )

  # Mosek erwartet untere Dreiecksmatrix zur Optimierung (Matrix wird als symmetrisch angenommen)
  adap_kernel_mat_as_triplet <- adap_kernel_mat_as_triplet[!(adap_kernel_mat_as_triplet$i<adap_kernel_mat_as_triplet$j),]

  cond_mat <- optimizazion_data$c_mat
  weight_vec <- optimizazion_data$w_vec

  prob <- list(sense="max")
  prob
  prob$c <- as.vector(risk_vec)
  prob$A <- Matrix(tail(cond_mat, n=num_event_times), sparse = TRUE)
  prob$bc <- rbind(blc=rep(0, num_event_times),
                   buc=rep(0, num_event_times))
  prob$bx <- rbind(blx=rep(0, length(risk_vec)),
                   bux=c(cost*weight_vec))

  prob$qobj$i <- adap_kernel_mat_as_triplet$i
  prob$qobj$j <- adap_kernel_mat_as_triplet$j
  prob$qobj$v <- -adap_kernel_mat_as_triplet$x
  prob

  capture.output(res_mosek <- mosek(prob))

  # Return the solution
  stopifnot(identical(res_mosek$response$code, 0))
  return(res_mosek$sol$itr$xx)
}


#' Optimal solution of SVHM
#'
#' Uses the Rmosek package to solve the quadratic optimization problem defined by SVHM.
#'
#'
#' @param optimization_data all values needed for optimization in a list with order (risk_vector, adapted_kernel_matrix, cond_mat, weight_vec)
#' @param num_event_times number of event times in the training dataset
#' @param cost cost parameter of the support vector machine of type numeric
#'
#' @return optimal solution for the SVHM
#'
#'
#' @import Matrix
#' @import Rmosek
#'


opt_time_sol_mosek <- function(e_vec, k_mat, w_vec, cost) {
  stopifnot(isSymmetric(k_mat) == TRUE)
  k_mat <- as(k_mat, "dgTMatrix")


  k_mat_as_triplet <- data.frame(
    i = k_mat@i + 1,  # m@i fängt bei 0 an mit der Nummerrierung
    j = k_mat@j + 1,  # m@j fängt bei 0 an mit der Nummerrierung
    x = k_mat@x
  )

  # Mosek erwartet untere Dreiecksmatrix zur Optimierung (Matrix wird als symmetrisch angenommen)
  k_mat_as_triplet <- k_mat_as_triplet[!(k_mat_as_triplet$i<k_mat_as_triplet$j),]

  prob <- list(sense="max")
  prob$A <- Matrix(rep(0, length(e_vec)), nrow = 1, ncol =length(e_vec) , sparse = TRUE)
  prob$c <- as.vector(rep(1, length(e_vec)))
  prob$bc <- rbind(blc=-Inf,
                   buc=Inf)
  prob$bx <- rbind(blx=rep(0, length(e_vec)),
                   bux=cost*w_vec)

  prob$qobj$i <- k_mat_as_triplet$i
  prob$qobj$j <- k_mat_as_triplet$j
  prob$qobj$v <- -k_mat_as_triplet$x

  capture.output(res_mosek <- mosek(prob))

  # Return the solution
  stopifnot(identical(res_mosek$response$code, 0))
  return(res_mosek$sol$itr$xx)
}
