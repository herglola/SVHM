
#' @title Constraint matrix in quadratic optimization problem
#'
#' calculates the matrix which defines the constraints in the SVHM algorithm
#'
#'
#' @param event_vec vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m
#' @param num_event_time number of event
#'
#' @return matrix
#'
#' @import Matrix
#'

library(Matrix)

condition_mat <- function(event_vec, num_event_time) {
  m <- Matrix(0, nrow = length(event_vec), ncol = num_event_time, sparse=TRUE)

  for (i in 1:num_event_time) {
    col_vec <- event_vec
    col_vec[c(rep(TRUE, i-1), FALSE, rep(TRUE, num_event_time-i))] <- 0
    m[,i] <- col_vec
  }
  a <- rbind(Diagonal(length(event_vec)), t(m))
  return(Matrix(a,sparse = TRUE))
}
