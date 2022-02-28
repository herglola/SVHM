
library(Rmosek)
library(osqp)
library(Matrix)

#' Train SVHM
#'
#' Uses the Rmosek or osqp package to train the time dependent SVHM on a given training and test set. The training calculates the risk scores and
#' the optimal decision function values for each individual at every event time of the training set.
#' The columns of the dataset must contain
#' \code{id, futime, death, covariates, lt, rt}
#' where \code{lt} and \code{rt} are the start and end times of each time interval. the \code{death} column must also be logical values.
#'
#'
#' @param train training dataset
#' @param test test dataset
#' @param cost cost parameter of the support vector machine of type numeric
#' @param opt which quadratic optimization is used (\code{opt='mosek'} or \code{opt='osqp'})
#' @param gamma_squared width of gaussian kernel
#'
#' @return {trained model with
#'          \code{$e_vec} vector indicating if an event happens at each event time
#'          \code{$sol} calculated optimal solution for each event time
#'          \code{$train} train dataset with risk scores
#'          \code{$test} test dataset with risk scores
#' }
#'
#'
#' @import Rmosek
#' @import osqp
#' @import Matrix
#'
#' @examples {
#' library(timereg)
#' library(SVHM)
#'
#' data(csl)
#'
#' df <- csl
#'
#' names(df)[names(df) == "dc"] <- "death"
#'
#' names(df)[names(df) == "eventT"] <- "futime"
#'
#' df <- transform(df,
#'                death = as.logical(death))
#'
#'
#' df<-split(df, df$id)
#'
#' df[sample(1:length(df))]
#'
#' partition <- SVHM:::createListPartition(df, 1, test_size=.3)
#'
#' df_test <- partition$"test"
#'
#' df_train <- partition[["1"]]
#'
#' trained_model <- train_time_svhm(df_train, df_test, c("sex"), 10, opt="osqp", gamma_squared=100)
#' }
#'
#' @export
train_time_svhm <-function(train, test, covariates, cost, opt='osqp', gamma_squared=.5){
  n <- length(train)

  train_join <- bind_rows(train)
  ordered_event_times <- data.frame('futime' = unique(sort(train_join$futime[train_join$death==TRUE])))
  print(ordered_event_times)
  num_event_times <- nrow(ordered_event_times)

  train_arr <- apply(expand.grid(1:nrow(ordered_event_times), 1:length(train)), MARGIN = 1,  function(x) {data_at_time(x[[2]], x[[1]], train, ordered_event_times)})
  train_arr <- aperm(array(unlist(train_arr, recursive = FALSE), c(ncol(train[[1]])+4, num_event_times, length(train))), c(3,1,2))
  dimnames(train_arr)[[2]] <- c(names(train[[1]]), 'Y', 'risk_at_j', 'risk_to_j', 'f')

  test_arr <- apply(expand.grid(1:nrow(ordered_event_times), 1:length(test)), MARGIN = 1,  function(x) {data_at_time(x[[2]], x[[1]], test, ordered_event_times)})
  test_arr <- aperm(array(unlist(test_arr, recursive = FALSE), c(ncol(test[[1]])+4, num_event_times, length(test))), c(3,1,2))
  dimnames(test_arr)[[2]] <- c(names(test[[1]]), 'Y', 'risk_at_j', 'risk_to_j', 'f')

  gamma_sol <- Matrix(rep(0, num_event_times*dim(train_arr)[[1]]), nrow=dim(train_arr)[[1]], ncol=num_event_times, sparse = TRUE)
  event_mat <- Matrix(rep(0, num_event_times*dim(train_arr)[[1]]), nrow=dim(train_arr)[[1]], ncol=num_event_times, sparse = TRUE)

  for (j in 1:num_event_times) {
    m<- train_arr[,,j]
    m<- matrix(m[m[,'Y'] == 1, ], ncol = ncol(m))
    dimnames(m)[[2]] <- c(names(train[[1]]), 'Y', 'risk_at_j', 'risk_to_j', 'f')
    m_covariates <- as(m[, covariates], 'matrix')

    opt_data <- SVHM:::optimization_time_data(matrix(as.numeric(m_covariates), nrow = nrow(m_covariates), ncol = ncol(m_covariates)), m, ordered_event_times$futime[j], gamma_squared=gamma_squared)

    kernel_mat <- opt_data$k_mat
    adap_kernel_mat <- opt_data$adap_k_mat
    event_vec <- opt_data$e_vec
    weight_vec <- opt_data$w_vec

    #############################
    # Loese Optimierungsproblem #
    #############################
    if (opt=='osqp'){
      gamma_sol_j <-  SVHM:::opt_time_sol_osqp(event_vec, adap_kernel_mat, weight_vec, cost)
    } else if (opt == 'mosek'){
      gamma_sol_j <-  SVHM:::opt_time_sol_mosek(event_vec, adap_kernel_mat, weight_vec, cost)
    } else {
      stop("Invalid optimization method!")
    }
    ind_at_risk <- which(train_arr[,'Y',j] == 1)
    gamma_sol[ind_at_risk,j] <- gamma_sol_j
    event_mat[ind_at_risk,j] <- event_vec

    f_at_j <- colSums(gamma_sol_j*event_vec*kernel_mat)
    risk_at_j <- risk_time_score_training(kernel_mat, event_vec, weight_vec, f_at_j, n)
    train_arr[ind_at_risk, 'risk_at_j', j] <- risk_at_j
    train_arr[ind_at_risk, 'risk_to_j', j] <- rowSums(cbind(unlist(train_arr[ind_at_risk, 'risk_to_j', j]), risk_at_j), na.rm=TRUE)
    train_arr[ind_at_risk, 'f', j] <- f_at_j



    ######################
    # Berechne Riskwerte #
    ######################
    m_test<- test_arr[,,j]
    m_test <- matrix(m_test[m_test[,'Y'] == 1, ], ncol = ncol(m_test))
    dimnames(m_test)[[2]] <- c(names(test[[1]]), 'Y', 'risk_at_j', 'risk_to_j', 'f')
    m_test_covariates <- matrix(m_test[, covariates], ncol = length(covariates))

    test_ind_at_risk <- which(test_arr[,'Y',j] == 1)
    if(length(test_ind_at_risk) != 0){
      test_risk_and_f <- unlist(apply(m_test_covariates, MARGIN = 1, function(x){risk_time_score(gamma_sol_j, event_vec, weight_vec, x, m_covariates, n, gamma_squared=gamma_squared, d=d)}), recursive = TRUE)
      test_f_at_j <- test_risk_and_f[c(TRUE,FALSE)]
      test_risk_at_j <- test_risk_and_f[c(FALSE,TRUE)]
      test_arr[test_ind_at_risk, 'risk_at_j', j] <- test_risk_at_j
      test_arr[test_ind_at_risk, 'risk_to_j', j] <- rowSums(cbind(unlist(test_arr[test_ind_at_risk, 'risk_to_j', j-1]),test_risk_at_j), na.rm=TRUE)
      test_arr[test_ind_at_risk, 'f', j] <- test_f_at_j
    }

  }


  return(list('e_vec' = event_mat,
              'sol' = gamma_sol,
              'train' = train_arr,
              'test' = test_arr
  ))
}
