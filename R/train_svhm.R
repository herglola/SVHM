
library(Rmosek)
library(osqp)
library(survival)

#' Train SVHM
#'
#' Uses the Rmosek or osqp package to train the SVHM on a given training and test set. Names of the cencoring variable and event variable mus be \code{death} and \code{futime}
#'
#'
#' @param train training dataset
#' @param test test dataset
#' @param covariates vector of name of covariates
#' @param cost cost parameter of the support vector machine of type numeric
#' @param k integer of how many nearest event times are used to predict the event time (default is 3)
#' @param opt which quadratic optimization is used (\code{opt='mosek'} or \code{opt='osqp'})
#' @param gamma_squared width of gaussian kernel
#'
#' @return {trained model with
#'          \code{$e_vec} vector indicating vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m,
#'          \code{$k_mat} kernel matrix,
#'          \code{$sol} calculated optimal solution,
#'          \code{$t_predict} test dataset with risk scores \code{risk} and \code{t_predict},
#'          \code{$p_corr} pearson correlation of the predicted times
#'          \code{$C_indes} C-
#' }
#'
#' @note The mosek package requires a license
#'
#' @import Rmosek
#' @import osqp
#' @import survival
#'
#' @examples {
#'
#' library(KMsurv)
#' library(SVHM)
#'
#' data(bmt)
#' df<-bmt[1:40,]
#'
#' # shuffle data
#' rows <- sample(nrow(df))
#' df <- df[rows, ]
#'
#' covariates <- c('z3', 'z4')
#'
#' # censoring variable and event variable need to have names "death" and "futime"
#' names(df)[names(df) == "d3"] <- "death"
#' names(df)[names(df) == "t2"] <- "futime"
#'
#' n<-floor(nrow(df)/2)
#' train<- df[(1:n), ]
#' test<- df[-(1:n), ]
#'
#' train_svhm(train, test, covariates, 10, .5, k=1, opt='osqp')
#' }
#' @export
train_svhm <-function(train, test, covariates, cost, k=3, opt='osqp', gamma_squared=.5){
  train <- transform(train, training_id = 1:nrow(train))

  #sortiert Datensatz nach ?berlebenszeit
  train <- train[order(train$futime),]
  train <- transform(train, Y = nrow(train):1)

  train_covariates <- subset( train, select = covariates )
  # Na Eintr?ge werden durch 0 ersetzt. Verf?lscht Ergebnis??? Besser Imputation??
  train_covariates[is.na(train_covariates)] <- 0

  ordered_event_times <- with(train,
                              data.frame(
                                futime = sort(train$futime[train$death == TRUE]),
                                training_id = train$training_id[train$death == TRUE])
  )
  num_event_times <- nrow(ordered_event_times)

  #############################
  # Loese Optimierungsproblem #
  #############################

  # Erstellt Daten f?r Optimierungsproblem der Form max(gamma^t*risk_vec -1/2* gamma^t*kernel_mat*gamma) mit Nebenbedingungen lower_bound \leq cond_mat*gamma \leq upper_bound
  opt_data <- SVHM:::optimization_data(train_covariates, train, ordered_event_times, gamma_squared=gamma_squared)

  kernel_mat <- opt_data$k_mat
  event_vec <- opt_data$e_vec

  if (opt=='osqp'){
    gamma_sol <-  SVHM:::opt_sol_osqp(opt_data, num_event_times, cost)
  } else if (opt == 'mosek'){
    gamma_sol <-  SVHM:::opt_sol_mosek(opt_data, num_event_times, cost)
  } else {
    stop("Invalid optimization method!")
  }

    ######################
  # Berechne Riskwerte #
  ######################
  risk_scores_training <-  SVHM:::risk_score_training(gamma_sol, kernel_mat, event_vec, num_event_times, trainig_set_size)
  train <- transform(train, risk = risk_scores_training)

  test <- transform(test,
                    risk  = sapply( 1:nrow(test), function(x) {
                                                SVHM:::risk_score(gamma_sol,
                                                event_vec,
                                                as.numeric(test[x, covariates]),
                                                train_covariates,
                                                num_event_times,
                                                gamma_squared=gamma_squared,
                                                d=d)
                                                })
  )

  train_eventOnly <- train[(train$death==TRUE),]
  test <- transform(test,
                    t_predict  =  SVHM:::predict_event_time_vec(train_eventOnly, test$risk, k)
  )

  xBar <- sum(test$futime)/nrow(test)
  yBar <- sum(test$t_predict)/nrow(test)
  up <- sum((test$futime-xBar)*(test$t_predict-yBar))
  downx <- sqrt(sum((test$futime-xBar)^2))
  downy <- sqrt(sum((test$t_predict-yBar)^2))
  if (downy == 0){
    warning('Pearson correlation not defined and was artificially set to 0')
    pearson_corr <- 0
  } else{
    pearson_corr <- up/(downx*downy)
  }

  C_statistic <- survival::concordancefit(test$futime, test$risk, timewt="n")

  return(list('e_vec' = event_vec,
              'k_mat' = kernel_mat,
              'sol' = gamma_sol,
              't_predict' = test,
              'p_corr' = pearson_corr,
              'C_index' = C_statistic
              ))
}
