
library(Rmosek)
library(osqp)

train_svhm <-function(train, test, covariates, cost, gamma_squared, k=3, opt='osqp'){
  train <- transform(train, training_id = 1:nrow(train))

  #sortiert Datensatz nach ?berlebenszeit
  train <- train[order(train$futime),]
  train <- transform(train, Y = nrow(train):1)

  train_covariates <- subset( train, select = covariates )
  # Na Eintr?ge werden durch 0 ersetzt. Verf?lscht Ergebnis??? Besser Imputation??
  train_covariates[is.na(train_covariates)] <- 0

  ordered_event_times <- with(train,
                              data.frame(
                                futtime = sort(train$futime[train$death == TRUE]),
                                training_id = train$training_id[train$death == TRUE])
  )
  num_event_times <- nrow(ordered_event_times)

  #############################
  # Loese Optimierungsproblem #
  #############################

  # Erstellt Daten f?r Optimierungsproblem der Form max(gamma^t*risk_vec -1/2* gamma^t*kernel_mat*gamma) mit Nebenbedingungen lower_bound \leq cond_mat*gamma \leq upper_bound
  opt_data <- SVHM:::optimization_data(train_covariates, train, ordered_event_times, gamma_squared)
  kernel_mat <- opt_data[[5]]
  event_vec <- opt_data[[6]]

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
  train <- transform(train, g = risk_scores_training)


  test <- transform(test,
                    g  = sapply( 1:nrow(test), function(x) {
                                                SVHM:::risk_score(gamma_sol,
                                                event_vec,
                                                as.numeric(test[x, covariates]),
                                                train_covariates,
                                                num_event_times,
                                                gamma_squared)
                                                })
  )

  train_eventOnly <- train[(train$death==TRUE),]
  test <- transform(test,
                    t_predict  =  SVHM:::predict_event_time_vec(train_eventOnly, test$g, k)
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

  return(list('e_vec' = event_vec,
              'k_mat' = kernel_mat,
              'sol' = gamma_sol,
              't_predict' = test,
              'p_corr' = pearson_corr
              ))
}
