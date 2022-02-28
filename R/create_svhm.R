
library(dplyr)
#' Train SVHM
#'
#' predicts the event times of a given dataset using cross validation for the cost parameter. Final model includes predicted event times as well as parameters to predict new event times from subject who are not in df.
#' All values of covariates are first normalized to the intervall [0,1] before the SVHM algorithm is applied. The cost parameter for the final model is chosen with the best pearson correlation.
#'
#'
#' @param df data frame
#' @param covariates vector of name of covariates
#' @param cross_validation_val number of subset to use for cost optimization
#' @param cost_grid grid of all cost parameter to be optoimzed uponl
#' @param varName_cencored name of variable in df that indicates cencoring
#' @param varName_futime name of variable in df that indicates event time
#' @param k integer of how many nearest event times are used to predict the event time (default is 3)
#' @param test_size size of final test set in precent
#' @param opt which quadratic optimization is used (\code{opt='mosek'} or \code{opt='osqp'})
#' @param gamma_squared width of gaussian kernel
#' @param choose optional parameter which decides if the C-index or the pearson correlation is used to determine the optimal cost parameter. Values are either \code{'c'} for the C-Index or \code{'p'} for the pearson correlation
#'
#' @return {trained model with
#'          \code{$e_vec} vector indicating vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m,
#'          \code{$k_mat} kernel matrix,
#'          \code{$sol} calculated optimal solution,
#'          \code{$t_predict} test dataset with risk scores \code{risk} and \code{t_predict},
#'          \code{$p_corr} pearson correlation of the predicted times
#'          \code{$C_indes} C-Index
#' }
#'
#' @note The mosek package requires a license
#'
#' @import dplyr
#'
#' @examples {
#'
#' library(KMsurv)
#' library(SVHM)
#'
#' ##############
#' # Parameters #
#' ##############
#'
#' gamma_squared <- 100
#' k <- 1
#' cross_validation_val <- 3
#' test_size=.3
#' cost_grid <- 2^c(-6:6)
#'
#' covariates <- c('z7')
#'
#' ######################
#' #  Model prediction  #
#' ######################
#'
#' data(bmt)
#'
#' model <- create_svhm(bmt, covariates, cross_validation_val, cost_grid, varName_cencored="d3", varName_futime = "t2", k=k, test_size=test_size, opt='osqp', gamma_squared=gamma_squared)
#' }
#'
#' @export
create_svhm <- function(df, covariates, cross_validation_val, cost_grid, varName_cencored, varName_futime, k=3, test_size=.2, opt='osqp', gamma_squared=.5, choose='c'){

  names(df)[names(df) == varName_cencored] <- "death"

  names(df)[names(df) == varName_futime] <- "futime"

  df <- transform(df,
                death = as.logical(death),
                id = 1:nrow(df))

  df[covariates] <- SVHM:::normalize(df, covariates)

  rows <- sample(nrow(df))
  df <- df[rows, ]

  partition <- SVHM:::createDataPartition(df, cross_validation_val, test_size=test_size)

  df_test <- partition$'test'

  cross_validation_sets <- partition[names(partition) != "test"]

  ################################################
  #        train model for cost parameter        #
  ################################################
  start_time <- Sys.time()


  training_sets <- list()
  validation_sets <- list()
  for (i in 1:cross_validation_val) {
    train_set <- bind_rows(cross_validation_sets[-i])
    training_sets[[i]] <- train_set
    validation_sets[i] <- cross_validation_sets[i]
  }

  mean_person_of_grid <-  list()
  mean_C_of_grid <- list()
  for (j in 1:length(cost_grid)) {
    cost <- cost_grid[[j]]
    mean_pearson_corr <- 0
    mean_C_index <- 0
    cat("Current cost paramter for which training is performed:", cost, "\n")
    for (i in 1:cross_validation_val) {
      training_set <- training_sets[[i]]
      validation_set <- validation_sets[[i]]
      model <- train_svhm(training_set, validation_set, covariates, cost, k=k, opt=opt, gamma_squared=gamma_squared)
      mean_pearson_corr <- mean_pearson_corr + model$p_corr

      mean_C_index <- mean_C_index + model$C_index$concordance
    }

    mean_person_of_grid[j] <- mean_pearson_corr/cross_validation_val
    mean_C_of_grid[j] <- mean_C_index/cross_validation_val
  }

  if (choose == 'p'){
    best_cost_ind <- which.max(mean_person_of_grid)
  } else if (choose == 'c'){
    best_cost_ind <- which.max(mean_C_of_grid)
  } else {
    stop('Wrong choice of choose parameter!')
  }

  best_cost <- cost_grid[best_cost_ind]


  end_time <- Sys.time()
  time_cross_val <- end_time - start_time

  cat('The best cost parameter is ', best_cost, '\n',
      'Total time for training the cost parameter was ', time_cross_val)

  ################################################
  #       train model for final prediction       #
  ################################################
  start_time <- Sys.time()


  df_train <- bind_rows(cross_validation_sets)

  trained_model <- train_svhm(df_train, df_test, covariates, best_cost, k, opt, gamma_squared=gamma_squared)
  trained_model['cost'] <- best_cost


  end_time <- Sys.time()
  time_train <- end_time - start_time

  cat('
    #######################################\n
    #            model results            #\n
    #######################################\n',
    'Time to find cost parameter was ', time_cross_val, '\n',
    'Time of the training was ', time_train, '\n',
    'Total time was ', time_cross_val + time_train, '\n',
    'optimal costparamter is', best_cost, '\n',
    'The pearson correlation is', trained_model$p_corr, '\n',
    'The C-Index is', trained_model$C_index$concordance ,'\n'
    )

  return(trained_model)
}
