
library(dplyr)

create_svhm <- function(df, covariates, cross_validation_val, cost_grid, gamma_squared, k=3, test_size=.2, varName_cencored="death", varName_futime = "futime", opt='osqp'){

  names(df)[names(df) == varName_cencored] <- "death"

  names(df)[names(df) == varName_futime] <- "futime"

  df[covariates] <- SVHM:::normalize(df, covariates)
  df <- transform(df,
                  death = as.logical(death),
                  id = 1:nrow(df))

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
  for (j in 1:length(cost_grid)) {
    cost <- cost_grid[[j]]
    mean_pearson_corr <- 0
    for (i in 1:cross_validation_val) {
      cat('####################\n',
          'RUNNING: ', j, i, '\n')
      training_set <- training_sets[[i]]
      validation_set <- validation_sets[[i]]

      model <- train_svhm(training_set, validation_set, covariates, cost, gamma_squared, k=k, opt=opt)
      mean_pearson_corr <- mean_pearson_corr + model$p_corr
    }

    mean_person_of_grid[j] <- mean_pearson_corr/cross_validation_val
  }

  best_cost_ind <- which.max(mean_person_of_grid)
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

  trained_model <- train_svhm(df_train, df_test, covariates, best_cost, gamma_squared, k, opt)
  trained_model['cost'] <- best_cost


  end_time <- Sys.time()
  time_train <- end_time - start_time

  cat('
    #######################################\n
    #            model results            #\n
    #######################################\n',
    'Time to find cost paramet er was ', time_cross_val, '\n',
    'Time of the training was ', time_train, '\n',
    'Total time was ', time_cross_val + time_train, '\n',
    'optimal costparamter is', best_cost, '\n',
    'The pearson correlation is', trained_model$p_corr, '\n'
    )

  return(trained_model)
}
