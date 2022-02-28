
library(dplyr)
#' Train Time Dependent SVHM
#'
#' Calculates the Risk score and the value of the prediction function for each individual in the data set.
#'
#' @note In contrast to the \code{create_svhm()} function this function does not predict event times!
#'
#' @param df data frame
#' @param covariates vector of name of covariates
#' @param cost cost parameter to be used
#' @param varName_cencored name of variable in df that indicates cencoring
#' @param varName_futime name of variable in df that indicates event time
#' @param start_interval name of variable that indicates when the interval starts
#' @param end_interval name of variable that indicates when the interval ends
#' @param test_size size of final test set in precent
#' @param opt which quadratic optimization is used (\code{opt='mosek'} or \code{opt='osqp'})
#' @param gamma_squared width of gaussian kernel
#'
#' @return {trained model with
#'          \code{$e_vec} vector indicating if an event happens at each event time
#'          \code{$sol} calculated optimal solution for each event time
#'          \code{$train} train dataset with risk scores
#'          \code{$test} test dataset with risk scores
#'          \code{cost} cost parameter
#' }
#'
#'
#' @import dplyr
#'
#' @examples {
#'
#' library(timereg)
#' library(SVHM)
#'
#' ##############
#' # Parameters #
#' ##############
#'
#' opt <- "osqp"
#' gamma_squared <- 200
#' test_size=.3
#' cost <- 16
#'
#' ######################
#' #  Model prediction  #
#' ######################
#'
#' data(csl)
#'
#' time_model <- create_time_svhm(csl, c("sex", "age"), cost, varName_cencored='dc', varName_futime='eventT', start_interval='lt, end_interval='rt, test_size=test_size, opt=opt, gamma_squared=gamma_squared)
#' }
#'
#' @export
create_time_svhm <- function(df, covariates, cost, varName_cencored, varName_futime, start_interval, end_interval, test_size=.3, opt='osqp', gamma_squared=.5){

  names(df)[names(df) == varName_cencored] <- "death"

  names(df)[names(df) == varName_futime] <- "futime"

  names(df)[names(df) == start_interval] <- "lt"

  names(df)[names(df) == end_interval] <- "rt"

  df <- transform(df,
                  death = as.logical(death))

  df[covariates] <- SVHM:::normalize(df, covariates)

  df<-split(df, df$id)

  df[sample(1:length(df))]

  partition <- SVHM:::createListPartition(df, 1, test_size=test_size)

  df_test <- partition$'test'

  df_train <- partition[['1']]

  ################################################
  #         train model for risk scores          #
  ################################################
  start_time <- Sys.time()


  trained_model <- train_time_svhm(df_train, df_test, covariates, cost, opt, gamma_squared=gamma_squared)
  trained_model['cost'] <- cost


  end_time <- Sys.time()
  time_train <- end_time - start_time

  cat('
    #######################################\n
    #            model results            #\n
    #######################################\n',
      'Calculation finished! \n',
      'The total time was ', time_train, '\n'
  )

  return(trained_model)
}
