
#' Predict Event Time of a subject
#'
#' calculate the predicted event time of an individual
#'
#' This function predicts the event time of a subject based on the k closest risks subjects
#' in the training dataset of non cencored individuals df. The risks in df are ranked and the
#' predicted event time is the average of the k event times that coincide with the rank of the
#' k closest risks. The predicted event time is rounded up to integers by default.
#' A vectorized version predict_event_time_vec() for the parameter x exists.
#'
#' @param df dataframe of non censored subjects in the training set
#' @param x  Risk score of the individual which will be predicted upon
#' @param k  integer of how many nearest event times are used to predict the event time (default is 3)
#' @param rounding Options are 'ceil', 'floor' and 'no'. (default is 'ceil')
#'
#' @return predicted event time
#'
#' @references Wang, Y., Chen, T., and Zeng, D. Support vector hazards machine: A counting process framework for learning risk scores for censored outcomes. Journal of Machine Learning Research, 17(167):1-37, 2016
#'
#'
predict_event_time <- function(df, x, k=3, rounding='ceil') {
  stopifnot(rounding=='ceil' | rounding == 'floor' | rounding == 'no')
  event_times <- sort(df$futime)
  df1 <- transform(df, dist = abs(df$g-x))
  df1 <- df1[order(df1$g, decreasing = TRUE),]
  df1 <- transform(df1, rank = rank(df1$dist, ties.method = 'first'))
  rank_indices <- which(df1$rank < k+1)
  predicted_time <- sum(event_times[rank_indices])/length(rank_indices)
  if (rounding == 'ceil'){
    predicted_time <- ceiling(predicted_time)
  } else if (rounding == 'floor'){
    predicted_time <- floor(predicted_time)
  }
  return(predicted_time)
}

predict_event_time_vec <- Vectorize(predict_event_time, c("x"))
