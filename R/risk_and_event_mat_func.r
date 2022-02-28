library(Matrix)


#' Risk and Event Matrix
#'
#' calculates two matrices of length n*m, if n are the number of subjects and m the number of event times.
#' The Risk Matrix indicates for every subject in a dataset if the subject is still at risk at every event time.
#' The Event Matrix is equal to the Risk Matrix but if a subject experiences an event at an event time the entrie is set to -1
#'
#'
#' @param training_dataset data frame representing the training data
#' @param ordered_event_times data frame of all event times ordered in ascending order
#'
#' @return {List
#' \code{$r_mat}  matrix indicating at risk at every event time for subjects
#' \code{$e_mat}  matrix indicating if subjects are at risk and if they are experiencing an event at any event time
#' }
#'
#' @import Matrix
#'
#' @examples {
#' # Create random data
#' train <- data.frame(futime = sample.int(10,6),
#'                     death = sample(c(TRUE,FALSE), 6, replace=TRUE),
#'                     training_id=1:6)
#' ordered_event_times <- with(train,
#'                            data.frame(
#'                               futime = sort(train$futime[train$death == TRUE]),
#'                               training_id = train$training_id[train$death == TRUE])
#' )
#'
#' SVHM:::create_risk_and_evemt_matrix(train, ordered_event_times)
#' }
#'

create_risk_and_event_matrix <- function(training_dataset, ordered_event_times) {
  risk_matrix <- Matrix(0, nrow = nrow(training_dataset),
                        ncol = nrow(ordered_event_times),
                        sparse = TRUE
  )
  event_matrix <- Matrix(0, nrow = nrow(training_dataset),
                         ncol = nrow(ordered_event_times),
                         sparse = TRUE
  )
  for(row in 1:nrow(risk_matrix)){
    for (col in 1:ncol(risk_matrix)){
      if (training_dataset$futime[row] >= ordered_event_times$futime[col]) {
        risk_matrix[row,col] <- 1
        if ((training_dataset$futime[row] == ordered_event_times$futime[col]) & (training_dataset$death[row]==TRUE)){
          event_matrix[row,col] <- 1
        } else {
          event_matrix[row,col] <- -1
        }
      }
    }
  }
  return(list('r_mat'=risk_matrix,
              'e_mat'=event_matrix))
}
