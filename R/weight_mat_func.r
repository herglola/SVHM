
#' Weight matrix
#'
#' calculates the weights for every individual in the training dataset at every event time. If an individual experiences an event the weight is given by the ratio of at risk subjects with no event to all at risk subjects.
#' If no event is experienced the weight is given by the ratio of one over all at risk subjects.
#'
#' @param training_dataset data frame representing the training data
#' @param ordered_event_times data frame of all event times ordered in ascending order
#'
#' @return matrix storing all weights for every individual
#'
#'

weight_mat <- function(training_dataset, ordered_event_times) {
  weight_matrix <- matrix(nrow = nrow(training_dataset),
                          ncol = length(training_dataset$futime[training_dataset$death == TRUE]),
                    )
  for(row in 1:nrow(weight_matrix)){
    for (col in 1:ncol(weight_matrix)){
      if ((training_dataset$futime[row] == ordered_event_times$futime[col]) & (training_dataset$death[row]==TRUE)){
        weight_matrix[row,col] <- 1-1/(training_dataset$Y[row])
      } else {
        weight_matrix[row,col] <- 1/(training_dataset$Y[row])
      }
    }
  }
  return (Matrix(weight_matrix, sparse=TRUE))
}
