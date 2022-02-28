

#' createDataPartition
#'
#' partitions a dataset into a test set and cross validation sets. createDataPartition() is not randomized, therefore df should be randomized before creating the partition!
#'
#' @param df data frame
#' @param cross_validation_val number of sets for k-fold cross validation
#' @param test_size size of test set (default=.8)
#'
#' @return {partitioned dataset
#' }
#'
#' @examples {
#' # Example with the preloaded mtcars dataset
#' df<-mtcars
#' partition <- SVHM:::createDataPartition(mtcars, 4, .2)
#' }
#'
createDataPartition <- function(df, cross_validation_val, test_size=.2){
  vals <- c(rep((1-test_size)/cross_validation_val, cross_validation_val), test_size)
  names <- c(c(1:cross_validation_val), 'test')
  spec <- setNames(vals,names)

  g = cut(
      seq(nrow(df)),
      nrow(df)*cumsum(c(0,spec)),
      labels = names(spec)
      )

  partition = split(df, g)
  return(partition)
}

#' createListPartition
#'
#' partitions a List into a test set and cross validation sets. createDataPartition() is not randomized, therefore the list should be randomized before creating the partition!
#'
#' @param l list of data
#' @param cross_validation_val number of sets for k-fold cross validation
#' @param test_size size of test set (default=.8)
#'
#' @return {partitioned list
#' }
#'
#' @examples {
#' # Example with the preloaded mtcars dataset
#' l<-list("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' partition <- SVHM:::createListPartition(l, 3, .1)
#' }
#'
createListPartition <- function(l, cross_validation_val, test_size=.2){
  vals <- c(rep((1-test_size)/cross_validation_val, cross_validation_val), test_size)
  names <- c(c(1:cross_validation_val), 'test')
  spec <- setNames(vals,names)

  g = cut(
    seq(length(l)),
    length(l)*cumsum(c(0,spec)),
    labels = names(spec)
  )

  partition = split(l, g)
  return(partition)
}
