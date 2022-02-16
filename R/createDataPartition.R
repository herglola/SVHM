

#' createDataPartition
#'
#' partitions a dataset into a train set and cross validation sets. createDataPartition() is not randomized, therefore df should be randomized before creating the partition!
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
#' partition <- createDataPartition(mtcars, 4, .2)
#' }
#'
createDataPartition <- function(df, corss_validation_val, test_size=.2){
  vals <- c(rep((1-test_size)/corss_validation_val, corss_validation_val), test_size)
  names <- c(c(1:corss_validation_val), 'test')
  spec <- setNames(vals,names)

  g = cut(
      seq(nrow(df)),
      nrow(df)*cumsum(c(0,spec)),
      labels = names(spec)
      )

  partition = split(df, g)
  return(partition)
}
