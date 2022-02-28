
#' data_at_time
#'
#' Retrieves the relevant data for each individual at every event time.
#'
#' @param i index of individuals
#' @param j index of event time
#' @param df list of dataframes of the individuals
#' @param times dataframe of event times
#'
#' @return row of data for i-th individual if the individual is still under risk, otherwise return row of NA.
#'
#'
data_at_time<-function(i,j,df,times){
  v<-df[[i]]$'lt'
  if(df[[i]]$futime[1] < times$futime[j]){
    ind<-which.max(v[v<=times$futime[j]])
    r <- c(df[[i]][ind,], 0, NA, NA, NA)
  } else{
    ind<-which.max(v[v<=times$futime[j]])
    r <- c(df[[i]][ind,], 1, NA, NA, NA)
  }
  return(r)
}
