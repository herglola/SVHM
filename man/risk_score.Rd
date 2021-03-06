% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risk_score_func.r
\name{risk_score}
\alias{risk_score}
\title{risk scores}
\usage{
risk_score(
  gamma_sol,
  event_vec,
  v,
  covariates_train,
  num_event_times,
  gamma_squared = 0.5,
  d = 1
)
}
\arguments{
\item{gamma_sol}{optimal solution of the SVHM}

\item{event_vec}{vector containing information of the training if a subject in the training dataset is at risk or if an event happens.
If n are the number of subjects in the training dataset and m the number of event times in the training dataset, then event_vec has length n*m}

\item{v}{covariates of the individual for which the risk is to be calculated}

\item{covariates_train}{dataset of covariates of the subjects in the training dataset}

\item{num_event_times}{number of event times that occour in the training data set}

\item{gamma_squared}{width of gaussian kernel}

\item{d}{degree of polynomial kernel}

\item{type}{Type of kernel, either 'gauss' or 'poly' for gaussian or polynomial kernel}
}
\value{
risk score of the individual
}
\description{
calculates the risk scores for one individual with the help of the calculated optimal solution to the quadratic programming problem of SVHM and the kernel matrix of the covariates of the test dataset.
}
\note{
The calculated risk score is not the actual risk scores defined by the Risk function but it induce an ordering of the risk scores. For detailed information see reference
}
\references{
Wang, Y., Chen, T., and Zeng, D. Support vector hazards machine: A counting process framework for learning risk scores for censored outcomes. Journal of Machine Learning Research, 17(167):1-37, 2016
}
