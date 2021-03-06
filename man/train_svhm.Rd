% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_svhm.R
\name{train_svhm}
\alias{train_svhm}
\title{Train SVHM}
\usage{
train_svhm(
  train,
  test,
  covariates,
  cost,
  k = 3,
  opt = "osqp",
  gamma_squared = 0.5
)
}
\arguments{
\item{train}{training dataset}

\item{test}{test dataset}

\item{covariates}{vector of name of covariates}

\item{cost}{cost parameter of the support vector machine of type numeric}

\item{k}{integer of how many nearest event times are used to predict the event time (default is 3)}

\item{opt}{which quadratic optimization is used (\code{opt='mosek'} or \code{opt='osqp'})}

\item{gamma_squared}{width of gaussian kernel}
}
\value{
{trained model with
         \code{$e_vec} vector indicating vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m,
         \code{$k_mat} kernel matrix,
         \code{$sol} calculated optimal solution,
         \code{$t_predict} test dataset with risk scores \code{risk} and \code{t_predict},
         \code{$p_corr} pearson correlation of the predicted times
         \code{$C_indes} C-Index
}
}
\description{
Uses the Rmosek or osqp package to train the SVHM on a given training and test set. Names of the cencoring variable and event variable mus be \code{death} and \code{futime}
}
\note{
The mosek package requires a license
}
\examples{
{

library(KMsurv)
library(SVHM)

data(bmt)
df<-bmt[1:40,]

# shuffle data
rows <- sample(nrow(df))
df <- df[rows, ]

covariates <- c('z3', 'z4')

# censoring variable and event variable need to have names "death" and "futime"
names(df)[names(df) == "d3"] <- "death"
names(df)[names(df) == "t2"] <- "futime"

n<-floor(nrow(df)/2)
train<- df[(1:n), ]
test<- df[-(1:n), ]

train_svhm(train, test, covariates, 10, .5, k=1, opt='osqp')
}
}
