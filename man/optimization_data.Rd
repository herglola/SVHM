% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opt_data_func.R
\name{optimization_data}
\alias{optimization_data}
\title{Optimization Data}
\usage{
optimization_data(
  covariates,
  training_dataset,
  ordered_event_times,
  gamma_squared = 0.5,
  d = 1
)
}
\arguments{
\item{covariates}{dataset of covariates of the subjects in a dataset}

\item{training_dataset}{data frame representing the trainings dataset}

\item{ordered_event_times}{data frame of all event times ordered in ascending order}

\item{gamma_squared}{width of gaussian kernel}

\item{d}{degree of polynomial kernel}

\item{type}{Type of kernel, either 'gauss' or 'poly' for gaussian or polynomial kernel}
}
\value{
{List
\code{$r_vec}            vector representing at which event times the subjects are under risk
\code{$adap_kernel_mat}  matrix on which quadratic optimization will be performed
\code{$c_mat}            matrix representing the constraints of the optimization problem
\code{$w_vec}            vector of weights at any event time for all subjects
\code{$kernel_mat}       Gram matrix of covariates
\code{$e_vec}            vector indicating vector containing information if a subject is at risk or if an event happens. If n are the number of subjects and m the number of event times, then event_vec has length n*m,
}
}
\description{
calculates all needed values to execute quadratic optimization in SVHM
}
