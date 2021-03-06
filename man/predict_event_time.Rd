% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_event_time.R
\name{predict_event_time}
\alias{predict_event_time}
\title{Predict Event Time of a subject}
\usage{
predict_event_time(df, x, k = 3, rounding = "ceil")
}
\arguments{
\item{df}{dataframe of non censored subjects in the training set}

\item{x}{Risk score of the individual which will be predicted upon}

\item{k}{integer of how many nearest event times are used to predict the event time (default is 3)}

\item{rounding}{Options are 'ceil', 'floor' and 'no'. (default is 'ceil')}
}
\value{
predicted event time
}
\description{
calculate the predicted event time of an individual
}
\details{
This function predicts the event time of a subject based on the k closest risks subjects
in the training dataset of non cencored individuals df. The risks in df are ranked and the
predicted event time is the average of the k event times that coincide with the rank of the
k closest risks. The predicted event time is rounded up to integers by default.
A vectorized version predict_event_time_vec() for the parameter x exists.
}
\references{
Wang, Y., Chen, T., and Zeng, D. Support vector hazards machine: A counting process framework for learning risk scores for censored outcomes. Journal of Machine Learning Research, 17(167):1-37, 2016
}
