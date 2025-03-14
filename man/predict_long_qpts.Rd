% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jmwiv_functions.R
\name{predict_long_qpts}
\alias{predict_long_qpts}
\title{Predict longitudinal biomarker value at quadrature points}
\usage{
predict_long_qpts(
  standata,
  stanfit,
  times = qpts,
  ind,
  id_selected,
  psb,
  assoc_code
)
}
\arguments{
\item{standata}{Stan data input}

\item{times}{Quadrature times}

\item{ind}{Index for the individual}

\item{psb}{Dataset for individual}

\item{assoc_code}{Association code (0 for 'RE', 1 for 'LP', 2 for 'CV').}
}
\value{
Prediction of longitudinal biomarkers at quadrature points (to be used in \code{pred_plot_jmwiv}).
}
\description{
Internal function
}
