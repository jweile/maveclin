% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{flipScores}
\alias{flipScores}
\title{Flip transformation}
\usage{
flipScores(xs)
}
\arguments{
\item{xs}{numerical vector}
}
\value{
the flipped scores as a numerical vector
}
\description{
Applies the "flip" transformation to a numerical vector. I.e. all values 
greater than one are transformed by 1/x
}
