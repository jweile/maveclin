% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{calc.posterior}
\alias{calc.posterior}
\title{Calculate posterior probability via Bayesian Inference}
\usage{
calc.posterior(llrs, prior)
}
\arguments{
\item{llrs}{Vector of log likelihood ratios (i.e. log(Bayes-Factor))}

\item{prior}{The prior probability}
}
\value{
the posterior probability
}
\description{
Calculate posterior probability via Bayesian Inference
}
