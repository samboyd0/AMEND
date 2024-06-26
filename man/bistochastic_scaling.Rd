% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandomWalk.R
\name{bistochastic_scaling}
\alias{bistochastic_scaling}
\title{Bistochastic Scaling}
\usage{
bistochastic_scaling(trans_mat)
}
\arguments{
\item{trans_mat}{transition matrix}
}
\value{
A modified transition matrix that is approximately bistochastic and thus has been adjusted for degree bias.
}
\description{
Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
This is done by scaling the transition matrix to be approximately bistochastic (all row & column sums equal 1).
}
