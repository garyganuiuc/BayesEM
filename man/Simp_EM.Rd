\name{Simp_EM}
\alias{Simp_EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	A simplified and faster version of Bayesian EM algorithm for the Sparsity learning in GGMs
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Simp_EM(S, n, p_n, a0, b0,alpha, beta, v0, v1, maxiter, w, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
  	Sample Covariance Matrix 
%%     ~~Describe \code{S} here~~
}
  \item{n}{
  	Number of observations
%%     ~~Describe \code{n} here~~
}
  \item{p_n}{
  	Size of the variables
%%     ~~Describe \code{p_n} here~~
}
  \item{a0}{
  	Hyperparameter on the Beta Prior
%%     ~~Describe \code{a0} here~~
}
  \item{b0}{
  	Hyperparameter on the Beta Prior
%%     ~~Describe \code{b0} here~~
}
  \item{alpha, beta}{
  	Hyperparameters on the Gamma Prior of tau
%%     ~~Describe \code{lambda} here~~
}
  \item{v0}{
  	The spike component in the spike and slab prior
%%     ~~Describe \code{v0} here~~
}
  \item{v1}{
  	The slab component in the spike and slab prior
%%     ~~Describe \code{v1} here~~
}
  \item{maxiter}{
  	Maximium rounds of iterations
%%     ~~Describe \code{maxiter} here~~
}
  \item{w}{
  	Hyperparameter on the Inverse Gamma Prior
%%     ~~Describe \code{w} here~~
}
  \item{l}{
  	Hyperparameter on the Inverse Gamma Prior
%%     ~~Describe \code{l} here~~
}
}

\value{
	\item{P}{Marginal posterior probabilty matrix of each entry being sparse}
	\item{Theta}{An MAP estimate of the precision matrix}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}