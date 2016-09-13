\name{Q}
\alias{Q}
\alias{Q4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculating the Q function
%%  ~~function to do ... ~~
}
\description{
	Calculating the Q function. Q function is used to ensure updating in the "right direction", i.e the objective function increases, for every update.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Q(Theta, tau, alpha, beta,  a0, b0, v0, v1, p, P, sigma2, id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Theta}{
  	An estimate of the precision matrix
%%     ~~Describe \code{Theta} here~~
}
  \item{tau}{
  	An MAP estimate of tau
%%     ~~Describe \code{tau} here~~
}
  \item{alpha, beta}{
  	An hyperparamter of the Exponential Prior on tau
%%     ~~Describe \code{lambda} here~~
}
  \item{a0}{
  	An hyperparameter of the Beta Prior 
%%     ~~Describe \code{a0} here~~
}
  \item{b0}{
  	An hyperparameter of the Beta Prior
%%     ~~Describe \code{b0} here~~
}
  \item{v0}{
  	The spike component in the spike and slab prior
%%     ~~Describe \code{v0} here~~
}
  \item{v1}{
  	The slab component in the spike and slab prior
%%     ~~Describe \code{v1} here~~
}
  \item{sigma2}{
  	The component that controls the scale in the spike and slab prior
%%     ~~Describe \code{sigma2} here~~
}

  \item{p}{
  	An MAP estimate on p
%%     ~~Describe \code{p} here~~
}
  \item{P}{
  	Marginal posterior probabilty matrix of each entry being sparse
%%     ~~Describe \code{P} here~~
}
  \item{id}{
  	The id of the upper triangular of the matrix
%%     ~~Describe \code{id} here~~
}
}

\value{
	\item{Q}{The value of Q function, the objective function}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
