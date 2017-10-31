\name{EM}
\alias{EM}
\title{
	EM algorithm
}
\description{
The EM algorithm
}
\usage{
EM(S, n, p_n, alpha, beta, v0, v1, maxiter, Ra)
}

\arguments{
  \item{S}{
  	Sample Covariance Matrix

}
  \item{n}{
  	Sample Size
}
  \item{p_n}{
  	Variable Size
}
  \item{alpha}{
  	Tuning Parameter alpha
}
  \item{beta}{
  	Tuning Parameter beta
}
  \item{v0}{
  	Spike and slab prior v0
}
  \item{v1}{
  	Spike and slab prior v1
}
  \item{maxiter}{
  	Max iteration
}
  \item{Ra}{
  	Constraint on Matrix Range
}
}
