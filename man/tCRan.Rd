\name{twoCrack}
\alias{twoCrack}
\title{The Two-parameter Crack Distribution}
\usage{
dTwoCrack(X,lambda,theta)
pTwoCrack(X,lambda,theta)
rTwoCrack(n,lambda,theta)
dBs(X,lambda,theta)
rBs(n,lambda,theta)
dIG(X,lambda,theta)
rIG(n,lambda,theta)
dLBIG(X,lambda,theta)
rLBIG(n,lambda,theta)
}
\arguments{
\item{X}{vector of data.}   
\item{n}{a number of observations.}   
\item{lambda}{a value of the parameter lambda.}   
\item{theta}{a value of the parameter theta.}   
}
\value{
\strong{Two-parameter Crack (TwoCrack) distribution}
\item{dTwoCrack}{gives the density of Two-parameter Crack distribution}
\item{pTwoCrack}{gives the distribution function of Two-parameter Crack distribution}
\item{rTwoCrack}{generates random numbers of Two-parameter Crack distribution}

\strong{Birnbaum-Saunders (BS) distribution}
\item{dBs}{gives the density of Birnbaum-Saunders distribution}
\item{rBs}{generates random numbers of Birnbaum-Saunders distribution}

\strong{Inverse Gaussian (IG) Distribution}
\item{dIG}{gives the density of Inverse Gaussian Distribution}
\item{rIG}{generates random numbers of Inverse Gaussian Distribution}

\strong{Length-Biased Inverse Gaussian (LBIG) distribution}
\item{dLBIG}{gives the density of Length-Biased Inverse Gaussian distribution}
\item{rLBIG}{generates random numbers of Length-Biased Inverse Gaussian distribution}
}
\description{
\strong{Inverse Gaussian (IG) Distribution}

Chikara and Folk studied the variables of the two-parameter inverse Gaussian distribution which is the
continuous probability distribution in 1978. This package is interested in the inverse Gaussian distribution
based on Ahmed et al. (2008). Suppose X is a random variable with an inverse Gaussian distribution.
Consequently, a density function can be written in this formula:

\code{fIG(x)=(lambda/(theta*sqrt(2*pi)))*((theta/x)^(3/2))*exp(-0.5*(sqrt(x/theta)-(lambda*sqrt(theta/x))))},
where lambda, theta>0.

The random numbers generation procedure of an inverse Gaussian distribution uses the composition
method (Ngamkham, 2019).

\strong{Length-Biased Inverse Gaussian (LBIG) distribution}

Define X as a positive random variable with the length-biased inverse gaussian distribution. A density
function can be written in this form:

\code{fLBIG=(1/(theta*sqrt(2*pi)))*((theta/X)^(1/2))*exp(-0.5*(sqrt(X/theta)-(lambda*sqrt(theta/X))))},
where lambda, theta>0.

Also, this package uses the composition method for the random numbers generation procedure proposed
by Ngamkham (2019).

\strong{Birnbaum-Saunders (BS) distribution}

The Birnbaum-Sanders distribution was first introduced by Birnbaum and Saunders (1969) as a
combination of inverse Gaussian distributions with a length-biased inverse Gaussian distribution. Later on,
Ahmed et al. (2008) introduced a new parametrization of the Birnbaum-Saunders distribution.
Consequently, the density function of the Birnbaum-Saunders distribution can be written in this form:

\code{fBS(x)=0.5*dIG(X,theta,lambda)+0.5*dLBIG(X,theta,lambda)},
where lambda, theta>0.

This package simulating random numbers of the Birnbaum-Sanders distribution with the composition
method (Ngamkham, 2019).

\strong{Two-parameter Crack (TwoCrack) distribution}

The two-parameter crack distribution is a positive skewness model, which is extensively used to model
failure times of fatiguing materials. This distribution proposed by Saengthong and Bodhisuwan (2014).
The two-parameter crack distribution has density:

\code{fTwoCrack(x)=(theta/(theta+1))*dIG(X,theta,lambda)+(1/(theta+1))*dLBIG(X,theta,lambda)},
where lambda, theta>0.

Hence, we recommend the acceptance-rejection method for the random numbers generation procedure
proposed by Phaphan and Pongsart (2019).
}

\examples{

\strong{the density function}
lambda <- 2
theta <- 2
X <- c(7,-2)
dIG(X,lambda,theta)
dLBIG(X,lambda,theta)
dBS(X,theta,lambda)
dTwoCrack(X,theta,lambda)

\strong{the random numbers generation}
n <- 30 #sample size
lambda <- 2
theta <- 2
X1 <- rBS(n,lambda,theta)
X2 <- rTwoCrack(n,lambda,theta)
X3 <- rIG(n,lambda,theta)
X4 <- rLBIG(n,lambda,theta)
}

\references{
Folks, J. L. & Chhikara, R. S. (1978). The inverse Gaussian distribution and its statistical application - a
review. Journal of the Royal Statistical Society.Serie B, 40, 263–289.

Ahmed, S.E., Budsaba, K., Lisawadi, S., & Volodin, A. (2008). Parametric estimation for the BirnbaumSaunders lifetime distribution based on new parametrization. Thailand Statistician, 6(2), 213-240.

Ngamkham, T. (2019). On the crack random numbers generation procedure. Lobachevskii Journal of
Mathematics, 40(8), 1204–1217.

Birnbaum, Z. W., & Saunders, S. C., A. (1969). New Family of Lifetime Distribution. Applied Probability, 6,
319-327. 

Saengthong, P. & Bodhisuwan, W. (2014). A new two-parameter crack distribution. Journal of Applied
Sciences, 14(8), 758-766.

Phaphan, W. & Pongsart, T. (2019). Estimating the Parameters of a Two-parameter Crack Distribution.
Applied Science and Engineering Progress, 12(1), 29-36.
}