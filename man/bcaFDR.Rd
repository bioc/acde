\name{bcaFDR}
\alias{bcaFDR}
\title{
BCa Confidence Upper Bound for the FDR
}
\description{
For internal use in function \code{stp}. Computes a BCa confidence 
upper bound for the FDR following \emph{Algorithm 2} in the vignette.
}
\usage{
bcaFDR(Z, design, th = NULL, B = 100, 
    lambda = 0.5, PER = FALSE, R = 1000, 
    gamma = 0.95, Q = NULL, ...)
}
\arguments{
\item{Z}{
a matrix or data.frame representing genes' expression levels. The rows of 
\eqn{Z} correspond to the genes in the experiment, and the columns 
correspond to the replicates. Treatment replicates are to the left, 
control replicates to the right.
}
\item{design}{
a vector of length equal to the number of columns in \code{Z} with 1's for 
the treatment replicates and 2's for the control replicates \eqn{(1, \ldots, 
1, 2, \ldots, 2)}{(1, \ldots, 1, 2, \ldots, 2)}.
}
\item{th}{
Threshold values for estimating the FDR. If \code{NULL}, the values from 
\code{abs(ac2(Z,design))} are used.
}
\item{B}{
Number of bootstrap or permutation replications for estimating the FDR at 
each iteration (as passed from \code{stp}).
}
\item{lambda}{
Parameter for the estimation of \eqn{\pi_0}{pi0} and the FDR as passed 
from \code{stp} (see Storey, 2002).
}
\item{R}{
Number of bootstrap replications for the computation of the FDR's BCa 
confidence upper bound (as passed from \code{stp}).
}
\item{gamma}{
Confidence level for the FDR's BCa upper confidence bound (as passed 
from \code{stp}).
}
\item{PER}{
If \code{FALSE} (default), bootstrap replications are used to estimate 
the FDR. If \code{TRUE}, permutation replications are used instead 
(as passed from \code{stp}). 
}
\item{Q}{
Estimated FDR as returned in object \code{\$Q} from \code{fdr} function 
(passed from call to \code{stp}). For internal use. 
}
\item{\dots}{
additional arguments for parallel computation in \code{boot} function as 
passed from \code{stp} (see \code{\link{stp}} help page for details).
}
}
\value{
\item{cbound }{
BCa upper confidence bound for the FDR for each threshold value in \code{th}.
}
\item{warnings }{
warning messages generated from use of \code{boot.ci} function from package 
\code{boot}.
}
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification of 
Differentially Expressed Genes in Microarray Data}. Unpublished MS thesis. 
Universidad Nacional de Colombia, Bogot\'a.

Storey, J. D. (2002) \emph{A direct approach to false discovery rates.} 
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 
\bold{64(3): }479--498.

Efron B. and Tibshirani R. J. (1994) \emph{An Introduction to the Bootstrap.} 
Chapman & Hall/CRC, 1993.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{stp}}.
}
\examples{
## Single time point analysis for 50 genes with 10 treatment 
## replicates and 10 control replicates
n <- 50; p <- 20; p1 <- 10
des <- c(rep(1, p1), rep(2, (p-p1)))
mu <- as.matrix(rexp(n, rate=1))
Z <- t(apply(mu, 1, function(mui) rnorm(p, mean=mui, sd=1)))
### 5 up regulated genes
Z[1:5,1:p1] <- Z[1:5,1:p1] + 5
### 10 down regulated genes
Z[6:15,(p1+1):p] <- Z[6:15,(p1+1):p] + 5

resFdr <- fdr(Z, des)
bca <- bcaFDR(Z, des, Q=resFdr$Q, B=50, R=500)
plot(resFdr$th, resFdr$Q, type="l", col="blue")
lines(resFdr$th, bca$cbound, col="green")
legend(x="topright", legend=c("FDR", "BCa upper bound"), 
    lty=c(1,1), col=c("blue", "green"))
## Note: Discontinuities in the BCa upper bound are due to warnings
## generated during computations with function \code{boot.ci} 
## from package \code{boot}.
}