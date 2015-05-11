\name{fdr}
\alias{fdr}
\title{
False Discovery Rate Computation
}
\description{
For internal use in functions \code{stp} and \code{bcaFDR}. 
Computes steps 2.1 to 2.4 from \emph{Algorithm 1} in the vignette.
}
\usage{
fdr(Z, design, th = NULL, B = 100, lambda = 0.5, PER = FALSE, ...)
}
\arguments{
\item{Z}{
a matrix or data.frame representing genes' expression levels. The rows 
of \eqn{Z} correspond to the genes in the experiment, and the columns 
correspond to the replicates. Treatment replicates are to the left, 
control replicates to the right.
}
\item{design}{
a vector of length equal to the number of columns in \code{Z} with 1's 
for the treatment replicates and 2's for the control replicates 
\eqn{(1, \ldots, 1, 2, \ldots, 2)}{(1, \ldots, 1, 2, \ldots, 2)}.
}
\item{th}{
Threshold values for estimating the FDR. If \code{NULL}, the values 
from \code{abs(ac2(Z,design))} are used.
}
\item{B}{
Number of bootstrap or permutation replications for estimating the FDR 
(as passed from \code{stp} and \code{bcaFDR}).
}
\item{lambda}{
Parameter for the estimation of \eqn{\pi_0}{pi0} and the FDR as passed 
from \code{stp} and \code{bcaFDR} (see Storey, 2002).
}
\item{PER}{
If \code{FALSE} (default), bootstrap replications are used to estimate 
the FDR. If \code{TRUE}, permutation replications are used instead 
(as passed from \code{stp} and \code{bcaFDR}). 
}
\item{\dots}{
additional arguments for parallel computation in \code{boot} function 
as passed from \code{stp} (see \code{\link{stp}} help page for details).
}
}
\value{
\item{Q }{
Estimations of the FDR using each value in \code{th} as the threshold.
}
\item{th }{
Threshold values used for estimating the FDR.
}
\item{pi0 }{
Estimation of \eqn{\pi_0}{pi0}, the true proportion of non differentially 
expressed genes in the experiment.
}
\item{B }{
Number of bootstrap or permutation replications used for estimating the FDR.
}
\item{lambda }{
Parameter used for the estimation of \eqn{\pi_0}{pi0} and the FDR.
}
\item{call }{The matched call.}
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification of 
Differentially Expressed Genes in Microarray Data}. Unpublished MS thesis. 
Universidad Nacional de Colombia, Bogot\'a.

Storey, J. D. (2002) \emph{A direct approach to false discovery rates.} 
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 
\bold{64(3): }479--498.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{stp}}.
}
\examples{
## Single time point analysis for 500 genes with 10 treatment 
## replicates and 10 control replicates
n <- 500; p <- 20; p1 <- 10
des <- c(rep(1, p1), rep(2, (p-p1)))
mu <- as.matrix(rexp(n, rate=1))
Z <- t(apply(mu, 1, function(mui) rnorm(p, mean=mui, sd=1)))
### 5 up regulated genes
Z[1:5,1:p1] <- Z[1:5,1:p1] + 5
### 10 down regulated genes
Z[6:15,(p1+1):p] <- Z[6:15,(p1+1):p] + 4

res <- fdr(Z, des)
plot(res$th, res$Q, type="l", col="blue")
legend(x="topright", legend="FDR", lty=1, col="blue")
}