\name{qval}
\alias{qval}
\title{
Q-Values Computation
}
\description{
For internal use in function \code{stp}. Computes the genes' Q-Values 
in the Single Time Point Analysis according to \emph{Algorithm 3} 
in the vignette.
}
\usage{
qval(Q, psi2)
}
\arguments{
\item{Q}{
vector with the estimated FDRs when the threshold values 
used are \code{abs(ac2(Z, design))}.
}
\item{psi2}{
vector with the second artificial component as returned by \code{ac2}.
}
}
\value{
returns a vector with the computed Q-Values for each gene in the experiment.
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification 
of Differentially Expressed Genes in Microarray Data}. Unpublished MS 
thesis. Universidad Nacional de Colombia, Bogot\'a.

Storey, J. D. (2002) \emph{A direct approach to false discovery rates.} 
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 
\bold{64(3): }479--498.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{stp}}, \code{\link{fdr}}, \code{\link{ac2}}.
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
qValues <- qval(res$Q, ac2(Z, des))
plot(res$th, res$Q, type="l", col="blue")
lines(res$th, qValues[order(abs(ac2(Z, des)))], col="green")
legend(x="topright", legend=c("FDR", "Q Values"), lty=c(1,1), 
    col=c("blue", "green"))
}