\name{stp}
\alias{stp}
\title{
Single Time Point Analysis for Detecting Differentially Expressed Genes
}
\description{
Performs the Single Time Point Analysis for detecting differentially 
expressed genes following Acosta (2015).  
}
\usage{
stp(Z, design, alpha = 0.05, B = 100, lambda = 0.5, 
    th = NULL, PER = FALSE, BCa = FALSE, gamma = 0.95, 
    R = 1000, ...)
}
\arguments{
\item{Z}{
a matrix or data.frame representing genes' expression levels. The 
rows of \eqn{Z} correspond to the genes in the experiment, and the 
columns correspond to the replicates. Treatment replicates are to 
the left, control replicates to the right.
}
\item{design}{
a vector of length equal to the number of columns in \code{Z} with 
1's for the treatment replicates and 2's for the control replicates 
\eqn{(1, \ldots, 1, 2, \ldots, 2)}{(1, \ldots, 1, 2, \ldots, 2)}.
}
\item{alpha}{
between 0 and 1. Desired level for controlling the false discovery 
rate (FDR).
}
\item{B}{
Number of bootstrap or permutation replications for estimating the FDR.
}
\item{lambda}{
Parameter for the estimation of \eqn{\pi_0}{pi0} and of the FDR 
(see Storey, 2002).
}
\item{th}{
Threshold values for estimating the FDR. If \code{NULL}, the values 
from \code{abs(ac2(Z,design))} are used.
}
\item{PER}{
If \code{FALSE} (default), bootstrap replications are used to estimate 
the FDR. If \code{TRUE}, permutation replications are used instead. 
}
\item{BCa}{
If \code{TRUE}, a BCa confidence upper bound for the FDR is computed 
(see Efron and Tibshirani, 1994).
}
\item{gamma}{
Confidence level for the FDR's BCa confidence upper bound.
}
\item{R}{
Number of bootstrap replications for the computation of the FDR's 
BCa confidence upper bound.
}
\item{\dots}{
additional arguments for parallel computation in \code{boot} 
function (see Details).
}
}
\details{
For details on the computations performed in this function, 
see Acosta (2015).

Additional parameters in the '\code{...}' argument are used for 
parallel computation in bootstrap calculations. These are supplied 
to calls to the \code{boot} function in package \code{boot}. With 
this in mind, the use of additional arguments must be restricted to 
arguments \code{parallel} and \code{ncpus} from function \code{boot}.
}
\section{Warning }{
If argument \code{BCa=TRUE}, computations may take a considerable 
amount of time.
}
\value{
\code{stp} returns an object of class '\code{STP}', which is a 
list with components:

\item{dgenes }{
factor with the classification of each gene in \code{Z}. Classes: 
"up-reg.", "down-reg.", "no-diff.".
}
\item{tstar }{
Threshold value used to identify differentially expressed genes.
}
\item{astar }{
Achieved FDR level.
}
\item{Q }{Estimations of the FDR using each value in \code{th} as 
threshold.
}
\item{th }{
Threshold values used for estimating the FDR.
}
\item{qvalues }{
Estimated Q-Values for the genes in the analysis. If argument 
\code{th!=NULL}, these are not computed.
}
\item{pi0 }{
Estimation of \eqn{\pi_0}{pi0}, the true proportion of non 
differentially expressed genes in the experiment.
}
\item{B }{
Number of bootstrap or permutation replications used for estimating 
the FDR.
}
\item{lambda }{
Parameter used for the estimation of \eqn{\pi_0}{pi0} and the FDR.
}
\item{ac }{
Artificial components of \code{Z}.
}
\item{gNames }{
Gene names (by default the row names in \code{Z}).
}
\item{iRatio }{
Inertia ratio \eqn{Var(\psi_2) / \lambda_1}{Var(\psi[2]) / \lambda[1]}, 
where \eqn{\lambda_1}{\lambda[1]} is the first eigenvalue of \code{Z}'s 
Principal Components Analysis.
}
\item{bca }{
BCa upper confidence bounds for the FDR using each value in \code{th} 
as the threshold.
}
\item{gamma }{
Confidence level used in the computation of the BCa upper bounds.
}
\item{alpha }{
Desired FDR level.
}
\item{call }{
The matched call.
}
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification 
of Differentially Expressed Genes in Microarray Data}. Unpublished MS 
thesis. Universidad Nacional de Colombia, Bogot\'a.

Storey, J. D. (2002) \emph{A direct approach to false discovery rates.} 
Journal of the Royal Statistical Society: Series B (Statistical 
Methodology), \bold{64(3): }479--498.

Efron B. and Tibshirani R. J. (1994) \emph{An Introduction to the 
Bootstrap.} Chapman & Hall/CRC, 1993.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{tc}} for Time Course Analysis; \code{\link{plot.STP}}, 
\code{\link{print.STP}}.
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

resSTP <- stp(Z, des)
resSTP
plot(resSTP)


\dontrun{
## Phytophthora Infestans Single Time Point Analysis (takes time...)
dataPI <- phytophthora[[4]]
desPI <- c(rep(1,8), rep(2,8))
resPI <- stp(dataPI, desPI)
resPI
plot(resPI, tp="60 hai")
}
}