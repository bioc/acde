\name{tc}
\alias{tc}
\title{
Time Course Analysis for Detecting Differentially Expressed Genes
}
\description{
Performs the Time Course Analysis from Acosta (2015) for detecting 
differentially expressed genes in time course experiments for gene 
expression data.  
}
\usage{
tc(data, designs, tPoints = NULL, 
    method = c("active vs complementary", "groups conformation"), 
    activeTP = NULL, alpha = 0.05, B = 100, lambda = 0.5, 
    PER = FALSE, BCa = FALSE, gamma = 0.95, R = 1000, ...)
}
\arguments{
\item{data}{
a list with matrices or data.frames representing genes' expression levels 
at each time point. The rows of each matrix correspond to the genes in 
the experiment, and the columns correspond to the replicates. Treatment 
replicates are to the left, control replicates to the right.
}
\item{designs}{
a list with a vector for each time point, of length equal to the number 
of columns in the respective matrix or data.frame in \code{data}, with 
1's for the treatment replicates and 2's for the control replicates.
}
\item{tPoints}{
a character vector with the names of the timepoints.
}
\item{method}{
if "active vs complementary", an analysis following the \emph{active vs 
complementary time points} approach (Acosta, 2015) is performed. If "groups 
conformation", an analysis following the \emph{groups conformation through 
time} approach (Acosta, 2015) is performed. The default is both.
}
\item{activeTP}{
numeric. The index of the active timepoint in \code{tPoints} for the 
\emph{active vs complementary time points} approach. If \code{NULL} 
(default), the active time point is selected via inertia ratios. 
}
\item{alpha}{
between 0 and 1. Desired level for controlling the false discovery 
rate (FDR).
}
\item{B}{
Number of bootstrap or permutation replications for estimating the FDR.
}
\item{lambda}{
Parameter for the estimation of \eqn{\pi_0}{pi0} and, further, the 
estimation of the FDR (see Storey, 2002).
}
\item{PER}{
If \code{FALSE} (default), bootstrap replications are used to estimate 
the FDR. If \code{TRUE}, permutation replications are used instead. 
}
\item{BCa}{
If TRUE, a BCa confidence upper bound for the FDR is computed 
(see Efron and Tibshirani, 1994).
}
\item{gamma}{
Confidence level for the FDR's BCa upper confidence bound.
}
\item{R}{
Number of bootstrap replications for the computation of the FDR's 
BCa upper confidence bound.
}
\item{\dots}{
additional arguments for parallel computation in \code{boot} 
function (see Details).
}
}
\details{
In the \emph{active vs complementary time points} approach, the 
time point that maximizes the inertia ratio is selected as the 
\emph{active} time point. Then, a Single Time Point Analysis (\code{stp}) 
is performed on this time point and plots of the behavior throughout the 
time course of the differentially expressed genes identified in this time 
point are displayed.

In the \emph{groups conformation through time} approach, a Single Time 
Point Analysis (\code{stp}) is performed at each time point and plots 
are displayed showing the behaviour of the differential expression process 
throughout the time course; that is, how many genes are differentially 
expressed and how strong is the differential expression at each time point.

For details on the computations performed in this function, see Acosta 
(2015).

Additional parameters in the '\code{...}' argument are used for parallel 
computation in bootstrap calculations. These are supplied to calls to the 
\code{boot} function in package \code{boot}. With this in mind, the use 
of additional arguments must be restricted to arguments \code{parallel} 
and \code{ncpus} from function \code{boot}.
}
\section{Warning }{
If argument \code{BCa=TRUE}, computations may take a considerable 
amount of time.
}
\value{
\code{tc} returns an object of class '\code{TC}', which is a 
list with components:

\item{iRatios }{
inertia ratios for each time point.
}
\item{gct }{
results for the \emph{groups conformation through time} approach. 
A list with an object of class '\code{STP}' for each timepoint.
}
\item{act }{
results for the \emph{active vs complementary time points} approach. 
A list with an object of class '\code{STP}' for each timepoint.
}
\item{activeTP }{
the index of the active timepoint in \code{tPoints} used in the 
\emph{active vs complementary time points} approach.
}
\item{tPoints}{
a character vector with the names of the timepoints.
} 
\item{call }{The matched call.}
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification of 
Differentially Expressed Genes in Microarray Data}. Unpublished MS 
thesis. Universidad Nacional de Colombia, Bogot\'a.

Storey, J. D. (2002) \emph{A direct approach to false discovery rates.} 
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 
\bold{64(3): }479--498.

Efron B. and Tibshirani R. J. (1994) \emph{An Introduction to the 
Bootstrap.} Chapman & Hall/CRC, 1993.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{stp}} for Single Time Point Analysis; 
\code{\link{plot.TC}}, \code{\link{print.TC}}, \code{\link{summary.TC}}.
}
\examples{
## Time course analysis for 500 genes with 10 treatment 
## replicates and 10 control replicates
tPts <- c("h0", "12h", "24h")
n <- 500; p <- 20; p1 <- 10
Z <- vector("list", 3)
des <- vector("list", 3)
for(tp in 1:3){ des[[tp]] <- c(rep(1, p1), rep(2, (p-p1))) }
mu <- as.matrix(rexp(n, rate=1))
### h0 time point (no diff. expr.)
Z[[1]] <- t(apply(mu, 1, function(mui) rnorm(p, mean=mui, sd=1)))
### h12 time point (diff. expr. begins)
Z[[2]] <- t(apply(mu, 1, function(mui) rnorm(p, mean=mui, sd=1)))
#### Up regulated genes
Z[[2]][1:5,1:p1] <- Z[[2]][1:5,1:p1] + 
    matrix(runif(5*p1, 1, 3), nrow=5)
#### Down regulated genes
Z[[2]][6:15,(p1+1):p] <- Z[[2]][6:15,(p1+1):p] + 
    matrix(runif(10*(p-p1), 1, 2), nrow=10)
### h24 time point (maximum differential expression)
Z[[3]] <- t(apply(mu, 1, function(mui) rnorm(p, mean=mui, sd=1)))
#### 5 up regulated genes
Z[[3]][1:5,1:p1] <- Z[[3]][1:5,1:p1] + 5
#### 10 down regulated genes
Z[[3]][6:15,(p1+1):p] <- Z[[3]][6:15,(p1+1):p] + 4

resTC <- tc(Z, des)
resTC
summary(resTC)
plot(resTC)

\dontrun{
## Phytophthora Infestans Time Course Analysis (takes time...)
dataPI <- phytophthora
desPI <- vector("list", 4)
for(tp in 1:4){ desPI[[tp]] <- c(rep(1, 8), rep(2, 8)) }
resPI <- tc(dataPI, desPI)
resPI
summary(resPI)
plot(resPI)
}
}
