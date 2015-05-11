\name{plot.STP}
\alias{plot.STP}
\title{
Plot Method for Single Time Point Analysis
}
\description{
a method for the plot generic. It is designed for displaying plots of the 
estimated FDR and the genes' classification when performing a Single Time 
Point Analysis for detecting differentially expressed genes in gene 
expression data.
}
\usage{
\method{plot}{STP}(x, FDR=TRUE, AC=TRUE, WARNINGS=FALSE, tp=NULL, ...)
}
\arguments{
\item{x}{
an object of class '\code{STP}' as returned by function \code{stp}.  
}
\item{FDR}{
if \code{TRUE}, a plot of the estimated FDR is displayed.  
}
\item{AC}{
if \code{TRUE}, a plot of the differentially expressed genes in the 
artificial components is displayed.  
}
\item{WARNINGS}{
if \code{TRUE} and if a BCa confidence upper bound was computed for 
obtaining \code{x}, the threshold values for which an extreme order 
statistic was used in the BCa computations are shown (these warnings 
are produced in calls to \code{\link{boot.ci}}).  
}
\item{tp}{
a character string to be added at the end of the plot's title 
(used for adding time points in \code{\link{plot.TC}}).
}
\item{\dots}{
further arguments passed to or from other methods.
}
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{stp}}, \code{\link{print.STP}}.
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
}
