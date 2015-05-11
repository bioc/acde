\name{plot.TC}
\alias{plot.TC}
\title{
Plot Method for Time Course Analysis
}
\description{
a method for the plot generic. It is designed for displaying plots of 
the estimated FDR and the genes' classification when performing a Time 
Course Analysis for detecting differentially expressed genes in gene 
expression data.
}
\usage{
\method{plot}{TC}(x, iRatios=TRUE, FDR = TRUE, AC = TRUE, 
    WARNINGS = FALSE, ...)
}
\arguments{
\item{x}{
if \code{TRUE}, a plot of inertia ratios for all the time points 
is displayed.  
}
\item{iRatios}{
an object of class '\code{TC}' as returned by function \code{tc}.  
}
\item{FDR}{
if \code{TRUE}, a plot of the estimated FDRs are displayed for each 
time point.  
}
\item{AC}{
if \code{TRUE}, a plot of the differentially expressed genes in the 
artificial components is displayed for each time point.  
}
\item{WARNINGS}{
if \code{TRUE} and if a BCa confidence upper bound was computed for 
obtaining \code{x}, the threshold values for which an extreme order 
statistic was used in the BCa computations are shown (these warnings 
are produced in calls to \code{\link{boot.ci}}).  
}
\item{\dots}{
further arguments passed to or from other methods.
}
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\seealso{
\code{\link{tc}}, \code{\link{print.TC}}, \code{\link{summary.TC}}.
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
