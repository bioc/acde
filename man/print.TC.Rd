\name{print.TC}
\alias{print.TC}
\alias{summary.TC}
\title{
Print Method for Time Course Analysis
}
\description{
methods for the print and summary generics that print relevant 
results when performing a Time Course Analysis for detecting 
differentially expressed genes in gene expression data.
}
\usage{
\method{print}{TC}(x, ...)

\method{summary}{TC}(object, ...)
}
\arguments{
\item{x}{
an object of class '\code{TC}' as returned by function \code{tc}.
}
\item{object}{
an object of class '\code{TC}' as returned by function \code{tc}.
}
\item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
With \code{print}, at each time point, if the desired FDR level 
was achieved (i.e. \code{x$astar <= x$alpha}), the results are 
printed for the differentially expressed genes and 10 more rows 
only. If the desired FDR level was not achieved, only ten rows 
are displayed.

\code{summary} prints a more concise version of the results.
}
\seealso{
\code{\link{tc}}, \code{\link{plot.TC}}.
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
