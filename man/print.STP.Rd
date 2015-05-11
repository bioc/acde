\name{print.STP}
\alias{print.STP}
\title{
Print Method for Single Time Point Analysis
}
\description{
a method for the print generic. It prints relevant results when 
performing a Single Time Point Analysis for detecting differentially 
expressed genes in gene expression data.
}
\usage{
\method{print}{STP}(x, headerSTP = TRUE, ...)
}
\arguments{
\item{x}{
an object of class '\code{STP}' as returned by function \code{stp}.
}
\item{headerSTP}{
if \code{FALSE}, the header is omitted (used for 
\code{\link{print.TC}}).
}
\item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
If the desired FDR level was achieved (i.e. \code{x$astar <= x$alpha}), 
the results are printed for the differentially expressed genes and 10 
more rows only. If the desired FDR level was not achieved, only ten 
rows are displayed.
}
\seealso{
\code{\link{stp}}, \code{\link{plot.STP}}.
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
