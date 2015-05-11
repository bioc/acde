\name{acde-package}
\alias{acde-package}
\alias{acde}
\docType{package}
\title{
Artificial Components Detection of Differentially Expressed Genes
}
\description{
This package provides a multivariate inferential analysis method for 
detecting differentially expressed genes in gene expression data. It 
uses artificial components, close to the data's principal components 
but with an exact interpretation in terms of differential genetic 
expression, to identify differentially expressed genes while 
controlling the false discovery rate (FDR). The methods on this 
package are described in the article \emph{Multivariate Method for 
Inferential Identification of Differentially Expressed Genes in Gene 
Expression Experiments} by Acosta (2015).
}
\details{
\tabular{ll}{
Package: \tab acde\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-25\cr
License: \tab GLP-3\cr
LazyData: \tab yes\cr
Depends: \tab R(>= 3.1), ade4(>= 1.6), boot(>= 1.3)\cr
Encoding: \tab UTF-8\cr
Built: \tab R 3.1.2; 2015-05-01; unix\cr
}

Index:
\preformatted{
ac              Artificial Components for Gene 
                Expression Data
acde-package    Artificial Components Detection of
                Differentially Expressed Genes
bcaFDR          BCa Confidence Upper Bound for the FDR.
fdr             False Discovery Rate Computation
phytophthora    Gene Expression Data for Tomato Plants
                Inoculated with _Phytophthora infestans_
plot.STP        Plot Method for Single Time Point Analysis
plot.TC         Plot Method for Time Course Analysis
print.STP       Print Method for Single Time Point Analysis
print.TC        Print Method for Time Course Analysis
qval            Q-Values Computation
stp             Single Time Point Analysis for Detecting
                Differentially Expressed Genes
tc              Time Course Analysis for Detecting
                Differentially Expressed Genes
}
}
\author{
Juan Pablo Acosta, Liliana Lopez-Kleine

Maintainer: Juan Pablo Acosta <jpacostar@unal.edu.co>
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification of 
Differentially Expressed Genes in Microarray Data}. Unpublished MS 
thesis. Universidad Nacional de Colombia, Bogot\'a.
}
\keyword{ package }
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
}
