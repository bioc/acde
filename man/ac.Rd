\name{ac}
\alias{ac}
\alias{ac2}
\title{
Artificial Components for Gene Expression Data
}
\description{
Computes the artificial components for gene expression data between two 
conditions for a single time point.  
}
\usage{
ac(Z, design)

ac2(Z, design)
}
\arguments{
\item{Z }{
a numeric matrix or data.frame with \eqn{n} rows and \eqn{p} columns 
representing genes' expression levels. The rows of \eqn{Z} correspond 
to the genes in the experiment, and the columns correspond to the 
replicates. Treatment replicates are to the left, control replicates 
to the right.
}
\item{design }{
a vector of length \eqn{p} with 1's for the treatment replicates and 
2's for the control replicates 
\eqn{(1, \ldots, 1, 2, \ldots, 2)}{(1, \ldots, 1, 2, \ldots, 2)}.
}
}
\details{
This function computes the artificial components of Z, based on the 
specified design vector.  First, the function scales \eqn{Z} so that 
its columns have zero mean and unit variance.  Then computation of 
the artificial components \eqn{\mathbf{\psi}_1}{\psi[1]} and 
\eqn{\mathbf{\psi}_2}{\psi[2]} is performed as \eqn{\psi_1 = 
Z\mathbf{v_1}}{\psi[1] = Zv[1]}, where \eqn{\mathbf{v_1} = 
\left(1, \ldots, 1\right) / \sqrt{p}}{v[1] = (1, \ldots , 1) / sqrt(p)}, 
and \eqn{\psi_2 = Z\mathbf{v_2}}{\psi[2] = Zv[2]}, where \eqn{\mathbf{v_2} = 
\left(1, \ldots, 1, -1, \ldots, -1 \right) / \sqrt{pp_1(p-p_1)}}{v[2] = (1, 
\ldots , 1, -1, \ldots , -1 ) / sqrt(p*p[1]*(p-p[1]))}.  Here, 
\eqn{p_1}{p[1]} is the number of treatment replicates, and 
\eqn{\mathbf{v_2}}{v[2]} has \eqn{p_1}{p[1]} positive and 
\eqn{p-p_1}{p-p[1]} negative entries.   
}
\value{
\code{ac} returns a matrix with the artificial components 
\eqn{\psi_1}{\psi[1]} and \eqn{\psi_2}{\psi[2]} in the columns. 

\code{ac2} returns a matrix with the second artificial component 
\eqn{\psi_2}{\psi[2]} in the only column.
}
\references{
Acosta, J. P. (2015) \emph{Strategy for Multivariate Identification 
of Differentially Expressed Genes in Microarray Data}. Unpublished 
MS thesis. Universidad Nacional de Colombia, Bogot\'a.
}
\author{
Juan Pablo Acosta (\email{jpacostar@unal.edu.co}).
}
\examples{
## Computes the artificial components for the 
## phitophthora infestans data at 60 hai.
psi <- ac(phytophthora[[4]], c(rep(1,8), rep(2,8)))
plot(x=psi[,1], y=psi[,2])
}