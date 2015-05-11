\name{phytophthora}
\alias{phytophthora}
\docType{data}
\title{
Gene Expression Data for Tomato Plants Inoculated with 
\emph{Phytophthora infestans}
}
\description{
Gene expression data for 16 tomato plants (line IL6-2) in field conditions. 
8 of these plants were inoculated with \emph{Phytophthora infestans}, and 
the other 8 were mock-inoculated with sterile water. Leaf tissue samples 
from each replicate were taken at 12 hours before and 12, 36 and 60 hours 
after inoculation. We refer to 12 hours before inoculation as the h0 time 
point. Expression levels were obtained for 13440 genes.}
\usage{data("phytophthora")}
\format{
A list with four matrices representing expression levels for 13440 genes 
(rows) in 16 replicates (columns) at time points "h0", "h12", "h36" and 
"h60". At each time point, the first 8 columns correspond to treatment 
(inoculated) replicates and the last 8 columns correspond to control 
(mock-inoculated) replicates. The names of the genes are the names of the 
rows in each matrix.
}
\details{
For details about experimental conditions, see Restrepo et al. (2005) 
and Cai et al. (2013).
}
\source{
Tomato Expression Database website (\url{http://ted.bti.cornell.edu/}), 
experiment E022 (Restrepo et al., 2005).
}
\references{
Restrepo, S., Cai, G., Fry, W. E. and Smart, C. D. (2005) 
\emph{Gene expression profiling of infection of tomato by Phytophthora 
infestans in the field}. Phytopathology, \bold{95(S88)}.

Cai, G., Restrepo, S., Myers, K., Zuluaga, P., Danies, G., Smart, C. and 
Fry W.E. (2013) \emph{Gene profiling in partially resistant and susceptible 
near-isogenic tomatoes in response to late blight in the field}. Molecular 
plant pathology, \bold{14(2): }171--184.
}
\examples{
for(tp in 1:4){
    cat(paste("Time Point:", names(phytophthora[tp]), "\n"))
    print(phytophthora[[tp]][1:10,])
    cat("...\n \n")
}
}
\keyword{datasets}