print.TC <-
function (x, ...) 
{
    dec <- 3
    ACT <- as.logical(sum(apply(as.matrix(1:length(x$act)), 1, 
        function(tp) class(x$act[[tp]]) == "STP")))
    GCT <- sum(apply(as.matrix(1:length(x$gct)), 1, 
        function(tp) class(x$gct[[tp]]) == "STP")) == length(x$gct)
    cat("\n")
    cat("Time course analysis for detecting differentially expressed\n")
    cat("genes in microarray data.\n \n")
    cat("Inertia ratios (%): \n")
    print(round(100 * x$iRatios, 2))
    cat("\n")
    if (ACT) 
        cat(paste("Active timepoint:", names(x$act)[x$activeTP], 
            "\n"))
    if (GCT) {
        for (tp in 1:length(x$gct)) {
            cat(paste("\n\nTIME POINT: ", names(x$gct)[tp], ".\n", sep=""))
            print(x$gct[[tp]], headerSTP = FALSE)
        }
    }
    else print(x$act[[x$activeTP]])
}
