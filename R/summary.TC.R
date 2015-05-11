summary.TC <-
function (object, ...) 
{
    dec <- 3
    ACT <- as.logical(sum(apply(as.matrix(1:length(object$act)), 
        1, function(tp) class(object$act[[tp]]) == "STP")))
    GCT <- sum(apply(as.matrix(1:length(object$gct)), 1, 
        function(tp) class(object$gct[[tp]]) == 
        "STP")) == length(object$gct)
    cat("\n")
    cat("Time course analysis for detecting differentially\n")
    cat("expressed genes in microarray data.\n \n")
    cat("Inertia ratios (%): \n")
    print(round(100 * object$iRatios, 2))
    if (ACT) {
        cat("\n\nActive vs complementary time points analysis:\n\n")
        cat(paste("Active timepoint:", names(object$act)[object$activeTP], 
            "\n\n"))
        cat(paste("Achieved FDR:", 
            round(100*object$act[[object$activeTP]]$astar, 2), "%.\n\n"))
        cat("Differentially expressed genes:")
        print(table(object$act[[object$activeTP]]$dgenes))
    }
    if (GCT) {
        cat("\n\nGroups conformation through time analysis:\n\n")
        tpstar <- rep(FALSE, length(object$gct))
        for (tp in 1:length(object$gct)) {
            tpstar[tp] <- object$gct[[tp]]$astar <= object$gct[[tp]]$alpha
        }
        tpstar <- (1:length(object$gct))[tpstar]
        if (length(tpstar) == 1) {
            cat("Differentially expressed genes:")
            print(table(object$act[[tpstar]]$dgenes))
        }
        else {
            cat("Differentially expressed genes:")
            for (i in 1:(length(tpstar) - 1)) {
                for (j in (i + 1):length(tpstar)) {
                    cat(paste("\n", names(object$gct)[tpstar[i]], 
                        "vs", paste(names(object$gct)[tpstar[j]])))
                    print(addmargins(table(object$gct[[tpstar[i]]]$dgenes, 
                        object$gct[[tpstar[j]]]$dgenes)))
                }
            }
        }
    }
}
