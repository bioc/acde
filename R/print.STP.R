print.STP <-
function (x, headerSTP = TRUE, ...) 
{
    dec <- 3
    nde <- 1 * (x$dgenes == "up-reg.") + 2 * (x$dgenes == "down-reg.") + 
        3 * (x$dgenes == "no-diff.")
    o <- order(nde, x$qvalues)
    q <- if (is.character(x$qvalues)) 
        x$qvalues
    else round(x$qvalues[o], 3)
    df <- data.frame(round(x$ac[o, 1], dec), round(x$ac[o, 2], 
        3), q, x$dgenes[o])
    colnames(df) <- c("psi1", "psi2", "Q-value", "Diff. expr.")
    if (headerSTP) {
        cat("\nSingle time point analysis for detecting differentially\n")
        cat("expressed genes in microarray data.\n")
    }
    cat(paste("\nAchieved FDR: ", round(100 * x$astar, 1), "%.\n", 
        sep = ""))
    if (!is.character(x$bca)) {
        bcastar <- min(x$bca$cbound[x$th == x$tstar])
        cat(paste(round(x$gamma * 100, 1), "% BCa upper bound for the FDR: ", 
            round(100 * bcastar, 1), "%.\n", sep = ""))
        cat(paste("Warnings in BCa computation: ", (length(x$bca$warnings) - 
            sum(x$bca$warnings == "Ok")), ".\n \n", sep = ""))
    }
    else cat("\n")
    cat(paste("Inertia ratio: ", round(100 * x$iratio, 2), "%.\n", 
        sep = ""))
    cat(paste("tstar: ", round(x$tstar, dec), ", pi0: ", round(x$pi0, 
        dec), ", B: ", x$B, ".\n \n", sep = ""))
    cat("Differentially expressed genes:")
    print(table(x$dgenes))
    cat("\nResults: \n")
    rows <- min(sum(x$dgenes == "up-reg.") + sum(x$dgenes == 
        "down-reg.") + 10, nrow(df) * (x$astar <= x$alpha) + 
        10 * (x$astar > x$alpha))
    print(df[1:rows, ])
    cat("...\n")
    cat("\n*More results are available in the objects:\n")
    cat("$ac, $qvalues and $dgenes.\n")
}
