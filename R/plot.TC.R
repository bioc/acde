plot.TC <-
function (x, iRatios = TRUE, FDR = TRUE, AC = TRUE, WARNINGS = FALSE, 
    ...) 
{
    ACT <- sum(apply(as.matrix(1:length(x$act)), 1, 
        function(tp) class(x$act[[tp]]) == "STP")) == length(x$act)
    GCT <- sum(apply(as.matrix(1:length(x$gct)), 1, 
        function(tp) class(x$gct[[tp]]) == "STP")) == length(x$gct)
    cols <- FDR + AC
    rows <- length(x[[1]])
    tPoints <- x$tPoints
    if (iRatios) {
        par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2.5, 1), las = 1)
        barplot(height = 100 * x$iRatios, col = "grey50", 
            main = "Inertia ratios (%)")
    }
    if (ACT) {
        if (FDR) {
            par(oma = c(0, 0, 2, 0))
            plot(x$act[[x$activeTP]], FDR, AC = FALSE, WARNINGS, 
                tp = tPoints[[x$activeTP]])
            mtext(paste("Active vs Complementary time points -", 
                tPoints[x$activeTP]), outer = TRUE, cex = 1.5, 
                font = 2)
        }
        if (AC) {
            par(mfrow = c(ceiling(rows/2), 2), mar = c(4.5, 4.5, 
                2.5, 1), oma = c(0, 0, 2, 0), las = 1)
            for (tp in 1:length(tPoints)) {
                plot(x$act[[tp]], FDR = FALSE, tp = tPoints[tp])
            }
            mtext(paste("Active vs Complementary time points -", 
                tPoints[x$activeTP]), outer = TRUE, cex = 1.5, 
                font = 2)
        }
    }
    if (GCT) {
        mf <- if (AC && FDR) 
            c(1, cols)
        else c(ceiling(rows/2), 2)
        par(mfrow = mf, mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 
            0, 2, 0), las = 1)
        for (tp in 1:length(tPoints)) {
            auxSTP <- x$gct[[tp]]
            if (auxSTP$astar > auxSTP$alpha) {
                auxSTP$dgenes <- rep("no-diff.", length(auxSTP$dgenes))
                auxSTP$tstar <- NA
            }
            plot(auxSTP, FDR, AC, WARNINGS, tPoints[tp])
            if (AC && FDR) 
                mtext(paste("Group conformation through time -", 
                    tPoints[tp]), outer = TRUE, cex = 1.5, font = 2)
        }
        if (!(AC && FDR)) 
            mtext("Group conformation through time", outer = TRUE, 
                cex = 1.5, font = 2)
    }
}
