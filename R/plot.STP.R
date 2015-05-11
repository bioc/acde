plot.STP <-
function (x, FDR = TRUE, AC = TRUE, WARNINGS = FALSE, tp = NULL, 
    ...) 
{
    if (FDR && length(x$Q) > 1) {
        BCa <- !is.character(x$bca)
        par(mar = c(4.5, 3, 2.5, 1), las = 1)
        tit <- if (is.null(tp)) 
            "Single Time Point Analysis: FDR"
        else paste("Single Time Point A.: FDR -", tp)
        plot(x$th, x$Q, type = "n", xlab = "t", ylab = "", main = tit, 
            ylim = c(0, 1), xlim = c(0, max(x$th)))
        if (BCa) {
            noNA <- !is.na(x$bca$cbound)
            polygon(x = c(min(x$th[noNA]), x$th[noNA], max(x$th[noNA])), 
                y = c(0, x$bca$cbound[noNA], 0), col = "grey92", 
                border = NA)
        }
        abline(v = 0, lty = 1, col = "black")
        abline(h = 0, lty = 1, col = "black")
        lines(x$th, x$Q, col = "blue")
        if (BCa) 
            lines(x$th, x$bca$cbound, col = "green")
        abline(v = x$tstar, lty = 2, col = "black")
        if (!BCa) {
            legend("topright", legend = c(expression(hat(Q)[0.5](italic(t))), 
                paste("t* (", round(x$tstar, 2), ")", sep = "")), 
                col = c("blue", "black"), lty = c(1, 2), cex = 1, 
                bty = "n")
        }
        else if (!WARNINGS) {
            legend("topright", legend = c(expression(hat(Q)[0.5](italic(t))), 
                paste("Upper bound (", round(100 * x$gamma, 1), 
                    "%)", sep = ""), paste("t* (", round(x$tstar, 
                    2), ")", sep = "")), col = c("blue", "Green", 
                "black"), lty = c(1, 1, 2), cex = 1, bty = "n")
        }
        else {
            war <- x$bca$warnings != "Ok"
            points(x$th[war], x$bca$cbound[war], pch = 1, col = "red", 
                cex = 0.7)
            boundNA <- is.na(x$bca$cbound)
            points(x$th[boundNA], rep(0, sum(boundNA)), pch = 1, 
                col = "red", cex = 0.7)
            legend("topright", legend = c(expression(hat(Q)[0.5](italic(t))), 
                paste("Upper bound (", round(100 * x$gamma, 1), 
                    "%)", sep = ""), paste("t* (", round(x$tstar, 
                    2), ")", sep = ""), "Warnings generated", 
                    "in BCa computations."), 
                col = c("blue", "green", "black", "red"), lty = c(1, 
                    1, 2, 0, 0), pch = c(NA_integer_, NA_integer_, 
                    NA_integer_, 1, NA_integer_), cex = 1)
        }
    }
    if (AC) {
        up <- x$dgenes == "up-reg."
        down <- x$dgenes == "down-reg."
        nd <- x$dgenes == "no-diff."
        par(mar = c(4.5, 4.5, 2.5, 1), las = 1)
        tit <- if (is.null(tp)) 
            "Single Time Point Analysis: AC"
        else paste("Single Time Point A.: AC -", tp)
        plot(x$ac[, 1], x$ac[, 2], type = "n", 
            xlab = expression(psi[1]), 
            ylab = expression(psi[2]), main = tit)
        abline(v = 0, lty = 1, col = "black")
        abline(h = 0, lty = 1, col = "black")
        abline(h = c(-x$tstar, x$tstar), lty = 2)
        legend("topleft", legend = c(paste(" Up-regulated (", 
            sum(up), ")", sep = ""), paste(" Down-regulated (", 
            sum(down), ")", sep = ""), paste(" No diff. expr. (", 
            sum(nd), ")", sep = "")), pch = c(21, 21, 20), 
            col = c("darkgreen", "darkred", "grey60"), 
            pt.bg = c("green", "red"), bty = 1, bg = "white")
        points(x$ac[nd, 1], x$ac[nd, 2], pch = 20, col = "grey60")
        points(x$ac[up, 1], x$ac[up, 2], pch = 21, bg = "green", 
            col = "darkgreen")
        points(x$ac[down, 1], x$ac[down, 2], pch = 21, bg = "red", 
            col = "darkred")
    }
}
