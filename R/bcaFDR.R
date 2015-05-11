bcaFDR <-
function (Z, design, th = NULL,B = 100, lambda = 0.5, PER = FALSE, 
    R = 1000, gamma = 0.95, Q = NULL, ...) 
{
    p1 <- table(design)[1]
    p2 <- table(design)[2]
    R <- min(R, choose(2 * p1 - 1, p1) * choose(2 * p2 - 1, p2))
    fdrBoot <- boot(t(Z), function(trZ, ii) fdr(t(trZ[ii, ]), 
        design, th, B, lambda, PER)$Q, R, stype = "i", strata = design, 
        ...)
    if (is.null(th) || is.null(Q)) {
        resFDR <- fdr(Z, design, th, B, lambda, PER, ...)
        th <- if (is.null(th)) resFDR$th else th
        Q <- if (is.null(Q)) resFDR$Q else Q
    }
    fdrBoot$t0 <- Q
    W <- rep("Ok", length(th))
    tryCatchWE <- function(expr, i) {
        wHandler <- function(w) {
            assign("W[i]", w$message, envir = 
                parent.env(environment()))
            invokeRestart("muffleWarning")
        }
        withCallingHandlers(tryCatch(expr, error = function(e) {
            assign("W[i]", e$message, envir = 
                parent.env(environment()))
            return(NA)
        }), warning = wHandler)
    }
    funBCa <- function(i) {
        tryCatchWE(expr = boot.ci(fdrBoot, conf = (1 - 2 * (1 - 
            gamma)), type = "bca", index = i)$bca[5], i)
    }
    bca <- apply(as.matrix(1:length(th)), 1, funBCa)
    list(cbound = bca, warnings = W)
}
