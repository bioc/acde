fdr <-
function (Z, design, th = NULL, B = 100, lambda = 0.5, PER = FALSE, ...) 
{
    if (is.null(th)) {
        th <- abs(ac2(Z, design))
        th <- th[order(th)]
    }
    psi2 <- ac2(Z, design)
    n <- length(psi2)
    stat <- function(Ztr, ii) {
        s <- abs(ac2(t(Ztr[ii, ]), design))
        c(s, apply(as.matrix(th), 1, function(thi) sum(s >= thi)))
    }
    parSim <- if (PER) 
        "permutation"
    else "ordinary"
    aux <- boot(t(Z), stat, B, stype = "i", sim = parSim, ...)$t
    psi2Boot <- aux[, 1:n]
    rBoot <- as.matrix(aux[, -(1:n)])
    ER <- colMeans(rBoot)
    tLambda <- quantile(psi2Boot, lambda)
    wLambda <- sum(psi2 < tLambda)
    pi0 <- wLambda/(n * (1 - lambda))
    pi0 <- pi0 * (pi0 < 1) + (pi0 > 1)
    r <- apply(as.matrix(th), 1, function(thi) sum(abs(psi2) >= 
        thi))
    Q <- pi0 * (ER/r)
    Q <- Q * (Q < 1) + (Q >= 1)
    res <- list(Q = Q, th = th, pi0 = pi0, B = B, lambda = lambda, 
        call = match.call())
    res
}
