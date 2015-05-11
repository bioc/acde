stp <-
function (Z, design, alpha = 0.05, B = 100, lambda = 0.5, th = NULL, 
    PER = FALSE, BCa = FALSE, gamma = 0.95, R = 1000, ...) 
{
    Z <- as.matrix(Z)
    if (ncol(Z) != length(design)) 
        stop("Replicates in Z do not correspond to classes in design.")
    if (length(table(design)) != 2) 
        stop(paste("There are more than 2 classes in design.",
            "Only two class comparison is available."))
    if (BCa && min(table(design)) <= 1) 
        stop(paste("There must be at least two replicates",
            "in each condition for BCa computations"))
    if (BCa && R <= 1) 
        stop(paste("R must be greater than one and",
            "around 1000 for good results."))    
    if (B <= 1) 
        stop(paste("B must be greater than one and",
            "around 100 for good results."))
    if (length(lambda) != 1) 
        stop("lambda should be a scalar.")
    if (lambda <= 0 || lambda >= 1) 
        stop("lambda must be in (0,1).")
    if (length(alpha) != 1) 
        stop("alpha should be a scalar.")
    if (alpha <= 0 || alpha >= 1) 
        stop("alpha must be in (0,1).")
    if (length(gamma) != 1) 
        stop("gamma should be a scalar.")
    if (gamma <= 0 || gamma >= 1) 
        stop("gamma must be in (0,1).")
    if (!is.null(th)) {
        if (sum(th <= 0) > 0) {
            stop("Threshold values in th must be positive.")
        } else if (max(th) > max(abs(ac2(Z, design)))) {
            stop(paste("Threshold values must be less than",
                "max(abs(ac2(Z, design)))."))
        } else th <- th[order(th)]
    }
    QVal <- is.null(th)
    n <- nrow(Z)
    gNames <- if (is.null(rownames(Z))) {paste("g", 1:nrow(Z))
    } else rownames(Z)
    ac <- ac(Z, design)
    resFDR <- fdr(Z, design, th, B, lambda, PER, ...)
    Q <- resFDR$Q
    th <- resFDR$th
    psi2 <- ac[, 2]
    pi0 <- resFDR$pi0
    iRatio <- var(psi2)/eigen(cor(Z))$values[1]
    tstar <- if (sum(Q <= alpha) > 0) 
        min(th[Q <= alpha])
    else th[order(Q)[1]]
    astar <- if (sum(Q <= alpha) > 0) 
        Q[th == tstar]
    else min(Q)
    dgenes <- rep("no-diff.", n)
    if (sum(psi2 <= -tstar) > 0) 
        dgenes[psi2 <= -tstar] <- "down-reg."
    if (sum(psi2 >= tstar) > 0) 
        dgenes[psi2 >= tstar] <- "up-reg."
    dgenes <- as.factor(dgenes)
    qvalues <- if (QVal) 
        qval(Q, psi2)
    else rep("Not computed", n)
    bca <- if (BCa) 
        bcaFDR(Z, design, th, B, lambda, PER, R, gamma, Q, ...)
    else "Not computed"
    res <- list(dgenes = dgenes, tstar = tstar, astar = astar, 
        Q = Q, th = th, qvalues = qvalues, pi0 = pi0, B = B, 
        lambda = lambda, ac = ac, gNames = gNames, iRatio = iRatio, 
        bca = bca, gamma = gamma, alpha = alpha, call = match.call())
    class(res) <- "STP"
    res
}
