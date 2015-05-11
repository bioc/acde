ac2 <-
function (Z, design) 
{
    if (!(is.matrix(Z) || is.data.frame(Z))) 
        stop("Z must be a matrix or a data.frame")
    if (ncol(Z) < 2) 
        stop("Z must have at least 2 columns / variables.")
    if (ncol(Z) != length(design)) 
        stop("Replicates in Z do not correspond to classes in design.")
    if (length(table(design)) != 2) 
        stop("There must be two conditions in design.")
    if (sum(diag(var(Z) == 0)) >= 1) 
        stop("There are columns of Z with cero variance.")
    X <- if (max(abs(colMeans(Z))) > 0 || sum(diag(var(Z))) != 
        ncol(Z)) 
        scale(Z)
    else Z
    colnames(X) <- colnames(Z)
    X <- X[, order(design)]
    p <- length(design)
    p1 <- table(design)[1]
    p2 <- table(design)[2]
    v2 <- c(rep(p2, p1), rep(-p1, p2))/sqrt(p * p1 * p2)
    psi2 <- X %*% v2
    psi2
}
