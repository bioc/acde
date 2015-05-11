qval <-
function (Q, psi2) 
{
    n <- length(Q)
    q <- rep("NC", n)
    o <- order(abs(psi2))
    q <- rep(1, n)
    q[o[1]] <- Q[1]
    for (i in 2:n) {
        q[o[i]] <- min(Q[i], q[o[i - 1]])
    }
    q
}
