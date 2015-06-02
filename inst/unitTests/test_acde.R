# Tests for function ac
test_ac <- function ()
{
    # Correct cases
    ## 1
    Z1 <- matrix(c(1,0,0,1), ncol=2)
    des1 <- c(1,2)
    checkEqualsNumeric(ac(Z1, des1), matrix(c(0,0,1,-1), nrow=2))
    ## 2
    Z2 <- cbind(Z1, Z1)
    des2 <- c(des1, des1)
    checkEqualsNumeric(ac(Z2, des2), matrix(c(0,0,sqrt(2),-sqrt(2)), nrow=2))
    checkEqualsNumeric(ac(Z2, c(1,1,2,2)), matrix(0,2,2))
    
    # Errors
    checkException(ac(), silent=TRUE)
    checkException(ac(des2,des1), silent=TRUE)
    checkException(ac(Z2,des1), silent=TRUE)
    checkException(ac(Z1, c(1,1)), silent=TRUE)
    checkException(ac(Z2, c(1,1,2,3)), silent=TRUE)
    checkException(ac(matrix(0,2,2), des1), silent=TRUE)
}
# Tests for function ac2
test_ac2 <- function ()
{
    # Correct cases
    ## 1
    Z1 <- matrix(c(1,0,0,1), ncol=2)
    des1 <- c(1,2)
    checkEqualsNumeric(ac2(Z1, des1), matrix(c(1,-1), nrow=2))
    ## 2
    Z2 <- cbind(Z1, Z1)
    des2 <- c(des1, des1)
    checkEqualsNumeric(ac2(Z2, des2), matrix(c(sqrt(2),-sqrt(2)), nrow=2))
    checkEqualsNumeric(ac2(Z2, c(1,1,2,2)), matrix(0,2,1))
    
    # Errors
    checkException(ac2(), silent=TRUE)
    checkException(ac2(des2,des1), silent=TRUE)
    checkException(ac2(Z2,des1), silent=TRUE)
    checkException(ac2(Z1, c(1,1)), silent=TRUE)
    checkException(ac2(Z2, c(1,1,2,3)), silent=TRUE)
    checkException(ac2(matrix(0,2,2), des1), silent=TRUE)
}
# Tests for function stp
test_stp <- function ()
{
    # Correct cases
    ## First gene up-reg., second down-reg.
    Z <- matrix(runif(100), 10, 20); Z[1,1:10] <- 10; Z[2,11:20] <- 10
    des <- c(rep(1,10), rep(2,10))
    res <- stp(Z, des)
    checkEqualsNumeric(table(res$dgenes), c(1,8,1))
    checkTrue(res$astar <= 0.05)
    ## One value for th
    checkEqualsNumeric(table(stp(Z, des, th=2)$dgenes), c(1,8,1))
    checkException(table(stp(Z, des, th=10)$dgenes), silent=TRUE)

    # Errors
    ## Missing arguments
    checkException(stp(), silent=TRUE)
    checkException(stp(matrix(0,10,4)), silent=TRUE)
    checkException(stp(design=c(1,1,2,2)), silent=TRUE)
    ## Main arguments do not match
    checkException(stp(matrix(0,10,4), c(1,2,2)), silent=TRUE)
    checkException(stp(matrix(0,10,4), c(1,2,3,4)), silent=TRUE)
    checkException(stp(matrix(rnorm(40),10,4), c(1,1,2,2), th=(-10:10)), 
        silent=TRUE)
    checkException(stp(matrix(rnorm(30),10,3), c(1,2,2), BCa=TRUE), 
        silent=TRUE)    
    ## Main arguments do match, but no variability in Z
    checkException(stp(matrix(0,10,4), c(1,1,2,2)), silent=TRUE)
    ## Secondary arguments do not match
    Z <- matrix(rnorm(40),10,4)
    des <- c(1,1,2,2)
    checkException(stp(Z, des, B=1), silent=TRUE)
    checkException(stp(Z, des, lambda=0), silent=TRUE)
    checkException(stp(Z, des, lambda=c(-1,0.5,2)), silent=TRUE)
    checkException(stp(Z, des, lambda=-1), silent=TRUE)
    checkException(stp(Z, des, alpha=c(-1,0.5,2)), silent=TRUE)
    checkException(stp(Z, des, alpha=-1), silent=TRUE)
    checkException(stp(Z, des, gamma=c(-1,0.5,2)), silent=TRUE)
    checkException(stp(Z, des, gamma=-1), silent=TRUE)
    checkException(stp(Z, des, th=-1), silent=TRUE)
    checkException(stp(Z, des, th=(-10:10)), silent=TRUE)
}
# Tests for function tc
test_tc <- function ()
{
    # Correct cases
    ## t1, no diff expr.
    Z1 <- matrix(runif(100), 10, 10)
    ## t2, First gene up-reg., second down-reg.
    Z2 <- matrix(runif(100), 10, 10); Z2[1,1:5] <- 10; Z2[2,6:10] <- 10
    des <- list(t1=c(rep(1,5), rep(2,5)), t2=c(rep(1,5), rep(2,5)))
    dat <- list(t1=Z1, t2=Z2)
    res <- tc(dat, des)
    checkTrue(res$iRatios[1] < res$iRatios[2])
    
    # Errors
    ## Missing arguments
    checkException(tc(), silent=TRUE)
    checkException(tc(dat), silent=TRUE)
    checkException(tc(designs=des), silent=TRUE)
    ## Main arguments do not match
    checkException(tc(dat, list(t1=c(rep(1,5), rep(2,5)))), silent=TRUE)
    checkException(tc(list(t1=Z1), des, rep(2,5)), silent=TRUE)
    checkException(tc(dat, list(t1=c(rep(1,5), rep(2,5)), t2=c(rep(1,5), 
        rep(2,5)), t3=c(rep(1,5), rep(2,5)))), silent=TRUE)
    checkException(tc(dat, list(t1=c(rep(1,5), rep(2,5)), t2=c(rep(1,5), 
        rep(2,4)))), silent=TRUE)
    rownames(dat$t1) <- paste("g", 1:10, sep="")
    rownames(dat$t2) <- paste("g", 10:1, sep="")
    checkException(tc(dat, des), silent=TRUE)
    ## Secondary arguments do not match
    rownames(dat$t2) <- paste("g", 1:10, sep="")
    checkException(tc(dat, des, tPoints=1:2), silent=TRUE)
    checkException(tc(dat, des, tPoints=c("t1", "t3", "t3")), silent=TRUE)
    checkException(tc(dat, des, activeTP="t3"), silent=TRUE)
    checkException(tc(dat, des, method=1:2), silent=TRUE)
    checkException(tc(dat, des, method=c("A", "B")), silent=TRUE)
    checkException(tc(dat, des, method=c("active vs complementary",
        "groups conformation", "groups conformation")), silent=TRUE)
}