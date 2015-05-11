tc <-
function (data, designs, tPoints = NULL, 
    method = c("active vs complementary", "groups conformation"), 
    activeTP = NULL, alpha = 0.05, B = 100, lambda = 0.5, 
    PER = FALSE, BCa = FALSE, gamma = 0.95, R = 1000, ...) 
{
    if (length(data) == 1) 
        stop("For a single time point, use stp function.")
    if (length(data) != length(designs)) 
        stop("The number of time points must be equal in data, designs.")
    if (is.null(tPoints)) 
        tPoints <- if (is.null(names(data))) 
            paste("t", 1:length(data), sep = "")
        else names(data)
    else {
        if (length(data) != length(tPoints)) 
            stop(paste("Time points in tPoints do not",
                "correspond to data sets in data."))
        if (!is.character(tPoints)) 
            stop("Time points must be character names.")
    }
    tps <- length(tPoints)
    n <- nrow(data[[1]])
    gNames <- rownames(data[[1]])
    for (tp in 1:tps) {
        data[[tp]] <- as.matrix(data[[tp]])
        if (ncol(data[[tp]]) != length(designs[[tp]])) 
            stop("Replicates in Z do not correspond to classes in design.")
        if (length(table(designs[[tp]])) != 2) 
            stop(paste("There are more than 2 classes in design.",
                "Only two class comparison is available."))
        if (nrow(data[[tp]]) != n) 
            stop("Number of genes are not equal in all time points.")
        if (sum(rownames(data[[tp]]) != gNames) > 0) 
            stop("Different order or naming of genes between timepoints.")
    }
    if (!is.character(method) || length(method) > 2 || sum(method == 
        c("active vs complementary", "groups conformation")) == 
        0) 
        stop("Argument method not valid.")
    if (!is.null(activeTP)) {
        if (length(activeTP) != 1) 
            stop("There must be only one active time point.")
        if (sum(activeTP == tPoints) == 0) 
            stop(paste("The active time point does not correspond",
                "to any of the time points."))
    }
    ACT <- sum(method == "active vs complementary") > 0
    GCT <- sum(method == "groups conformation") > 0
    ac <- vector("list", tps)
    for (tp in 1:tps) {
        ac[[tp]] <- ac(data[[tp]], designs[[tp]])
    }
    iratios <- rep(0, tps)
    for (tp in 1:tps) {
        iratios[tp] <- var(ac[[tp]][, 2])/eigen(cor(data[[tp]]))$values[1]
    }
    iratios <- t(iratios)
    colnames(iratios) <- tPoints
    rownames(iratios) <- ""
    gct <- as.list(rep("Not computed", tps))
    names(gct) <- tPoints
    act <- as.list(rep("Not computed", tps))
    names(act) <- tPoints
    activeTP <- if (is.null(activeTP) && ACT) 
        order(iratios, decreasing = TRUE)[1]
    else if (ACT) 
        as.integer(activeTP)
    else NA
    if (ACT && GCT) {
        for (tp in 1:tps) {
            gct[[tp]] <- stp(data[[tp]], designs[[tp]], alpha, 
                B, lambda, th = NULL, PER, BCa, gamma, R, ...)
        }
        for (tp in 1:tps) {
            act[[tp]] <- gct[[activeTP]]
            act[[tp]]$ac <- gct[[tp]]$ac
        }
    }
    else if (ACT) {
        act[[activeTP]] <- stp(data[[activeTP]], designs[[activeTP]], 
            alpha, B, lambda, th = NULL, PER, BCa, gamma, R, 
            ...)
        for (tp in (1:tps)[-activeTP]) {
            act[[tp]] <- act[[activeTP]]
            act[[tp]]$ac <- ac(data[[tp]], designs[[tp]])
        }
    }
    else {
        for (tp in 1:tps) {
            gct[[tp]] <- stp(data[[tp]], designs[[tp]], alpha, 
                B, lambda, th = NULL, PER, BCa, gamma, R, ...)
        }
    }
    res <- list(iRatios = iratios, gct = gct, act = act, activeTP = activeTP, 
        tPoints = tPoints, call = match.call())
    class(res) <- "TC"
    res
}
