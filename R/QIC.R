
#' Quasi Information Criterion
#'
#' Implementation based on MES::QIC.geeglm
#'
#' @param object is a \code{aftgee} fit
#' 
#' @importFrom MASS ginv
#' @export
#' @importFrom stats predict resid
#' 
#' @example inst/examples/ex_aftgee_QIC.R
QIC <- function (object) {
    mu <- log(predict(object)$fit)
    y <- log(object$y)
    quasi <- sum(((y - mu)^2) / -2)
    AIinverse <- ginv(object$gee.vbeta.naiv)
    Vr <- object$gee.vbeta
    trace <- sum(diag(AIinverse %*% Vr))
    params <- length(object$coef.res)
    kpm <- params + length(object$alpha)
    QIC <- -2 * (quasi - trace)
    QICu <- -2 * (quasi - params)
    QICC <- QIC + (2 * kpm * (kpm + 1))/(length(resid(object)) - kpm - 1)
    output <- c(QIC, QICu, quasi, trace, params, QICC)
    names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
    return(output)
}
