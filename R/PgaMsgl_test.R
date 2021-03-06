#' Proximal Gradient Algorithm for Multi-variate Sparse Group Lasso
#' 
#' Proximal Gradient Algorithm for Multi-variate Sparse Group Lasso.
#' 
#' @usage PgaMsgl_test(XX, YY, B0, model = c("L020v1", "L121"), Gm, mi = 1000, mg, mc, Beta, minlambda = 1e-5, rlambda = 0.98, mintau = 1e-5, rtau = 0.98, cutoff)
#' 
#' @param XX Matrix \code{X} in the model \code{Y=XB}.
#' @param YY Matrix \code{Y} in the model \code{Y=XB}.
#' @param B0 Initial matrix of the coefficient matrix \code{B} in the model \code{Y=XB}.
#' @param model The model for Sparse Group Lasso, \code{L121}, \code{L020v1}, or \code{L020v2}.
#' @param Gm Matrix of the group structure of coefficient matrix \code{B}. It is the a matrix of group boundaries, with each row indicating a group, four columns indicate the row-start, row-end, column-start and column-end of the group. The row/column index is 1-based.
#' @param mi Maximum number of iterations allowed, default value is 1000.
#' @param mg An interger indicates maximum number of groups in matrix \code{B} to be reserved.
#' @param mc An interger indicates maximum number of single coefficients in matrix \code{B} to be reserved.
#' @param Beta The real (simulated) coefficient matrix.
#' @param cutoff The cutoff of convergence.
#' @param minlambda Minimum value of lambda. Only used when \code{model = "L020v2"}, default value is 1e-5.
#' @param rlambda Rate of lambda decrease. Only used when \code{model = "L020v2"}, defalult value is 0.98.
#' @param mintau Minimum value of tau. Only used when \code{model = "L020v2"}, default value is 1e-5.
#' @param rtau Rate of tau decrease. Only used when \code{model = "L020v2"}, default value is 0.98.
#' 
#' @return \item{Beta}{The estimated coefficient matrix \code{B} in the model \code{Y=XB}.}
#' @return \item{Rss}{A vector with length \code{mi} (maximum number of iterations), which is the \code{residual sum of squares (RSS)} at each iteration.}
#' @return \item{Tau}{A vector with length \code{mi} (maximum number of iterations), which is the value of parameter \code{tau} at each iteration step.}
#' @return \item{Lambda}{A vector with length \code{mi} (maximum number of iterations), which is the value of parameter \code{lambda} at each iteration step.}
#' @return \item{iteration.time}{An interger, which is the times of iterations in practice.}
#' @return \item{Rss_relative}{A vector with length \code{mi} (maximum number of iterations), which is \code{Rss} divided by square of \code{Y}'s norm.}
#' @return \item{Error_relative}{A vector with length \code{mi} (maximum number of iterations), which is \code{Error} divided by square of \code{Beta}'s norm.}
#' 
#' @examples
#' data(lowD)
#' result <- PgaMsgl_test(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, lowD$mg, lowD$mc, lowD$Beta, 1e-3)
#' 
#' @useDynLib PgaMsgl
#' 
#' @import Rcpp RcppEigen
#' 
#' @export
PgaMsgl_test <- function(XX, YY, B0, model=c("L020v1", "L121"), Gm, mi=1000, mg, mc, Beta, cutoff, minlambda = 1e-5, rlambda = 0.98, mintau = 1e-5, rtau = 0.98)
{
  if(model=="L020v1"){
    L020v1_test(XX, YY, B0, Gm, mi, mg, mc, Beta, cutoff)
  } else if (model=="L121") {
    L121_test(XX, YY, B0, Gm, mi, mg, mc, Beta, cutoff)
  } else {
    stop("Please indicate the model: L121 or L020v1")
  }
}

