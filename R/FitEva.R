#' Evaluate the fitting result.
#' 
#' Evaluate the fitting result by calculating the \strong{\code{TPR}} (ratio of true positives to real positives), \strong{\code{FPR}} (ratio of false positives to real negatives), \strong{\code{precision}} (ratio of true positives to predicted positives) , \strong{\code{accuracy}} (ratio of the sum of true positives and true negatives to total number of predictions) and \strong{\code{F1 score}} (the harmonic average of the precision and recall (TPR)). 
#' 
#' @usage FitEva(Beta.simulate, Beta.fit, Gm, cutoff.sc, cutoff.g)
#' 
#' @param Beta.simulate The simulated or golden standarded \code{B} matrix.
#' @param Beta.fit The fitted \code{B} matrix using \code{\link{PgaMsgl}} or other algorithms.
#' @param Gm Matrix of the group structure of coefficient matrix \code{B}. It is the a matrix of group boundaries, with each row indicating a group, four columns indicate the row-start, row-end, column-start and column-end of the group. The row/column index is 1-based.
#' @param cutoff.sc A number indicates the cutoff of a single coefficient in \code{B} to be true positive.
#' @param cutoff.g A number indicates the cutoff of the norm of a group (a sub-matrix in truth) in \code{B} to be true positive.
#' 
#' @return \item{sc.TPR}{The true positive rate of single coefficients prediction.}
#' @return \item{sc.FPR}{The false positive rate of single coefficients prediction.}
#' @return \item{sc.TTPR}{The "true" true positive rate of single coefficients prediction.}
#' @return \item{sc.precision}{The precision of single coefficients prediction.}
#' @return \item{sc.accuracy}{The accuracy of single coefficients prediction.}
#' @return \item{sc.F1}{The F1 score of single coefficients prediction.}
#' @return \item{g.TPR}{The true positive rate of group prediction.}
#' @return \item{g.FPR}{The false positive rate of group prediction}
#' @return \item{g.precision}{The precision of goup prediction.}
#' @return \item{g.accuracy}{The accuracy of group prediction.}
#' @return \item{g.F1}{The F1 score of group prediction.}
#'      
#' @examples
#' 
#' data(lowD)
#' result <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, lowD$mg, lowD$mc)
#' FitEva_result <- FitEva(lowD$Beta, result$Beta, lowD$Gm, cutoff.sc=0, cutoff.g=0)
#' 
#' @useDynLib PgaMsgl
#' 
#' @export
FitEva <- function(Beta.simulate, Beta.fit, Gm, cutoff.sc, cutoff.g)
{
  if(cutoff.sc < 0) {
    stop("Please indicate a non-negative cutoff of coefficient.")
  }
  TP.sc = sum( which( abs(Beta.fit) > cutoff.sc ) %in% which( abs(Beta.simulate) > cutoff.sc ))
  TTP.sc = sum( Beta.fit*Beta.simulate > cutoff.sc^2 )
  FP.sc = sum(! which( abs(Beta.fit) > cutoff.sc ) %in% which( abs(Beta.simulate) > cutoff.sc ))
  P.sc = sum( abs(Beta.fit) > cutoff.sc )
  
  TN.sc = sum(which( abs(Beta.fit) <= cutoff.sc ) %in% which( abs(Beta.simulate) <= cutoff.sc ))
  
  CP.sc = sum(abs(Beta.simulate) > cutoff.sc)
  CN.sc = sum(abs(Beta.simulate) <= cutoff.sc)
  
  TPR.sc = TP.sc / CP.sc
  TTPR.sc = TTP.sc / CP.sc
  FPR.sc = FP.sc / CN.sc
  precision.sc = TP.sc / P.sc
  accuracy.sc = (TP.sc + TN.sc) / (CP.sc + CN.sc)
  F1.sc = 2 * precision.sc * TPR.sc/ (precision.sc + TPR.sc)
  
  if(cutoff.g < 0) {
    stop("Please indicate a non-negative cutoff of group norm.")
  }
  n.g = dim(Gm)[1]
  simulate.g.norm = rep(0,n.g)
  try.g.norm = rep(0,n.g)
  
  for (i in 1:n.g)
  {
    simulate.g.norm[i]=norm(Beta.simulate[Gm[i,1]:Gm[i,2],Gm[i,3]:Gm[i,4]] , 'F')
    try.g.norm[i]=norm(Beta.fit[Gm[i,1]:Gm[i,2],Gm[i,3]:Gm[i,4]] , 'F')
  }
  
  TP.g = sum(which(try.g.norm > cutoff.g) %in% which(simulate.g.norm > cutoff.g))
  FP.g = sum(! which(try.g.norm > cutoff.g) %in% which(simulate.g.norm > cutoff.g))
  P.g = sum(try.g.norm > cutoff.g)
  TN.g = sum(which(try.g.norm <= cutoff.g) %in% which(simulate.g.norm <= cutoff.g))
  CP.g = sum(simulate.g.norm > cutoff.g)
  CN.g = sum(simulate.g.norm <= cutoff.g)
  
  TPR.g = TP.g / CP.g
  FPR.g = FP.g / CN.g
  precision.g = TP.g / P.g
  accuracy.g = (TP.g + TN.g) / (CP.g + CN.g)
  F1.g = 2 * precision.g * TPR.g/ (precision.g + TPR.g)
  
  return(list(sc.TPR=TPR.sc, 
              sc.FPR=FPR.sc, 
              sc.TTPR=TTPR.sc, 
              sc.pricision=precision.sc, 
              sc.accuracy=accuracy.sc, 
              sc.F1=F1.sc,
              g.TPR=TPR.g, 
              g.FPR=FPR.g, 
              g.pricision=precision.g, 
              g.accuracy=accuracy.g, 
              g.F1=F1.g))
}
