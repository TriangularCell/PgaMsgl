#' Generate an initial Beta matrix for iteration with PgaMsgl.
#' 
#' Generate an initial Beta matrix with simple linear regression of each dependent and independent variable. 
#' 
#' @usage initialBeta(XX, YY, method=p.adjust.methods, cutoff=p.adjust.cutoff)
#' 
#' @param XX Matrix \code{X} in the model \code{Y=XB}.
#' @param YY Matrix \code{Y} in the model \code{Y=XB}.
#' @param method Method for p-value adjustment, acceptable values c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). More details see \code{\link{p.adjust}}. 
#' @param cutoff A number between 0 and 1 indicates the cutoff of adjusted p-value.
#' 
#' @return \item{B0}{The initial Beta matrix generated with simple linear regression with binary elements. The element in it is 1 if the simple linear regression of the corresponding dependent and independent variable passed the \code{p.adjust.cutoff}, and 0 if not.}
#' @return \item{B0_coeff}{The initial Beta matrix generated with simple linear regression. The element in it is the coefficient of the linear regression between variables if it passed the \code{p.adjust.cutoff}, and 0 if not.}
#'      
#' @examples
#' 
#' data(lowD)
#' B0 <- initialBeta(lowD$X, lowD$Y, "BH", 0.01)$B0
#' result <- PgaMsgl(lowD$X, lowD$Y, B0, model="L121", lowD$Gm, lowD$mi, lowD$mg, lowD$mc)
#' FitEva_result <- FitEva(lowD$Beta, result$Beta, lowD$Gm, cutoff.sc=0, cutoff.g=0)
#' 
#' @useDynLib PgaMsgl
#' 
#' @export
initialBeta <- function(XX, YY, method=p.adjust.methods, cutoff=p.adjust.cutoff)
{
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  dim.r <- dim(XX)[2]
  dim.c <- dim(YY)[2]
  
  B0p <- matrix(rep(0,dim.r*dim.c), nrow=dim.r)
  B0b <- matrix(rep(0,dim.r*dim.c), nrow=dim.r)
  
  for(i in 1:dim.r)
  {
    for(j in 1:dim.c)
    {
      B0p[i,j] <- lmp(lm(YY[,j] ~ XX[,i]))
      B0b[i,j] <- summary(lm(YY[,j] ~ XX[,i]))$coefficients[2]
    }
  }
  
  B0 <- matrix(rep(0,dim.r*dim.c), nrow=dim.r)
  B0[which(p.adjust(B0p, method) < cutoff)] <- 1
  
  B0_coeff <- matrix(rep(0,dim.r*dim.c), nrow=dim.r)
  B0_coeff[which(p.adjust(B0p, method) < cutoff)] <- B0b[which(p.adjust(B0p, method) < cutoff)]
  
  return(list(B0=B0, B0_coeff=B0_coeff))
}
