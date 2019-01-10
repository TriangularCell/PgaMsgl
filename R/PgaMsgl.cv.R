#' Cross validation of hyperparameters of PgaMsgl
#'
#' Cross validation of hyperparameters \code{mg} and \code{mc} of \code{\link{PgaMsgl}}.
#' 
#' @usage PgaMsgl.cv(XX, YY, B0, model = c("L020v1", "L020v2", "L121"), Gm, mi = 1000, mg.v, mc.v, minlambda = 1e-5, rlambda = 0.98, mintau = 1e-5, rtau = 0.98, fold = 5, seed = 1, ncores)
#' 
#' @param XX Matrix \code{X} in the model \code{Y=XB}. Remember to use \code{as.matrix()} if your input is a data frame.
#' @param YY Matrix \code{Y} in the model \code{Y=XB}. Remember to use \code{as.matrix()} if your input is a data frame.
#' @param B0 Initial matrix of the coefficient matrix \code{B} in the model \code{Y=XB}. Remember to use \code{as.matrix()} if your input is a data frame.
#' @param model The model for Sparse Group Lasso, \code{L121}, \code{L020v1}, or \code{L020v2}.
#' @param Gm Matrix of the group structure of coefficient matrix \code{B}. It is the a matrix of group boundaries, with each row indicating a group, four columns indicate the row-start, row-end, column-start and column-end of the group. The row/column index is 1-based.
#' @param mi Maximum number of iterations allowed, default value is 1000.
#' @param mg.v A vector indicates maximum number of groups in matrix \code{B} to be reserved to be cross validated.
#' @param mc.v A vector indicates maximum number of single coefficients in matrix \code{B} to be reserved to be cross validated.
#' @param minlambda Minimum value of lambda. Only used when \code{model = "L020v2"}, default value is 1e-5.
#' @param rlambda Rate of lambda decrease. Only used when \code{model = "L020v2"}, defalult value is 0.98.
#' @param mintau Minimum value of tau. Only used when \code{model = "L020v2"}, default value is 1e-5.
#' @param rtau Rate of tau decrease. Only used when \code{model = "L020v2"}, default value is 0.98.
#' @param fold Number of fold for k-fold cross validation.
#' @param seed Numeric value for \code{set.seed()} when generate test or evaluate samples for k-fold cross validation. 
#' @param ncores The number of cores to use for parallel execution. A parameter of \code{registerDoMC()} in \code{doMC} package.
#' 
#' @return \item{rss}{A vector of the \code{residual sum of squares (RSS)} of model fitting with each combination of parameters.}
#' @return \item{RMSE}{A vector of the \code{root-mean-square error (RMSE)} of model fitting with each combination of parameters.}
#' @return \item{Rsquare}{A vector of the \eqn{R^{2}}{R^2} of model fitting with each combination of parameters.}
#' @return \item{rss.matr}{A matrix of the \code{residual sum of squares (RSS)} of model fitting with each combination of parameters.}
#' @return \item{RMSE.matr}{A matrix of the \code{root-mean-square error (RMSE)}E of model fitting with each combination of parameters.}
#' @return \item{Rsquare.matr}{A matrix of the \eqn{R^{2}}{R^2} of model fitting with each combination of parameters.}
#' @return \item{mg.v}{A vector of values of the parameter \code{mg} tested.}
#' @return \item{mc.v}{A vector of values of the parameter \code{mc} tested.}
#' 
#' @author Yiming Qin
#' 
#' @examples
#' data(lowD)
#' 
#' mg.v <- seq(from=0.01*100, to=0.8*100, by=0.02*100)
#' mc.v <- seq(from = 0.01*2500, to = 0.5*2500, by = 0.01*2500)
#' 
#' result.cv <- PgaMsgl.cv(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, mg.v, mc.v, fold=5, seed=1, ncores=4)
#' 
#' grp.max.cv.rss <- result.cv$mg.v[which.min(result.cv$rss)]
#' coe.max.cv.rss <- result.cv$mc.v[which.min(result.cv$rss)]
#' 
#' grp.max.cv.rmse <- result.cv$mg.v[which.min(result.cv$RMSE)]
#' coe.max.cv.rmse <- result.cv$mc.v[which.min(result.cv$RMSE)]
#' 
#' grp.max.cv.r2 <- result.cv$mg.v[which.max(result.cv$Rsquare)]
#' coe.max.cv.r2 <- result.cv$mc.v[which.max(result.cv$Rsquare)]
#' 
#' OR: One can select hyperparameters following the one standard error rule.
#' 
#' index.minrss <- which.min(result.cv$rss)
#' grp.max.cv.rss <- result.cv$mg.v[ min( which( result.cv$rss < result.cv$rss[index.minrss] + result.cv$rss.se[index.minrss] ) ) ]
#' coe.max.cv.rss <- result.cv$mc.v[ min( which( result.cv$rss < result.cv$rss[index.minrss] + result.cv$rss.se[index.minrss] ) ) ]
#' 
#' index.minrmse <- which.min(result.cv$RMSE)
#' grp.max.cv.rmse <- result.cv$mg.v[ min( which( result.cv$RMSE <= result.cv$RMSE[index.minrmse] + result.cv$RMSE.se[index.minrmse] ) ) ]
#' coe.max.cv.rmse <- result.cv$mc.v[ min( which( result.cv$RMSE <= result.cv$RMSE[index.minrmse] + result.cv$RMSE.se[index.minrmse] ) ) ]
#' 
#' index.minr2 <- which.max(result.cv$Rsquare)
#' grp.max.cv.r2 <- result.cv$mg.v[ min( which( result.cv$Rsquare >= result.cv$Rsquare[index.minr2] - result.cv$Rsquare.se[index.minr2] ) ) ]
#' coe.max.cv.r2 <- result.cv$mc.v[ min( which( result.cv$Rsquare >= result.cv$Rsquare[index.minr2] - result.cv$Rsquare.se[index.minr2] ) ) ]
#' 
#' system.time(try1 <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, 10, 120))
#' system.time(try2 <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, grp.max.cv.rss, coe.max.cv.rss))
#' system.time(try3 <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, grp.max.cv.rmse, coe.max.cv.rmse))
#' system.time(try4 <- PgaMsgl(lowD$X, lowD$Y, lowD$B0, model="L121", lowD$Gm, lowD$mi, grp.max.cv.r2, coe.max.cv.r2))
#' 
#' @useDynLib PgaMsgl
#' 
#' @import Rcpp RcppEigen doParallel foreach parallel
#' 
#' @export
PgaMsgl.cv <- function(XX, YY, B0, model=c("L020v1", "L020v2", "L121"), Gm, mi=1000, mg.v, mc.v, minlambda = 1e-5, rlambda = 0.98, mintau = 1e-5, rtau = 0.98, fold=5, seed=1, ncores=NULL)
{
  requireNamespace("doParallel", quietly = TRUE)
  if(is.null(ncores)) {
    registerDoParallel()
  } else if (ncores > detectCores()) {
    registerDoParallel(cores=detectCores())
  } else {
    registerDoParallel(cores=ncores)
  }
  
  N=NROW(XX) # number of samples
  mg.v <- sort(mg.v)
  mc.v <- sort(mc.v)
  
  index.cv <- NULL
  set.seed(seed)
  sample.cv = sample(1:N, N, replace = F)
  N.cv=floor(N/fold)
  for (i in 1:fold)
  {
    index.cv[[i]] <- sample.cv[((i - 1) * N.cv + 1) : (i * N.cv)]
  }
  
  n.mg.v <- length(mg.v)
  n.mc.v <- length(mc.v)
  grp.max.all <- rep(mg.v, each=n.mc.v)
  coe.max.all <- rep(mc.v, n.mg.v)
  
  #	rss.test <- matrix(rep(0, fold * n.mg.v * n.mc.v), nrow=fold)
  cs <- ceiling((fold*n.mg.v*n.mc.v) / getDoParWorkers())
  opt <- list(chunkSize=cs)
  
  rcom <- function(x, ...)
  {
    mapply(rbind,x,...,SIMPLIFY=FALSE)
  }
  ccom <- function(x,...)
  {
    mapply(c,x,...,SIMPLIFY=FALSE)
  }
  cv <-
    foreach(f=1:fold, .combine='rcom') %:%
    foreach(j=1:(n.mg.v*n.mc.v), .combine='ccom', .options.mpi=opt) %dopar% 
    {
      X.train <- XX[-index.cv[[f]],]
      Y.train <- YY[-index.cv[[f]],]
      X.test <- XX[index.cv[[f]],]
      Y.test <- YY[index.cv[[f]],]
      
      grp.max = grp.max.all[j]
      coe.max = coe.max.all[j]
      if (model == "L020v1"){
        result <- L020v1(X.train, Y.train, B0, Gm, mi, grp.max, coe.max)
      } else if (model == "L020v2") {
        result <- L020v2(X.train, Y.train, B0, Gm, mi, grp.max, coe.max, minlambda, rlambda, mintau, rtau)
      } else if (model == "L121") {
        result <- L121(X.train, Y.train, B0, Gm, mi, grp.max, coe.max)
      } else {
        stop("Please indicate the model: L121, L020v1, or L020v2")
      }
      beta <- result$Beta
      rss.temp <- norm(as.matrix(X.test%*%beta-Y.test),'F')^2
      RMSE.temp <- sqrt(rss.temp/(NROW(Y.test)*NCOL(Y.test)))
      Rsquare.temp <- 1-rss.temp/(norm(as.matrix(Y.test-mean(Y.test)),'F'))^2
      list(rss.temp, RMSE.temp, Rsquare.temp)
    }
  rss.test <- cv[[1]]
  RMSE.test <- cv[[2]]
  Rsquare.test <- cv[[3]]
  
  rss <- colSums(rss.test)/fold
  RMSE <- colSums(RMSE.test)/fold
  Rsquare <- colSums(Rsquare.test)/fold
  
  fun.se <- function(x) sd(x)/sqrt(length(x))
  rss.se <- apply(rss.test, 2, fun.se)
  RMSE.se <- apply(RMSE.test, 2, fun.se)
  Rsquare.se <- apply(Rsquare.test, 2, fun.se)
  
  return(list(rss=rss, RMSE=RMSE, Rsquare=Rsquare, rss.se=rss.se, RMSE.se=RMSE.se, Rsquare.se=Rsquare.se, rss.matr = rss.test, RMSE.matr = RMSE.test, Rsquare.matr = Rsquare.test, mg.v=grp.max.all, mc.v=coe.max.all))
}