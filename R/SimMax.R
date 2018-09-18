#' Simulate matrices X, B, and Y in the model Y=XB.
#' 
#' Simulate matrices \code{X}, \code{B}, and \code{Y} in the model \code{Y=XB} with a certain group structure and sparsity of \code{B}. With this function, groups can only be simulated with equal sizes. 
#' 
#' @usage SimMax(N, p, q, n.g.row, n.g.col, mg, mc, seed)
#' 
#' @param N Number of samples to be simulated, which is also the number of rows in matrices \code{X} and \code{Y}.
#' @param p Number of covariates, which is also the number of columns of \code{X} or number of rows of \code{B}.
#' @param q Number of response variables, which is also the number of columns of \code{Y} and \code{B}.
#' @param n.g.row Number of groups of rows in \code{B}. It should be divisible by \code{p}.
#' @param n.g.col Number of groups of columns in \code{B}. It should be divisible by \code{q}.
#' @param mg An interger indicates maximum number of non-zero groups in matrix \code{B}.
#' @param mc An interger indicates maximum number of non-zero single coefficients in matrix \code{B}.
#' @param seed Numeric value for \code{set.seed()} for reproduce the random latter.
#' 
#' @return \item{X}{Simulated matrix \code{X} with dimension \code{Nxp} in the model \code{Y=XB}.}
#' @return \item{Beta}{Simulate coefficient matrix \code{B} with dimension \code{pxq} in the model \code{Y=XB}.}
#' @return \item{Y}{Simulated matrix \code{Y} with dimension \code{Nxq} in the model \code{Y=XB}.}
#' @return \item{Gm}{Matrix of the group structure of coefficient matrix \code{B}. It is the a matrix of group boundaries, with each row indicating a group, four columns indicate the row-start, row-end, column-start and column-end of the group. The row/column index is 1-based.}
#' @examples
#' 
#' N <- 25
#' p <- 50
#' q <- 50
#' n.g.row <- 5
#' n.g.col <- 5
#' mg <- 10
#' mc <- 100
#' seed <- 1
#' sim <- SimMax(N, p, q, n.g.row, n.g.col, mg, mc, seed)
#' 
#' @useDynLib PgaMsgl
#' 
#' @export
SimMax <- function (N, p, q, n.g.row, n.g.col, mg, mc, seed)
{
  set.seed(seed) 
  X0 <- matrix(rnorm(N*p, mean=0, sd=1), N, p)
  
  g.row.len = as.integer(p/n.g.row) # number of rows in each group
  g.col.len = as.integer(q/n.g.col) # number of columns in each group
  
  Beta0=matrix(rep(0, p*q), nrow=p, ncol=q, byrow=TRUE)
  n.g = n.g.row * n.g.col # Total number of groups.
  
  # random select rows of groups with non-zero groups 
  nz.g.index <- sort(sample(1:n.g, mg, replace=F), decreasing=FALSE)
  
  col1 <- rep(seq(1, p-g.row.len+1, g.row.len), each=n.g.col)
  col2 <- rep(seq(g.row.len, p, g.row.len),each=n.g.col)
  col3 <- rep(seq(1, q-g.col.len+1, g.col.len), n.g.row)
  col4 <- rep(seq(g.col.len, q, g.col.len), n.g.row)
  Gm <- matrix(c(col1, col2, col3, col4), nrow=n.g, ncol=4, byrow=F)
  
  for (i in nz.g.index)
  {
    grp.row.start=Gm[i,1]
    grp.col.start=Gm[i,3]
    
    # Random select non-zero elements in each group.
    nz.ce=sort(sample(1:(g.row.len*g.col.len), as.integer(mc/mg), replace=F), decreasing=FALSE)
    for (k in nz.ce)
    {
      Beta0[grp.row.start+(k-1)%/%g.col.len, grp.col.start+(k-1)%%g.col.len]=runif(1, -1, 1)
    }
  }
  #  E0=matrix(rnorm(N*q, mean=0, sd=1), N, q)
  
  Y0 <- X0%*%Beta0
  return(list(X=X0, Beta=Beta0, Y=Y0, Gm=Gm))
}