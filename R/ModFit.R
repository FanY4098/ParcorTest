##' Model Estimation with either Group or Individual Penalties
##'
##' This function fits the regression model by regularized method to get the fitted
##' residuals and coefficients.
##'
##' @param x Predictors, a matrix of dimensions N * p. There is a group structure
##' within x, and columns that belong to the same group should be adjacent to each other.
##' @param y Response variables, it can be either a vector (univariate) or a
##' matrix (multivariate).
##' @param xGsz Size of each sub-group in x, it is a vector, the length of which is the
##' number of groups. To be more specific, if \eqn{xGsz=(g_1,..,g_m)}, then it means
##' the first \eqn{g_1} columns forms a group, the next \eqn{g_2} colums form another group, and so on.
##' @param lam1 Individual penalty (a scalar), imposed to all coefficients.
##' @param lamG Group penalites (a vector), each entry is the group penalty imposed
##' on a sub-group of coefficients, corresponding to the sub-group in x.
##'
##' @return Returns an list with \item{res}{The fitted residuals. It is a matrix in the
##' case of multivariate regression. Otherwise it is just a vector.}
##' \item{Beta}{The fitted coefficients. It can be a matrix (multivariate case)
##' or a vector (univariate case)}
##' @export




ModFit=function(y,x,xGsz,lam1,lamG) {

  N=dim(x)[1] #sample size
  y=as.matrix(y) # in case y is N by 1

  GarrEnds <-cumsum(xGsz)-1
  GarrStarts <-cumsum(xGsz)+1
  GarrStarts=c(1,GarrStarts[-length(xGsz)])-1

  RarrStarts <-0
  RarrEnds <-dim(y)[2]-1

  P <- sum(xGsz)
  Q <- dim(y)[2]
  G<-length(xGsz) # number of groups in x
  R <- 1          # number of groups in y
  gmax <- 1
  cmax <- 400


  tmp <- MSGLasso::FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQgrps <- tmp$PQgrps
  tmp1 <- MSGLasso::Cal_grpWTs(P, Q, G, R, gmax, PQgrps)
  grpWTs <- tmp1$grpWTs
  tmp2 <- MSGLasso::FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps <- tmp2$GRgrps
  grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)

  P.x <- P
  Q.y <- Q
  Pen_L <- matrix(rep(1,P.x*Q.y),P.x,Q.y, byrow=TRUE)
  Pen_G<-matrix(rep(1,G),G,1, byrow=TRUE)



  ############################################
  MSGLassolam1=lam1
  MSGLassolamG.m <- matrix(lamG,G,R,byrow=TRUE)


  model<- MSGLasso::MSGLasso(x, y, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps, grp_Norm0, MSGLassolam1, MSGLassolamG.m, Beta0=NULL)

  res=y-x%*%model$Beta
  Beta=model$Beta

  return(list(res=res,Beta=Beta))

}
