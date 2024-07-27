##' \eqn{\ell_2} Norm and Sup Norm Based Tests
##'
##' This function use both \eqn{\ell_2} and Sup norm based statistics
##' to infer if covariance between two random vectors are zero or not.
##'
##' @param res1,res2  Two data matrice with the same sample size. The goal is to test
##' if the covariance between them are 0 or not.
##' @return Returns an list with \item{tstats_sup}{The Sup norm based test statistic.}
##' \item{pval_sup}{The corresponding p-value for the Sup norm based statistic. The
##' reference distribution is Gumbel.}
##' \item{tstats_chisq}{The \eqn{\ell_2} norm based test statistic.}
##' \item{pval_chisq}{The corresponding p-value for the \eqn{\ell_2} norm based statistic. The
##' reference distribution is weighted Chi-square.}
##' \item{cov_selfnorm}{The entry-wise self-normalized partial covariance matrix, based
##' on which, the \eqn{\ell_2} and Sup norm based statistics are constructed.}
##' \item{weights}{The estimated weights for the weighted Chi-square of reference.}
##' @export

TestLS=function(res1,res2){
##function used to calculate p-value for weighted chisq
##(see Liu et al. "A New Chi-Square Approximation to the Distribution of
##Non-Negative Definite Quadratic Forms in Non-Central Normal Variables" )
  chisqtail=function(lam,t) {
    c1=sum(lam);c2=sum(lam^2);c3=sum(lam^3);c4=sum(lam^4)
    s1=c3/c2^{1.5};s2=c4/c2^2
    a=1/(s1-sqrt(max(s1^2-s2,0)))
    del=round(s1*a^3-a^2,6)
    l=a^2-2*del

    ts=(t-c1)/sqrt(2*c2)
    ux=l+del
    sigx=sqrt(2)*a
    tail=pchisq(q=ts*sigx+ux,df=l,ncp=del,lower.tail=FALSE)

    return(tail)
  }
############################################
  p1=dim(res1)[2];p2=dim(res2)[2];N=dim(res1)[1]

  dd=NULL
  for (i in 1:dim(res1)[2]) {

    for (j in 1:dim(res2)[2]){
      v=res1[,i]*res2[,j]
      dd=cbind(dd,v)
    }

  }
  w=apply(dd,2,mean)
  sig=diag(coop::covar(dd))

  ## sup based stats
  wsq=w^2/(sig/N)
  tstats_sup=max(wsq)-4*log(max(p1,p2))+log(log(max(p1,p2)))
  pval_sup=1-exp(-1/sqrt(8*pi)*exp(-tstats_sup/2))

  ## L2 based stats
  ws=w/sqrt(sig/N)
  sqsum=sum(ws^2)
  lam=eigen(coop::pcor(dd))$values

  tstats_chisq=sqsum
  pval_chisq=chisqtail(lam,sqsum)

  ## Save the entry-wise self-normalized partial covariance matrix
  mtn_selfnorm=matrix(ws,nrow=p1,ncol=p2,byrow=T)

  return(list(tstats_sup=tstats_sup,pval_sup=pval_sup,tstats_chisq=tstats_chisq,
              pval_chisq=pval_chisq,cov_selfnorm=mtn_selfnorm,weights=lam))
}

