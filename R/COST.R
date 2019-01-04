#library(mvtnorm)
#library(fields)
#library(copula)
#library(MASS)
#' @importFrom mvtnorm rmvnorm
#' @importFrom mvtnorm rmvt
#' @importFrom copula normalCopula
#' @importFrom copula tCopula
#' @importFrom copula dCopula
#' @importFrom stats quantile
#' @importFrom stats ecdf
#' @importFrom stats mahalanobis
#' @importFrom stats optim
#' @importFrom stats pnorm
#' @importFrom stats pt
#' @importFrom stats qgamma
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @importFrom stats rt
#' @importFrom stats sd


#' @export
Data.COST = function(n,n.total,seed1,coord,par.t)
{
  set.seed(seed1)
  d = nrow(coord)
  V.theta = matrix(c(cos(par.t[1]),sin(par.t[1]),-sin(par.t[1]),cos(par.t[1])),2,2)
  V.lambda = diag(c(1,1/par.t[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }
  R.m = matrix(0,2*d,2*d)
  R.m[1:d,1:d] = R.m[(d+1):(2*d),(d+1):(2*d)] = exp(-par.t[3]*dd^(2*par.t[4]))
  R.m[1:d,(d+1):(2*d)] = exp(-par.t[3]*dd^(2*par.t[4])/(par.t[5]^par.t[4]))/par.t[5]
  R.m[(d+1):(2*d),1:d] = R.m[1:d,(d+1):(2*d)]
  R.11 = R.m[1:d,1:d]

  R.12 = R.m[1:d,(d+1):(2*d)]
  R.21 = R.m[(d+1):(2*d),1:d]
  R.11.inv = solve(R.11)
  B.m = R.21%*%R.11.inv
  Omega.m = R.11-R.21%*%R.11.inv%*%R.12

  dfs = par.t[6]

  if (dfs>=50)
  {
    X.var = matrix(0,d,n.total)
    X.var[,1] = rmvnorm(1,sigma=R.11)

    V.var = rmvnorm(n.total,sigma=Omega.m)
    V.var = t(V.var)

    for (i in 2:n.total)
    {
      X.var[,i] = B.m%*%X.var[,i-1]+V.var[,i]
    }
    alpha=2*coord[,1]+coord[,2]^2
    beta=coord[,1]+coord[,2]
    U.var = pnorm(X.var)
    ytrans.var=qgamma(U.var,shape=alpha,scale=beta)
    ytrans.var = t(ytrans.var)
    Y = ytrans.var[(n.total-n):n.total,]
    X.con = B.m%*%X.var[,n.total-1]
    m1 = 500
    xx = rmvnorm(m1,mean=X.con,sigma=Omega.m)
    yy = qgamma(pnorm(t(xx)),shape=alpha,scale=beta)
    mean.true = apply(yy,1,mean)
  }

  if (dfs<50)
  {
    X.var = matrix(0,d,n.total)
    X.var[,1] = as.vector(rmvt(1,sigma=R.11,df=dfs,delta=rep(0,nrow(R.11)),type = "shifted"))
    for (i in 2:n.total)
    {
      aa.t = t(X.var[,i-1])%*%R.11.inv%*%X.var[,i-1]
      aa.t = as.numeric(aa.t)
      R.11.t = (dfs+aa.t)/(dfs+d)*Omega.m
      V.t = rmvt(1, sigma = R.11.t, df = dfs+d, delta = rep(0, nrow(R.11.t)), type = "shifted")
      X.var[,i] = as.vector(V.t)+as.vector(B.m%*%X.var[,i-1])
    }
    alpha=2*coord[,1]+coord[,2]^2
    beta=coord[,1]+coord[,2]
    U.var = pt(X.var,df=dfs)
    ytrans.var=qgamma(U.var,shape=alpha,scale=beta)
    ytrans.var = t(ytrans.var)
    Y = ytrans.var[(n.total-n):n.total,]
    X.con = B.m%*%X.var[,n.total-1]
    m1 = 500
    aa = t(X.var[,n.total-1])%*%R.11.inv%*%X.var[,n.total-1]
    aa = as.numeric(aa)
    R.11.t = (dfs+aa)/(dfs+d)*Omega.m
    V.t = rmvt(m1, sigma = R.11.t, df = dfs+d, delta = rep(0, nrow(R.11.t)), type = "shifted")
    xx = t(V.t)+as.vector(B.m%*%X.var[,n.total-1])
    yy = qgamma(pt(xx,df=dfs),shape=alpha,scale=beta)
    mean.true = apply(yy,1,mean)
  }
  dat = list(Y.all=Y,mean.true=mean.true)
  return(dat)
}

#' @export
logL.COST.t <- function(par,Y,s.ob){
  d = ncol(Y)
  coord = s.ob
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  R.m = matrix(0,2*d,2*d)
  R.m[1:d,1:d] = R.m[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d,(d+1):(2*d)] = R.m[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  params = R.m[lower.tri(R.m)]
  myCop <- tCopula(param=params, dim = 2*d, dispstr = "un", df=par[6], df.fixed=FALSE)
  R.11 = R.m[1:d,1:d]
  params.marginal = R.11[lower.tri(R.11)]
  myCop.marginal = tCopula(param=params.marginal,dim=d,dispstr="un",df=par[6],df.fixed=FALSE)

  n = nrow(Y)
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)

  L = 0
  aa = log(dCopula(Gn[1,],myCop.marginal))
  L = L+aa
  for(i in 2:n)
  {
    c = log(dCopula(c(Gn[i-1,],Gn[i,]),myCop))
    c.marginal = log(dCopula(Gn[i-1,],myCop.marginal))
    L <- L + c - c.marginal
  }
  return(-L)
}

#' @export
logL.COST.G <- function(par,Y,s.ob){
  d = ncol(Y)
  coord = s.ob
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  R.m = matrix(0,2*d,2*d)
  R.m[1:d,1:d] = R.m[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d,(d+1):(2*d)] = R.m[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  params = R.m[lower.tri(R.m)]
  myCop <- normalCopula(param=params, dim = 2*d, dispstr = "un")
  R.11 = R.m[1:d,1:d]
  params.marginal = R.11[lower.tri(R.11)]
  myCop.marginal = normalCopula(param=params.marginal, dim = d, dispstr = "un")

  n = nrow(Y)
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)
  L = 0
  aa = log(dCopula(Gn[1,],myCop.marginal))
  L = L+aa
  for(i in 2:n)
  {
    c = log(dCopula(c(Gn[i-1,],Gn[i,]),myCop))
    c.marginal = log(dCopula(Gn[i-1,],myCop.marginal))
    L <- L + c - c.marginal
  }
  return(-L)
}

#' @export
logL.CF <- function(par,Yk,dfs){
  deltak = par
  n = length(Yk)
  myCop <- tCopula(deltak, df=dfs, df.fixed=TRUE)
  Gnk = ecdf(Yk)(Yk)*n/(n+1)
  L = 0
  for(i in 2:n)
  {
    L = L+log(dCopula(c(Gnk[i-1],Gnk[i]),myCop))
  }
  return(-L)
}

#' @export
logL.GP <- function(par,Y,s.ob)
  {
  coord = s.ob
  d = nrow(coord)
  n = nrow(Y)
  mu.est = apply(Y,2,mean)
  sigma.est = as.vector(apply(Y,2,sd))
  var.est = sigma.est%*%t(sigma.est)

  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  Sigma = matrix(0,2*d,2*d)
  Sigma[1:d,1:d] = Sigma[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))*var.est
  Sigma[1:d,(d+1):(2*d)] = Sigma[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]*var.est
  Sigma.11.inv = solve(Sigma[1:d,1:d])
  Sigma.12 = Sigma[1:d,(d+1):(2*d)]
  Sigma.21 = Sigma[(d+1):(2*d),1:d]
  B = Sigma.21%*%Sigma.11.inv
  Omega = Sigma[1:d,1:d]-Sigma.21%*%Sigma.11.inv%*%Sigma.12
  Omega.inv = solve(Omega)
  y.cen = t(Y[2:n,])-mu.est-B%*%(t(Y[1:(n-1),])-mu.est)
  L = -(n-1)/2*log(det(Omega))-sum(apply(y.cen,2,function(x){t(x)%*%Omega.inv%*%x/2}))
  return(-L)
}

#' @export
Forecasts.COST.G <- function(par,Y,s.ob,seed1,m,isotropic)
{
  coord = s.ob
  d = nrow(coord)
  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  n = nrow(Y)
  R.m = matrix(0,2*d,2*d)
  R.m[1:d,1:d] = R.m[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d,(d+1):(2*d)] = R.m[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  R.11 = R.m[1:d,1:d]
  R.11.inv = solve(R.11)
  R.12 = R.m[1:d,(d+1):(2*d)]
  R.21 = R.m[(d+1):(2*d),1:d]
  B.m = R.21%*%R.11.inv
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)
  x.n = qnorm(Gn[n,])
  qq = c(0.025,0.975,0.5)
  y.qq = matrix(NA,d,length(qq))
  for (k in 1:d)
  {
    R.m.k = R.m[c(1:d,d+k),c(1:d,d+k)]
    sigma.dk = R.m.k[1:d,d+1]
    mean.k = t(sigma.dk)%*%R.11.inv%*%x.n
    var.k = max(1-t(sigma.dk)%*%R.11.inv%*%sigma.dk,0.001)
    sd.k = var.k^0.5
    q.x.qq = pnorm(qnorm(qq,mean.k,sd.k))
    y.qq[k,] = as.numeric(quantile(Y[,k],q.x.qq,type=6))
  }
  y.draw.random = matrix(0,d,m)
  mean.fore = B.m%*%x.n
  sigma.fore = R.11-B.m%*%R.12
  set.seed(seed1)
  x.draw.random = t(rmvnorm(m,mean=mean.fore,sigma=sigma.fore))

  for (k in 1:d) y.draw.random[k,] = as.numeric(quantile(Y[,k],pmin(pnorm(x.draw.random[k,]),0.999),type=6))
  mean.est = apply(y.draw.random,1,mean)

  result = list(y.qq=y.qq,mean.est=mean.est,y.draw.random=y.draw.random)
  return(result)
}

#' @export
Forecasts.COST.t <- function(par,Y,s.ob,seed1,m,isotropic)
{
  coord = s.ob
  d = nrow(coord)
  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  n = nrow(Y)
  R.m = matrix(0,2*d,2*d)
  R.m[1:d,1:d] = R.m[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d,(d+1):(2*d)] = R.m[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  R.11 = R.m[1:d,1:d]
  R.11.inv = solve(R.11)
  R.12 = R.m[1:d,(d+1):(2*d)]
  R.21 = R.m[(d+1):(2*d),1:d]
  B.m = R.21%*%R.11.inv
  Omega.m = R.11-R.21%*%R.11.inv%*%R.12
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)

  x.n = as.vector(qt(Gn[n,],df=par[6]))
  qq = c(0.025,0.975,0.5)
  y.qq = matrix(NA,d,length(qq))
  for (k in 1:d)
  {
    R.m.k = R.m[c(1:d,d+k),c(1:d,d+k)]
    sigma.dk = R.m.k[1:d,d+1]
    mean.k = as.numeric(t(sigma.dk)%*%R.11.inv%*%x.n)
    Omega.k = max(1-t(sigma.dk)%*%R.11.inv%*%sigma.dk,0.001)
    aa = as.numeric(t(x.n)%*%R.11.inv%*%x.n)
    Omega.k.star = (par[6]+aa)*Omega.k/(par[6]+d)
    aaa = qt(qq,df=par[6]+d)*Omega.k.star^0.5+mean.k
    q.x.qq = pt(aaa,df=par[6])
    y.qq[k,] = as.numeric(quantile(Y[,k],q.x.qq,type=6))
  }
  aa = t(x.n)%*%R.11.inv%*%x.n
  aa = as.numeric(aa)
  R.11.t = (par[6]+aa)/(par[6]+d)*Omega.m
  set.seed(seed1)
  V.t = rmvt(m, sigma = R.11.t, df = par[6]+d, delta = rep(0, nrow(R.11.t)), type = "shifted")
  x.draw.random = t(V.t)+as.vector(B.m%*%x.n)

  y.draw.random = matrix(0,d,m)
  for (k in 1:d) y.draw.random[k,] = as.numeric(quantile(Y[,k],pmin(pt(x.draw.random[k,],df=par[6]),0.999),type=6))
  mean.est = apply(y.draw.random,1,mean)

  result = list(y.qq=y.qq,mean.est=mean.est,y.draw.random=y.draw.random)
  return(result)
}

#' @export
Forecasts.CF <- function(par,Y,seed1,m)
{
  d = ncol(Y)
  n = nrow(Y)
  y.draw.random = matrix(0,d,m)
  dfs = par[1]
  deltas = par[2:(d+1)]

  qq = c(0.025,0.975,0.5)
  y.qq = matrix(NA,d,length(qq))
  for (k in 1:d)
  {
    deltak = deltas[k]
    Yk = Y[,k]
    Gnk = ecdf(Yk)(Yk)*n/(n+1)
    x.kn = qt(Gnk[n],df=dfs)
    mean.k = deltak*x.kn

    var.k = max(1-deltak^2,0.001)
    var.k.star = (dfs+x.kn^2)*var.k/(dfs+1)
    aaa = qt(qq,df=dfs+1)*var.k.star^0.5+mean.k
    q.x.qq = pt(aaa,df=dfs)
    y.qq[k,] = as.numeric(quantile(Yk,q.x.qq,type=6))
    set.seed(seed1)
    x.draw.random = rt(m,df=dfs+1)*var.k.star^0.5+mean.k
    y.draw.random[k,] = as.numeric(quantile(Yk,pmin(pt(x.draw.random,df=dfs),0.999),type=6))
  }
  mean.est = apply(y.draw.random,1,mean)
  result = list(y.qq=y.qq,mean.est=mean.est,y.draw.random=y.draw.random)
  return(result)
}

#' @export
Forecasts.GP <- function(par,Y,s.ob,seed1,m,isotropic)
{
  coord = s.ob
  d = nrow(coord)
  n = nrow(Y)
  mu.est = apply(Y,2,mean)
  sigma.est = as.vector(apply(Y,2,sd))
  var.est = sigma.est%*%t(sigma.est)

  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d,d)
  for (i in 1:d)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  Sigma = matrix(0,2*d,2*d)
  Sigma[1:d,1:d] = Sigma[(d+1):(2*d),(d+1):(2*d)] = exp(-par[3]*dd^(2*par[4]))*var.est
  Sigma[1:d,(d+1):(2*d)] = Sigma[(d+1):(2*d),1:d] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]*var.est
  Sigma.11.inv = solve(Sigma[1:d,1:d])
  Sigma.12 = Sigma[1:d,(d+1):(2*d)]
  Sigma.21 = Sigma[(d+1):(2*d),1:d]
  B = Sigma.21%*%Sigma.11.inv
  Omega = Sigma[1:d,1:d]-Sigma.21%*%Sigma.11.inv%*%Sigma.12
  Omega.inv = solve(Omega)
  Y = t(t(Y)-mu.est)
  y.mean0 = B%*%as.vector(Y[n,])

  qq = c(0.025,0.975,0.5)
  y.qq = matrix(NA,d,length(qq))
  for (k in 1:d)
  {
    y.qq[k,] = qnorm(qq,mean=y.mean0[k],sd=Omega[k,k]^0.5)+mu.est[k]
  }
  set.seed(seed1)
  y.draw.random = t(rmvnorm(m,mean=y.mean0,sigma=Omega))+mu.est #d*m

  mean.est = apply(y.draw.random,1,mean)

  result = list(mean.est=mean.est,y.qq=y.qq,y.draw.random=y.draw.random)
  return(result)
}

#' @export
rank.multivariate = function(y.test,y.random,seed1)
{
  yy = cbind(y.test,y.random)
  mm = ncol(yy)
  pre.rank = rep(0,mm)
  for (j in 1:mm)
  {
    dif.j = (yy-yy[,j]<=0)*1
    pre.rank[j] = sum(apply(dif.j,2,min))
  }
  s.less = sum(pre.rank<pre.rank[1])
  s.eq = sum(pre.rank-pre.rank[1]==0)
  set.seed(seed1)
  rank.multi = s.less+sample(1:s.eq,1)
  return(rank.multi)
}

#' @export
example.forecast <- function(n,n.total,seed1)
{
  m = 500 #number of random draws for forecast
  s.ob = cbind(rep(c(1,3,5)/6,each=3),rep(c(1,3,5)/6,3))
  coord = s.ob
  d = nrow(coord) #number of locations, =9
  par.t = c(0,1,1,0.5,1.5,100) #the true angel, lambda12, c, gamma, eta and degrees of freedom
  #df=100, means data from Gaussian copula

  dat = Data.COST(n,n.total,seed1,coord,par.t)
  Y.all = dat$Y.all
  Y = Y.all[1:n,] #for parameter estimation
  Y.d.test = Y.all[n+1,] #test data for forecast

  #t copula
  #pars.t.COST<-optim(par=c(0,1,2,0.8,2,5),logL.COST.t,Y=Y,s.ob=s.ob,method="L-BFGS-B",
                     #lower=c(0,0.1,0.2,0.05,1.1,1),upper=c(pi/2,50,10,0.95,30,100))$par
  #pars.t.COST <- round(pars.t.COST,4)
  result.COST.t <- Forecasts.COST.t(par.t,Y,s.ob,seed1,m,isotropic=FALSE)
  low.up.COST.t = result.COST.t$y.qq
  COST.t.fore.ECP = (Y.d.test>=low.up.COST.t[,1])*(Y.d.test<=low.up.COST.t[,2])*1
  COST.t.fore.ML = low.up.COST.t[,2]-low.up.COST.t[,1]
  COST.t.fore.ML <- round(COST.t.fore.ML,4)
  COST.t.fore.rank = rank.multivariate(Y.d.test,result.COST.t$y.draw.random,seed1)

  #Gaussian copula
  #pars.G.COST<-optim(par=c(0,1,2,0.8,2),logL.COST.G,Y=Y,s.ob=s.ob,method="L-BFGS-B",
                     #lower=c(0,0.1,0.2,0.05,1.1),upper=c(pi/2,50,10,0.95,30))$par
  #pars.G.COST <- round(pars.G.COST,4)
  result.COST.G <- Forecasts.COST.G(par.t,Y,s.ob,seed1,m,isotropic=FALSE)
  low.up.COST.G = result.COST.G$y.qq
  COST.G.fore.ECP = (Y.d.test>=low.up.COST.G[,1])*(Y.d.test<=low.up.COST.G[,2])*1
  COST.G.fore.ML = low.up.COST.G[,2]-low.up.COST.G[,1]
  COST.G.fore.ML <- round(COST.G.fore.ML,4)
  COST.G.fore.rank = rank.multivariate(Y.d.test,result.COST.G$y.draw.random,seed1)

  #Temporal fit
  #dfs = pars.t.COST[6]
  #pars.CF = c(dfs,rep(0,d))
  #for (k in 1:d)
  #{
    #pars.CF[k+1]<-optim(par=0.5,logL.CF,Yk=Y[,k],dfs=dfs,method="L-BFGS-B",lower=0.05,upper=0.95)$par
  #}
  #pars.CF <- round(pars.CF,4)
  #result.CF <- Forecasts.CF(pars.CF,Y,seed1,m)
  #low.up.CF = result.CF$y.qq
  #CF.fore.ECP = (Y.d.test>=low.up.CF[,1])*(Y.d.test<=low.up.CF[,2])*1
  #CF.fore.ML = low.up.CF[,2]-low.up.CF[,1]
  #CF.fore.ML <- round(CF.fore.ML,4)
  #CF.fore.rank = rank.multivariate(Y.d.test,result.CF$y.draw.random,seed1)

  #Gaussian process
  #pars.GP <- optim(par=c(0,1,2,0.8,2),logL.GP,Y=Y,
                   #s.ob=s.ob,method="L-BFGS-B",lower=c(0,0.1,0.2,0.05,1.1),upper=c(pi/2,50,10,0.95,30))$par
  #pars.GP <- round(pars.GP,4)
  result.GP <- Forecasts.GP(par.t[-6],Y,s.ob,seed1,m,isotropic=FALSE)
  low.up.GP = result.GP$y.qq
  GP.fore.ECP = (Y.d.test>=low.up.GP[,1])*(Y.d.test<=low.up.GP[,2])*1
  GP.fore.ML = low.up.GP[,2]-low.up.GP[,1]
  GP.fore.ML <- round(GP.fore.ML,4)
  GP.fore.rank = rank.multivariate(Y.d.test,result.GP$y.draw.random,seed1)
  result.forecast = list(COST.t.fore.ECP=COST.t.fore.ECP,
                         COST.t.fore.ML=COST.t.fore.ML,COST.t.fore.rank=COST.t.fore.rank,
                         COST.G.fore.ECP=COST.G.fore.ECP,
                         COST.G.fore.ML=COST.G.fore.ML,COST.G.fore.rank=COST.G.fore.rank,
                         GP.fore.ECP=GP.fore.ECP,GP.fore.ML=GP.fore.ML,
                         GP.fore.rank=GP.fore.rank)
  return(result.forecast)
}

##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################

#' @export
Predictions.COST.t <- function(par,Y,s.ob,s.new,isotropic)
{
  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  coord = rbind(s.ob,s.new)
  d.all = nrow(coord)
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d.all,d.all)
  for (i in 1:d.all)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }
  dfs = par[6]
  d = ncol(Y)
  n = nrow(Y)
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)

  R.m = matrix(0,2*d.all,2*d.all)
  R.m[1:d.all,1:d.all] = R.m[(d.all+1):(2*d.all),(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d.all,(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  R.m[(d.all+1):(2*d.all),1:d.all] = R.m[1:d.all,(d.all+1):(2*d.all)]

  h = nrow(s.new)
  pre.CP = pre.ML = rep(0,h)
  qq = c(0.025,0.975,0.5)
  y.qq = x.qq = matrix(0,3,h)

  R.2d = R.m[c(1:d,(d.all+1):(d.all+d)),c(1:d,(d.all+1):(d.all+d))]
  R.2d.inv = solve(R.2d)
  R.2d.new = R.m[c(1:d,(d.all+1):(d.all+d)),c((d.all+d+1):(d.all*2))]
  B.m = t(R.2d.new)%*%R.2d.inv
  Omega.m = R.m[c((d.all+d+1):(d.all*2)),c((d.all+d+1):(d.all*2))]-B.m%*%R.2d.new

  x.c0 = c(qt(Gn[n-1,],df=dfs),qt(Gn[n,],df=dfs))
  x.mean0 = as.vector(B.m%*%x.c0)
  aa0 = t(x.c0)%*%R.2d.inv%*%x.c0
  aa0 = as.numeric(aa0)
  R.11.t0 = (dfs+aa0)/(dfs+2*d)*Omega.m

  for (k in 1:h)
  {
    x.qq[,k] = qt(qq,dfs+2*d)*R.11.t0[k,k]^0.5+x.mean0[k]
  }

  for (k in 1:h)
  {
    dd.2 = dd[c(1:d,d+k),c(1:d,d+k)]
    neighbor = order(dd.2[1:d,d+1])[1:4]
    weight.nei = rep(1/4,4)
    Gn.new.low = min(Y[,neighbor])
    Gn.new.up  = max(Y[,neighbor])
    Gn.new.gri = seq(Gn.new.low,Gn.new.up,length=n)
    Gn.new.all = matrix(0,n,4)
    for (j in 1:4) Gn.new.all[,j] = ecdf(Y[,neighbor[j]])(Gn.new.gri)*n/(n+1)
    Gn.new = Gn.new.all%*%weight.nei

    yy1 = pt(x.qq[,k],df=dfs)
    zz1 = pmin(ecdf(Gn.new)(yy1)*n+1,n)
    y.qq[,k] = Gn.new.gri[zz1]
  }

  return(y.qq)
}

#' @export
Predictions.COST.G <- function(par,Y,s.ob,s.new,isotropic)
{
  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  coord = rbind(s.ob,s.new)
  d.all = nrow(coord)
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d.all,d.all)
  for (i in 1:d.all)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }
  dfs = par[6]
  d = ncol(Y)
  n = nrow(Y)
  Gn = matrix(0,n,d)
  for (k in 1:d) Gn[,k] = ecdf(Y[,k])(Y[,k])*n/(n+1)

  R.m = matrix(0,2*d.all,2*d.all)
  R.m[1:d.all,1:d.all] = R.m[(d.all+1):(2*d.all),(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4]))
  R.m[1:d.all,(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]
  R.m[(d.all+1):(2*d.all),1:d.all] = R.m[1:d.all,(d.all+1):(2*d.all)]

  h = nrow(s.new)
  pre.CP = pre.ML = rep(0,h)
  qq = c(0.025,0.975,0.5)
  y.qq = x.qq = matrix(0,3,h)

  R.2d = R.m[c(1:d,(d.all+1):(d.all+d)),c(1:d,(d.all+1):(d.all+d))]
  R.2d.inv = solve(R.2d)
  R.2d.new = R.m[c(1:d,(d.all+1):(d.all+d)),c((d.all+d+1):(d.all*2))]
  B.m = t(R.2d.new)%*%R.2d.inv
  Omega.m = R.m[c((d.all+d+1):(d.all*2)),c((d.all+d+1):(d.all*2))]-B.m%*%R.2d.new

  x.c0 = c(qnorm(Gn[n-1,]),qnorm(Gn[n,]))
  x.mean0 = as.vector(B.m%*%x.c0)

  for (k in 1:h)
  {
    x.qq[,k] = qnorm(qq)*Omega.m[k,k]^0.5+x.mean0[k]
    dd.2 = dd[c(1:d,d+k),c(1:d,d+k)]
    neighbor = order(dd.2[1:d,d+1])[1:4]
    weight.nei = rep(1/4,4)
    Gn.new.low = min(Y[,neighbor])
    Gn.new.up  = max(Y[,neighbor])
    Gn.new.gri = seq(Gn.new.low,Gn.new.up,length=n)
    Gn.new.all = matrix(0,n,4)
    for (j in 1:4) Gn.new.all[,j] = ecdf(Y[,neighbor[j]])(Gn.new.gri)*n/(n+1)
    Gn.new = Gn.new.all%*%weight.nei

    yy1 = pnorm(x.qq[,k])
    zz1 = pmin(ecdf(Gn.new)(yy1)*n+1,n)
    y.qq[,k] = Gn.new.gri[zz1]
  }

  return(y.qq)
}

#' @export
Predictions.GP <- function(par,Y,s.ob,s.new,isotropic)
{
  h = nrow(s.new)
  if (isotropic==TRUE)
  {
    par = c(0,1,par)
  }
  coord = rbind(s.ob,s.new)
  d.all = nrow(coord)
  n = nrow(Y)
  V.theta = matrix(c(cos(par[1]),sin(par[1]),-sin(par[1]),cos(par[1])),2,2)
  V.lambda = diag(c(1,1/par[2]))
  V.m = V.theta%*%V.lambda%*%t(V.theta)
  dd = matrix(0,d.all,d.all)
  for (i in 1:d.all)
  {
    dd[,i] = sqrt(mahalanobis(x=coord,center=coord[i,],V.m))
  }

  d = ncol(Y)

  sigma.est = as.vector(apply(Y,2,sd))
  sigma.new = rep(0,h)
  for (k in 1:h)
  {
    dd.2 = dd[c(1:d,d+k),c(1:d,d+k)]
    neighbor = order(dd.2[1:d,d+1])[1:4]
    weight.nei = rep(1/4,4)
    sigma.new[k] = (t(weight.nei)%*%(sigma.est[neighbor])^2)^0.5
  }
  sigma.est.all = as.vector(c(sigma.est,sigma.new))

  var.est = sigma.est.all%*%t(sigma.est.all)
  Sigma = matrix(0,2*d.all,2*d.all)
  Sigma[1:d.all,1:d.all] = Sigma[(d.all+1):(2*d.all),(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4]))*var.est
  Sigma[1:d.all,(d.all+1):(2*d.all)] = exp(-par[3]*dd^(2*par[4])/(par[5]^par[4]))/par[5]*var.est
  Sigma[(d.all+1):(2*d.all),1:d.all] = Sigma[1:d.all,(d.all+1):(2*d.all)]

  Sigma.2d = Sigma[c(1:d,(d.all+1):(d.all+d)),c(1:d,(d.all+1):(d.all+d))]
  Sigma.2d.inv = solve(Sigma.2d)
  Sigma.2d.new = Sigma[c(1:d,(d.all+1):(d.all+d)),c((d.all+d+1):(d.all*2))]
  B.m = t(Sigma.2d.new)%*%Sigma.2d.inv
  Omega.m = Sigma[c((d.all+d+1):(d.all*2)),c((d.all+d+1):(d.all*2))]-B.m%*%Sigma.2d.new

  qq = c(0.025,0.975,0.5)

  mean.est = apply(Y,2,mean)
  Y = t(t(Y)-mean.est)

  y.mean0 = B.m%*%c(Y[n-1,],Y[n,])
  mean.new = rep(0,h)
  y.qq = matrix(0,3,h)
  for (k in 1:h)
  {
    dd.2 = dd[c(1:d,d+k),c(1:d,d+k)]
    neighbor = order(dd.2[1:d,d+1])[1:4]
    weight.nei = rep(1/4,4)
    mean.new[k] = weight.nei%*%mean.est[neighbor]
    y.qq[,k] = qnorm(qq,y.mean0[k],Omega.m[k,k]^0.5)+mean.new[k]
  }
  return(y.qq)
}


#########################################
##Example
#########################################

#' @export
example.prediction <- function(n,n.total,seed1)
{
  s.ob = cbind(rep(c(1,3,5)/6,each=3),rep(c(1,3,5)/6,3))
  s.new = cbind(rep(c(1,2)/3,each=2),rep(c(1,2)/3,2))
  coord = rbind(s.ob,s.new)
  d.ob = nrow(s.ob) #number of locations, =9

  par.t = c(0,1,1,0.5,1.5,100) #the true angel, lambda12, c, gamma, eta and degrees of freedom
  #df=100, means data from Gaussian copula

  Y.all = Data.COST(n,n.total,seed1,coord,par.t)$Y.all
  Y = Y.all[1:n,1:d.ob] #for parameter estimation
  Y.newloc.ob = as.vector(Y.all[n,-(1:d.ob)])

  #pars.t.COST<-optim(par=c(0.1,1,2,0.8,2,5),logL.COST.t,Y=Y,s.ob=s.ob,method="L-BFGS-B",
                     #lower=c(0,0.1,0.2,0.05,1.1,1),upper=c(pi/2,20,10,0.95,10,100))$par
  low.up.COST.t <- Predictions.COST.t(par.t,Y,s.ob,s.new,isotropic=FALSE)
  COST.t.pre.ECP = (Y.newloc.ob>=low.up.COST.t[1,])*(Y.newloc.ob<=low.up.COST.t[2,])*1
  COST.t.pre.ML = low.up.COST.t[2,]-low.up.COST.t[1,]
  COST.t.pre.med.error = low.up.COST.t[3,]-Y.newloc.ob

  #pars.G.COST<-optim(par=c(0.1,1,2,0.8,2),logL.COST.G,Y=Y,s.ob=s.ob,method="L-BFGS-B",
                     #lower=c(0,0.1,0.2,0.05,1.1),upper=c(pi/2,20,10,0.95,10))$par
  low.up.COST.G <- Predictions.COST.G(par.t[-6],Y,s.ob,s.new,isotropic=FALSE)
  COST.G.pre.ECP = (Y.newloc.ob>=low.up.COST.G[1,])*(Y.newloc.ob<=low.up.COST.G[2,])*1
  COST.G.pre.ML = low.up.COST.G[2,]-low.up.COST.G[1,]
  COST.G.pre.med.error = low.up.COST.G[3,]-Y.newloc.ob

  #pars.GP <- optim(par=c(0,1,2,0.8,2),logL.GP,Y=Y,
                   #s.ob=s.ob,method="L-BFGS-B",lower=c(0,0.1,0.2,0.05,1.1),upper=c(pi/2,50,10,0.95,30))$par
  low.up.GP <- Predictions.GP(par.t[-6],Y,s.ob,s.new,isotropic=FALSE)
  GP.pre.ECP = (Y.newloc.ob>=low.up.GP[1,])*(Y.newloc.ob<=low.up.GP[2,])*1
  GP.pre.ML = low.up.GP[2,]-low.up.GP[1,]
  GP.pre.med.error = low.up.GP[3,]-Y.newloc.ob

  result.prediction = list(COST.t.pre.ECP=COST.t.pre.ECP,
                           COST.t.pre.ML=COST.t.pre.ML,COST.t.pre.med.error=COST.t.pre.med.error,
                           COST.G.pre.ECP=COST.G.pre.ECP,
                           COST.G.pre.ML=COST.G.pre.ML,COST.G.pre.med.error=COST.G.pre.med.error,
                           GP.pre.ECP=GP.pre.ECP,GP.pre.ML=GP.pre.ML,
                           GP.pre.med.error=GP.pre.med.error)
  return(result.prediction)
}


