rbern = function(n,prob)
{
  x = runif(n,min=0,max=1)
  x.bern = ifelse(x<=prob,1,0)
  return(x.bern)
}
get.y = function(betas, design.x, ICCy, M.psu, psu.size){
  dim(betas) = c(length(betas), 1)   #make it matrix of ncol=1
  odds = exp( design.x %*% betas)
  N = nrow(design.x)
  p = odds/(1+odds) #Pr(y=1|x)
  if (ICCy==0) {
    y = ifelse(runif(N)<=p,1,0)
  }
  if (ICCy!=0){
    ei0 = rep(rnorm(M.psu),psu.size)
    eij = rnorm(N)
    Uij = rbern(N,sqrt(ICCy))
    thetaij = qnorm(p)
    threshold = Uij*ei0 + (1-Uij)*eij
    y = ifelse(threshold <= thetaij,1,0)
  }
  return(y)
}
#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps=function(popul,Msize, n){
  prob.s = Msize/sum(Msize)
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
    sam.pps.data = as.data.frame(popul[pps.samID,])
    names(sam.pps.data) = names(popul)
  }else{sam.pps.data = popul[pps.samID,]}
  sam.pps.data$wt = sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}
#################################################################################
samp.slct = function(seed, fnt.pop, n, psu.name=NULL, m.psu=NULL, dsgn, size = NULL, size.I = NULL){
  #fnt.pop = pop; seed = seed.sim[1]; psu.name="psu"; 
  #dsgn = "pps"; 
  #size = "odds.s"; 
  #dsgn = "srs-pps"; 
  #dsgn = "pps-pps"; 
  #size.I = aggregate(fnt.pop[,size], list(fnt.pop[,psu.name]), sum)[,2]
  set.seed(seed)
  N = nrow(fnt.pop)
  # one-ste sample design
  if(dsgn=="pps"){
    samp = sam.pps(fnt.pop,fnt.pop[,size], n)
  }
  # two-stage cluster sampling design (informative design at the second stage)
  if(dsgn == "srs-pps"){
    #-- first stage: select clusters by srs
    size.psu = as.numeric(table(fnt.pop[,psu.name]))
    M.psu = length(size.psu)
    # m.psu = 25
    # n = 1000
    index.psuI = sample(1:M.psu, m.psu, replace = F)
    sample.I = fnt.pop[fnt.pop[,psu.name] %in% index.psuI,]
    sample.I$wt.I = M.psu/m.psu
    #-- second stage: select subj.samp within selected psus by pps 
    # calcualte the size for the second stage
    samp=NULL
    for (i in 1: m.psu){
      popn.psu.i= sample.I[sample.I[,psu.name]==index.psuI[i],]
      samp.i = sam.pps(popn.psu.i,popn.psu.i[,size], n/m.psu)
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  if(dsgn == "pps-pps"){
    #-- first stage: select clusters by pps
    size.psu = as.numeric(table(fnt.pop[,psu.name]))
    M.psu = length(size.psu)
    # Clt.samp = 25
    # n = 1000
    index.psuI = sam.pps(matrix(1:M.psu,,1),size.I, m.psu)
    index.psuI = index.psuI[order(index.psuI[,1]),]  #sort selected psus
    sample.I = fnt.pop[fnt.pop[,psu.name]%in% index.psuI[,1],]
    sample.I = sample.I[order(sample.I[,psu.name]),]
    sample.I$wt.I = rep(index.psuI[,'wt'], size.psu[index.psuI[,1]])
    #-- second stage: select subj.samp within selected psus by pps 
    # calcualte the size for the second stage
    samp=NULL
    for (i in 1: m.psu){
      popn.psu.i = sample.I[sample.I$psu==index.psuI[i,1],    ]
      samp.i = sam.pps(popn.psu.i,popn.psu.i[,size], n/m.psu)
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  rownames(samp) = as.character(1:dim(samp)[1])
  return(samp)      
}



greg.f = function(samp, wt0, N.hat, aux.mtx=NULL, aux.tot, f_w=T){
  n=length(samp[,wt0]);n
  if(!is.null(aux.mtx)) samp = cbind(samp, aux.mtx)
  ds.wt = svydesign(ids=~1, data=samp, weights=as.formula(paste0("~",wt0)))
  calib.fm = as.formula(paste0("~", paste0(names(aux.tot)[-1],collapse="+")))
  V.hat = c(sum(samp[,wt0])/N.hat, svytotal(calib.fm,ds.wt))
  V.hat[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]=V.hat[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  V = aux.tot
  V[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]=V[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  v.mtx = model.matrix(calib.fm,samp)
  v.mtx[,1] = aux.mtx[,1]
  v.mtx[,!grepl( "Delta_beta", names(V.hat), fixed = TRUE)] = v.mtx[,!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  
  vWv_inv = solve(t(samp[,wt0]*v.mtx)%*%v.mtx)
  f = c(t(1+(V-V.hat)%*%vWv_inv%*%t(v.mtx)))
  #f[f<0]=1e-10
  if(f_w){
    f_w1  = (V-V.hat)%*%vWv_inv
    vWv_w = lapply(1:n, function(i) outer(v.mtx[i,], v.mtx[i,]))
    f_w2  = vWv_inv%*%t(v.mtx)
    f_w = -sapply(1:n, function(i) f_w1%*%vWv_w[[i]]%*%f_w2)-
      (v.mtx)%*%vWv_inv%*%t(v.mtx)
    return(list(f = c(f), f_w = f_w))
  }else{return(list(f=c(f)))}
  
}
################################################################################################
beta.est.var = function(ds, fit, fit.i, sub, sub.wt, structure="combine"){
  #ds = ds.all; fit = as.formula(fit.y); fit.i = as.formula(gsub("x4", "x4.i", fit.y))
  #cyc = "cycle"; sub.cyl=1; method="linear"
  # design variables
  samp = ds$variables
  samp$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt = as.numeric(unlist(1/ds$allprob))
  samp0 = samp[samp[,sub]==1,]
  samp0$wt = samp0[,sub.wt]
  rm(ds)
  # outcome variable
  fit = as.formula(fit); fit.i = as.formula(fit.i)
  y = all.vars(fit)[1]
  # sampling design for the combined sample
  ds.all = svydesign(ids=~psu, strata = ~strata, weights=samp$wt, data=samp)
  # outcome model using surrogate variables
  imp.mdl = svyglm(fit.i, ds.all, family="binomial")
  p.tilde = imp.mdl$fitted.values
  x.tilde = model.matrix(imp.mdl)
  beta.tilde = imp.mdl$coefficients
  # influence functions
  U.tilde_beta.tilde=-t(c(samp$wt*p.tilde*(1-p.tilde))*x.tilde)%*%x.tilde
  v.mtx = t(solve(U.tilde_beta.tilde)%*%
          t(c(samp[,y]-p.tilde)*x.tilde))
  colnames(v.mtx)[1] = "x0" 
  colnames(v.mtx) = paste0("delta.", gsub(":", ".", colnames(v.mtx))) 

  v.mtx0 = v.mtx[samp[,sub]==1,]
  # subsample estimate of influence total 
  #Vs0.hat = c(samp0$wt)%*%v.mtx0
  # combined sample estimate of influence total
  VS.hat  = c(samp$wt)%*%v.mtx
  # GREG adjustment factor
    greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=VS.hat#, Vs0.hat=Vs0.hat
                      )
    f = greg.out$f
    f[f<0]=0
    eta = greg.out$eta
  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  if(structure=="combine"){
    #partial derivative of the estimating equation system, w.r.t. the sample weights
    U_w = matrix(0,nrow(samp), ncol(X))
    U_w[samp[,sub]==1,] = samp0$wt.c*c(samp0[,y]-p)*X
    S1_w = matrix(0,nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[,sub]==1,] = samp0$wt.c*v.mtx0
    S_w  = S1_w-samp$wt*v.mtx
    U.tilde_w = samp$wt*c(samp[,y]-p.tilde)*x.tilde
    Phi_w = cbind(U_w, S_w, U.tilde_w)
    #partial derivative of the estimating equation system, w.r.t. the parameters
    #U_beta.fun = function(x){
    #  odds.p = exp(X%*%x)
    #  p1 = odds.p/(1+odds.p)
    #  c(samp0$wt.c*c(samp0[,y]-p1))%*%X
    #}
    #U_beta.fun(x=beta.est)
    #U_beta = jacobian(func=U_beta.fun, x = beta.est)
    #U_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    #}
    #U_eta.fun(x=eta)
    #U_eta = jacobian(func=U_eta.fun, x = eta)
    U_eta = t(t(c(samp0$wt*c(samp0[,y]-p))*v.mtx0)%*%X)
    
    U_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde%*%x)
      p.tilde1 = odds.p.tilde/(1+odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt*p.tilde1*(1-p.tilde1))*x.tilde)%*%x.tilde
      v.mtx1 = t(solve(U.tilde_beta)%*%
                   t(c(samp[,y]-p.tilde1)*x.tilde))
      v.mtx10 = v.mtx1[samp[,sub]==1,]
      f1 = c(t(1+eta%*%t(v.mtx10)))
      t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    }
    U_beta.tilde.fun(x=beta.tilde)
    U_beta.tilde = jacobian(func=U_beta.tilde.fun, x = beta.tilde)
    
    S_beta = matrix(0,length(beta.est),length(beta.est))
    #S_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  c(f1*samp0$wt)%*%v.mtx0-VS.hat
    #}
    #S_eta.fun(x=eta)
    #S_eta = jacobian(func=S_eta.fun, x = eta)
    S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
    S_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde%*%x)
      p.tilde1 = odds.p.tilde/(1+odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt*p.tilde1*(1-p.tilde1))*x.tilde)%*%x.tilde
      v.mtx1 = t(solve(U.tilde_beta)%*%
                   t(c(samp[,y]-p.tilde1)*x.tilde))
      v.mtx10 = v.mtx1[samp[,sub]==1,]
      f1 = c(t(1+eta%*%t(v.mtx10)))
      c(f1*samp0$wt)%*%v.mtx10-c(samp$wt)%*%v.mtx1
    }
    #S_beta.tilde.fun(x=beta.tilde)
    S_beta.tilde = jacobian(func=S_beta.tilde.fun, x = beta.tilde)
    
    U.tilde_beta = matrix(0,length(beta.tilde),length(beta.est))
    U.tilde_eta  = matrix(0,length(beta.tilde),length(beta.tilde))
    #U.tilde_beta.tilde.fun = function(x){
    #  odds.p.tilde = exp(x.tilde%*%x)
    #  p.tlde1 = odds.p.tilde/(1+odds.p.tilde)
    #  c(samp$wt*c(samp[,y]-p.tlde1))%*%x.tilde
    #}
    #U.tilde_beta.tilde.fun(x=beta.tilde)
    #U.tilde_beta.tilde = jacobian(func=U.tilde_beta.tilde.fun, x = beta.tilde)
    #-t(c(samp$wt*p.tilde*(1-p.tilde))*x.tilde)%*%x.tilde
    
    Phi_theta = rbind(cbind(U_beta,       U_eta,       U_beta.tilde),
                      cbind(S_beta,       S_eta,       S_beta.tilde),
                      cbind(U.tilde_beta, U.tilde_eta, U.tilde_beta.tilde))
    US.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%U_eta%*%solve(S_eta)),
                   cbind(S_beta, solve(S_eta)))
    Phi_theta.inv = rbind(cbind(US.inv, -US.inv%*%rbind(U_beta.tilde, S_beta.tilde)%*%solve(U.tilde_beta.tilde)),
                          cbind(U.tilde_beta, U.tilde_eta, solve(U.tilde_beta.tilde))
    )
    #Phi_theta.inv1 = solve(Phi_theta)
    #t(solve(U.tilde_beta.tilde)%*%
    #    t(c(samp[,y]-p.tilde)*x.tilde))
    beta_w = t(Phi_theta.inv%*%t(Phi_w))
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))
  }
  if(structure=="subsample"){
    samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
    U_w.1 = matrix(0,nrow(samp), ncol(X))
    U_w.1[samp[,sub]==1,] = samp0$wt0*(u-v.mtx0%*%B)
    U_w.1 = v.mtx%*%B+U_w.1
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    #U_w.1 = c(samp[,y]-predict(mdl, samp, type="response"))*model.matrix(lm(as.formula(gsub(y, "wt", fit.y)),data=samp))
    #colSums(samp0$wt.c*c(samp0[,y]-p)*X)
    #colSums(samp$wt*c(samp[,y]-predict(mdl, samp, type="response"))*as.matrix(cbind(1,samp$x1,samp$x3,samp$x4,samp$x3*samp$x4)))
    #colSums((samp$wt*v.mtx)%*%B)+colSums(samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B))
    #e = matrix(0,nrow(samp), ncol(X))
    #e[samp[,sub]==1,] = samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B)
    #var((samp$wt*v.mtx)%*%B+e)
    #S_w.1 = matrix(0,nrow(v.mtx), ncol(v.mtx))
    #S_w.1 = c(v.mtx%*%eta)*v.mtx
    #Phi_w.1 = cbind(U_w.1, S_w.1)
    #TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, Phi_w.1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD))
    #var.Phi.2 = t(c((samp0$wt0-1)*samp$wt[samp[,sub]==1]*samp0$wt)*cbind(u,v.mtx0))%*%cbind(u,v.mtx0)
    #var.all = diag(Phi_theta.inv%*%(var.Phi.1+var.Phi.2)%*%t(Phi_theta.inv))
    var.all = diag(solve(U_beta)%*%(var.Phi.1)%*%t(solve(U_beta)))
  }
  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              beta.tilde = beta.tilde, beta.tilde.var = vcov(imp.mdl),
              eta = eta,
              var.all = var.all))
}

################################################################################################
# ds: survey design of the full sample
# fit: outcome model
# v.mtx:matrix of auxiliary variables in the full sample
# sub: variable name of the indicator for subsample
# sub.wt: weight variable name for the sub sample
# structure: data structure "combine" for combining multiple cycles; "subsample" for a random sub sample 
beta.est.varPop = function(ds, fit, v.mtx0, X.pop){
  #ds = ds.all; fit = fit.y; v.mtx =v.mtx.i;  
  #sub = "cycle"; method="linear", sub.wt = "wt1", structure="comb"
  # design variables
  samp0 = ds$variables
  samp0$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp0$strata = as.numeric(unlist(ds$strata))
  samp0$wt = as.numeric(unlist(1/ds$allprob))
  rm(ds)
  # outcome variable
  fit = as.formula(fit); 
  y = all.vars(fit)[1]
  # sample estimate of totals
  X.hat  = c(c(samp0$wt)%*%v.mtx0)
  # GREG adjustment factor
  greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=X.pop#, Vs0.hat=X.hat
                    )
  f = greg.out$f
  f[f<0]=0
  eta = greg.out$eta

  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  U_eta = t(samp0$wt*u)%*%v.mtx0
  S_beta = matrix(0,ncol(v.mtx0),length(beta.est))
  S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
  Phi_theta = rbind(cbind(U_beta, U_eta),
                    cbind(S_beta, S_eta))
  Phi_theta.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%t(B)),
                        cbind(S_beta, solve(S_eta)))
  U_w = samp0$wt.c*u
  S_w = samp0$wt.c*v.mtx0
  Phi_w = cbind(U_w, S_w)
  beta_w = -t(Phi_theta.inv%*%t(Phi_w))
  TD.dat = as.data.frame(cbind(1, samp0$psu, samp0$strata, beta_w))
  names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
  ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
  var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))

  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              eta = eta,
              var.all = var.all))
}


beta.est.var0 = function(ds, fit, v.mtx, sub, sub.wt, structure){
  #ds = ds.all; fit = fit.y; v.mtx =v.mtx.i;  
  #sub = "cycle"; method="linear", sub.wt = "wt1", structure="comb"
  # design variables
  samp = ds$variables
  samp$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt = as.numeric(unlist(1/ds$allprob))
  samp0 = samp[samp[,sub]==1,]
  samp0$wt = samp0[,sub.wt]
  rm(ds)
  # outcome variable
  fit = as.formula(fit); 
  y = all.vars(fit)[1]
  v.mtx0 = v.mtx[samp[,sub]==1,]
  # subsample estimate of influence total 
  #Vs0.hat = c(c(samp0$wt)%*%v.mtx0)
  # combined sample estimate of influence total
  VS.hat  = c(c(samp$wt)%*%v.mtx)
  # GREG adjustment factor
  #if(method=="linear"){
  greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=VS.hat#, Vs0.hat=Vs0.hat
                    )
  f = greg.out$f
  f[f<0]=0
  eta = greg.out$eta
  #c_eta = v.mtx0
  #c_v = eta
  #}
  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  if(structure=="combine"){
    #partial derivative of the estimating equation system, w.r.t. the parameters
    #U_beta.fun = function(x){
    #  odds.p = exp(X%*%x)
    #  p1 = odds.p/(1+odds.p)
    #  c(samp0$wt.c*c(samp0[,y]-p1))%*%X
    #}
    #U_beta.fun(x=beta.est)
    #U_beta = jacobian(func=U_beta.fun, x = beta.est)
    #U_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    #}
    #U_eta.fun(x=eta)
    #U_eta = jacobian(func=U_eta.fun, x = eta)
    U_eta = t(samp0$wt*u)%*%v.mtx0
    S_beta = matrix(0,length(beta.est),length(beta.est))
    #S_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  c(f1*samp0$wt)%*%v.mtx0-VS.hat
    #}
    #S_eta.fun(x=eta)
    #S_eta = jacobian(func=S_eta.fun, x = eta)
    S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
    Phi_theta = rbind(cbind(U_beta, U_eta),
                      cbind(S_beta, S_eta))
    #t(solve(U.tilde_beta.tilde)%*%
    #    t(c(samp[,y]-p.tilde)*x.tilde))
    Phi_theta.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%t(B)),
                          cbind(S_beta, solve(S_eta)))
    #Phi_theta.inv = solve(Phi_theta)
    #round(Phi_theta.inv,5)==round(solve(Phi_theta),5)
    #partial derivative of the estimating equation system, w.r.t. the sample weights
    U_w = matrix(0,nrow(samp), ncol(X))
    U_w[samp[,sub]==1,] = samp0$wt.c*u
    S1_w = matrix(0,nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[,sub]==1,] = samp0$wt.c*v.mtx0
    S_w  = S1_w-samp$wt*v.mtx
    Phi_w = cbind(U_w, S_w)
    beta_w = -t(Phi_theta.inv%*%t(Phi_w))
    #eta_w = -t(solve(S_eta)%*%t(S_w))
    #beta_w.c = -t(solve(U_beta)%*%t(u))
    #U.dat = as.data.frame(cbind(samp0$wt.c, samp0$psu, samp0$strata, beta_w.c))
    #names(U.dat)=c("wt.c", paste0("V", 1:(ncol(U.dat)-1)))
    #ds.U = svydesign(ids=~psu, strata =~strat, weights=~wt.c, data=U.dat)
    #var.beta = diag(vcov(svytotal(as.formula(paste0("~", paste(names(U.dat)[-1], collapse="+"))), ds.U)))
    #diag(vcov(mdl))
    #S.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, eta_w))
    #names(S.dat)=c("wt", "cycle", paste0("V", 1:(ncol(S.dat)-2)))
    #ds.S = svydesign(ids=~psu, strata =~strata, weights=~wt, data=S.dat)
    #var.eta = diag(vcov(svytotal(as.formula(paste0("~", paste(names(S.dat)[-c(1,2)], collapse="+"))), ds.S)))
    
    #Phi.dat = as.data.frame(cbind(samp$wt, samp[,cyc], Phi_w))
    #names(Phi.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(Phi.dat)-3)))
    #ds.Phi = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=Phi.dat)
    #var.Phi = vcov(svytotal(as.formula(paste0("~", paste(names(Phi.dat)[-c(1,2)], collapse="+"))), ds.Phi))
    #diag((Phi_theta.inv%*%var.Phi)%*%t(Phi_theta.inv))
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))
  }
  if(structure=="subsample"){
    samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
    U_w.1 = matrix(0,nrow(samp), ncol(X))
    U_w.1[samp[,sub]==1,] = samp0$wt0*(u-v.mtx0%*%B)
    U_w.1 = v.mtx%*%B+U_w.1
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    #U_w.1 = c(samp[,y]-predict(mdl, samp, type="response"))*model.matrix(lm(as.formula(gsub(y, "wt", fit.y)),data=samp))
    #colSums(samp0$wt.c*c(samp0[,y]-p)*X)
    #colSums(samp$wt*c(samp[,y]-predict(mdl, samp, type="response"))*as.matrix(cbind(1,samp$x1,samp$x3,samp$x4,samp$x3*samp$x4)))
    #colSums((samp$wt*v.mtx)%*%B)+colSums(samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B))
    #e = matrix(0,nrow(samp), ncol(X))
    #e[samp[,sub]==1,] = samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B)
    #var((samp$wt*v.mtx)%*%B+e)
    #S_w.1 = matrix(0,nrow(v.mtx), ncol(v.mtx))
    #S_w.1 = c(v.mtx%*%eta)*v.mtx
    #Phi_w.1 = cbind(U_w.1, S_w.1)
    #TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, Phi_w.1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD))
    #var.Phi.2 = t(c((samp0$wt0-1)*samp$wt[samp[,sub]==1]*samp0$wt)*cbind(u,v.mtx0))%*%cbind(u,v.mtx0)
    #var.all = diag(Phi_theta.inv%*%(var.Phi.1+var.Phi.2)%*%t(Phi_theta.inv))
    var.all = diag(solve(U_beta)%*%(var.Phi.1)%*%t(solve(U_beta)))
  }
  #c(var.beta, var.eta)
  #var.all[1:5]
  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              eta = eta,
              var.all = var.all))
}

greg_f = function(wt0, v.mtx0, VS.hat#, Vs0.hat=NULL
                  ){
  #wt0 = samp0$wt
  vWv_inv = solve(t(wt0*v.mtx0)%*%v.mtx0)
  #if(is.null(Vs0.hat)) 
  Vs0.hat=c(c(wt0)%*%v.mtx0)
  eta = c(c(VS.hat-Vs0.hat)%*%vWv_inv)
  f = c(t(1+eta%*%t(v.mtx0)))
  return(list(eta = eta, f = f))    
}



raking_f = function(wt0, v.mtx0, VS.hat, N.S, Vs0.hat){
  wt0 = samp0$wt
  N.S=sum(samp$wt)
  dat = as.data.frame(cbind(wt=wt0, v.mtx0))
  names(dat)[-1] = paste0("V", 1:(ncol(dat)-1))
  ds = svydesign(ids=~1, data=dat, weights=~wt)
  names(VS.hat)=names(dat)[-1]
  f=weights(calibrate(ds, as.formula(paste0("~", paste(names(dat)[-1], collapse = "+"))), 
            c("(Intercept)"=N.S, VS.hat), calfun="linear"))/wt0
  solve(t(v.mtx0)%*%v.mtx0, t(v.mtx0)%*%f)
  
}

#weights(calibrate(ds.wt, as.formula(paste0("~", paste0(names(aux.tot)[-c(1,length(aux.tot))],collapse="+"))), 
#                  aux.tot[-c(length(aux.tot))], calfun="raking"))


#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps<-function(popul,Msize, n){
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
    sam.pps.data=as.data.frame(popul[pps.samID,])
    names(sam.pps.data) = names(popul)
  }else{sam.pps.data=popul[pps.samID,]}
  sam.pps.data$wt=sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}

#popul = s1; sel.p="p.delta"
pps.wor = function(popul, sel.p, n){
  N = nrow(popul)
  popul$id = c(1:N)
  if(sum(popul[,sel.p]*n>=1)>0){
    n1 = n1.tmp = sum((popul[,sel.p]*n)>=1)
    samp.indx = rep(0,N)
    samp.indx[popul$id[(popul[,sel.p]*n)>=1]]=1
    while(n1.tmp>0){
      popul_sub = popul[samp.indx==0,]
      popul_sub[,sel.p] = popul_sub[,sel.p]/sum(popul_sub[,sel.p])
      samp.indx[popul_sub$id[(popul_sub[,sel.p]*(n-n1))>=1]]=1
      n1.tmp=sum(popul_sub[,sel.p]*(n-n1)>=1)
      n1 = n1+n1.tmp
    }
    samp = popul[samp.indx==1,]
    samp$wt=1
    if(n1<n){
      samp = rbind(samp, sam.pps(popul = popul_sub, Msize = popul_sub[,sel.p], n = (n-n1)))
    }
  }else{
    samp = sam.pps(popul = popul, Msize = popul[,sel.p], n = n)
  }
  samp
}



###################################################################################################################
#Derivative of U and beta (cox regression) w.r.t. the weight, used for variance calculation                       #
# INPUT:                                                                                                          #
# surv.fit  - output of coxph or      #
# dat     - a dataframe including variables of sample weight (wt), time (t), event status (d)                     #
# wt      - name of the sample weight variable                                                                    #
# OUTPUT:                                                                                                         #
# beta_wt       - influence functions for beta coefficients (derivative of beta w.r.t. the sample weights)        #
# Delta_beta.wt - sample weighted influence functions for beta coefficients                                       #
###################################################################################################################
beta_wt.cox = function(surv.fit, calib.list=NULL){
  dat = eval(surv.fit$call$data)
  if(is.null(dat)){
    dat = eval(surv.fit$survey.design$variables)
    dat$wt = 1/surv.fit$survey.design$prob
  } else(dat$wt=1)  
  x.mtrx = model.matrix(surv.fit)
  rel_hzd = exp(x.mtrx%*%surv.fit$coefficients)
  fit.vars = all.vars(surv.fit$formula)
  t = fit.vars[1]; d = fit.vars[2]
  n = nrow(dat)
  p = ncol(x.mtrx)
  dat$id = 1:n
  dat$rel_hzd = rel_hzd
  dat$wt_e = dat$wt*dat$rel_hzd
  
  x.mtrx = as.matrix(x.mtrx[order(dat[,t]),])
  dat = dat[order(dat[,t]),]
  dat$H_dnom = rev(cumsum(rev(dat$wt_e)))# S0
  dat$H_num = as.matrix(cumsum(as.data.frame(c(dat$wt_e)*x.mtrx)[n:1,]))[n:1,] #S1
  
  ties = duplicated(dat[,t])
  if(sum(ties)>0){
    H_uniq = dat[!ties,c(t, "H_dnom", "H_num")]
    H_dat = H_uniq[rep(1:nrow(H_uniq), table(dat[,t])),]
    #h_dat[,t]==dat[,t]
    dat[, c("H_dnom", "H_num")] = H_dat[, c("H_dnom", "H_num")]
    
  }
  
  x.mtrx = as.matrix(x.mtrx[order(dat$id), ])
  dat = dat[order(dat$id), ]
  
  H = as.matrix(dat$H_num/dat$H_dnom)
  
  d_indx = which(dat[,d]==1)
  dat1 = dat[d_indx,]
  # create a two-column matrix sum_wt.t, with the first column being the unique event time, 
  # the second column being the sum of weights for individuals having the event time t
  if(sum(duplicated(dat1[, t]))>0){
    dat1$weight = dat1$wt
    
    dat1_tb <- data.table(dat1, key=t)
    names(dat1_tb)[ncol(dat1_tb)]="weight"
    names(dat1_tb)[names(dat1_tb)==t]="t_tb"
    #sum_wt.t <- as.data.frame(dat1_tb[, sum(weight), by=list(t_tb)])
    # Remove build warnings from above line
    weight   <- dat1_tb$weight
    sum_wt.t <- as.data.frame(dat1_tb[, sum(weight), by="t_tb"])
    
    
    #sum_wt.t = aggregate(dat1$wt*, by=list(rep(1, nrow(dat1)),dat1[,t]), sum)
    sum_wt.t = sum_wt.t[match(unique(dat1[, t]), sum_wt.t[,1]),]
    unq_t.d_indx = d_indx[!duplicated(dat1[, t])]
  }else{
    sum_wt.t=cbind(dat1[, t], dat1$wt)
    unq_t.d_indx = d_indx
  }
  
  U_w_2 = 0 # Ai 
  U_beta_1 = 0
  
  for(j in 1:nrow(sum_wt.t)){
    k=unq_t.d_indx[j]
    U_w_2.k = sum_wt.t[j, 2]*(c((dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx/dat$H_dnom[k]-
                                outer(c((dat[,t]>=dat[,t][k])*dat$rel_hzd),as.numeric(as.matrix(dat$H_num)[k,]))/(dat$H_dnom[k])^2)
    U_w_2 = U_w_2.k+U_w_2
    sum_j = t(as.matrix(x.mtrx*c(dat$wt_e*as.numeric(dat[,t]>=dat[,t][k]))))%*%as.matrix(x.mtrx)
    U_beta_1 =-sum_wt.t[j, 2]*sum_j/dat$H_dnom[k]+U_beta_1
  }
  
  Ui_wt= as.matrix(dat[,d]*(x.mtrx- H)-U_w_2)
  U_beta =  t(c(dat$wt[d_indx])*H[d_indx,])%*%H[d_indx,]+U_beta_1
  if(!is.null(calib.list)) {
    U_eta_1 = 0
    for(j in 1:nrow(sum_wt.t)){
      k=unq_t.d_indx[j]
      U_eta_1.k = sum_wt.t[j, 2]*(t(c(dat$wt/calib.list$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)*x.mtrx)%*%calib.list$A/dat$H_dnom[k]-
                                    c(as.matrix(dat$H_num)[k,])%o%c(c(dat$wt/calib.list$f*(dat[,t]>=dat[,t][k])*dat$rel_hzd)%*%calib.list$A)/(dat$H_dnom[k])^2)
      U_eta_1 = U_eta_1+U_eta_1.k
    }
    U_eta = t(as.matrix(c(dat$wt/calib.list$f*dat[,d])*(x.mtrx- H)))%*%calib.list$A - U_eta_1
    Ui_wt = c(calib.list$f)*Ui_wt
    beta_wt = -(t(U_eta%*%solve(calib.list$Q_eta)%*%t(c(calib.list$f)*calib.list$A)) + Ui_wt)%*%solve(U_beta)
    eta_wt =(c(calib.list$f)*calib.list$A)%*%solve(calib.list$Q_eta)
  }else{
    beta_wt = -Ui_wt%*%solve(U_beta)
    eta_wt = NULL
  }
  theta_wt = cbind(eta_wt, beta_wt)
  #print("beta_wt.cox OK")
  return(list(#U_beta_inv = solve(U_beta),
    beta_wt = beta_wt,
    ui = as.matrix(dat[,d]*(x.mtrx- H))
    #Delta_beta.wt = Delta_beta.wt,
    #U_eta = U_eta
  )
  )
}
