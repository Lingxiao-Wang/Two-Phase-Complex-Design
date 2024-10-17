#rm(list=ls())
if(!("biowulf"%in%ls())){
  #c=4
  fm.i=2
  NSIMU_tot=2
  NSIMU=2
  setwd("/Volumes/wangl29/Calibration_svy")
  seed = read.table("seed.txt", header = T)
  x1.mdl=T
  node.k=1
  path=paste0("/Volumes/wangl29/Calibration_svy/results_test/combine/x1",x1.mdl,"/", fm.i, "/")
}else{
  setwd("/home/wangl29/Calibration_svy")
}
library(survey)
library(MASS)
#library(writexl)
library(numDeriv)
library(ICC)
source("subfunctions.R")

N=2000000
M.psu = 10000
psu.size = N/M.psu; psu.size
n = 10000
m.psu = 100
n.i = n/m.psu;n.i
fc = m.psu/M.psu;fc*100
f.i = n.i/psu.size;f.i*100;
x_mu = c(0, 0, 0)
sd_x = c(1, 1, 1)
rho=c(0, .1, .5,.7, .9)
beta.pop=NULL
j=2 # Varying effect of log(x2)
c.all=c(2,4,6,8,10)
#for(j in 1:length(rho)){
  set.seed(8291.27)
  x_rho = matrix(c(1.0,    0  , rho[j],
                   0  ,    1.0, 0.0   ,
                   rho[j], 0.0, 1.0    ), 3, 3)
  x_sigma = x_rho*outer(sd_x, sd_x)
  x_mtrx = as.data.frame(mvrnorm(N, mu = x_mu, Sigma = x_sigma))
  x1 = x_mtrx[,1] # available in all cycles 
  x2 = x_mtrx[,2] # available in all cycles
  x5 = x_mtrx[,3] # available in all cycles
  x6 = rnorm(N) 
  #logx2 = log(abs(x2))
  x5_c=as.numeric(cut(x5,breaks=c(quantile(x5, probs = c(0,0.4,0.6,1))), include.lowest = T))
  cor(x5, x5_c)
  x3 =  as.numeric(runif(N)<0.3);mean(x3) # # independent minority group
  gamma1  =         c(0.0, 0.0, 0.08, 0.05, 0.25)
  odds.s1 = exp(cbind( x1,  x2,   x3,   x5,  x6)%*%gamma1) # informative weights? mis-specification? efficiency improvement
  
  gamma2 =          c(0.01, 0.05, 0.05)
  odds.s2 = exp(cbind( x1,    x2,   x3)%*%gamma2) # informative weights? mis-specification? efficiency improvement
  
  # x4 is the special measurement that is only available in one cycle
  a=c(0.5, 1.5)   # effect of log(x2), controlling the difference between the true and the fitted imputation model
  rlts=list()
  # true imputation model
  #x4 = 0.5*x2+a[j]*logx2+rnorm(N,mean=1,sd=0.5)
  x4 = 0.5*x2+a[2]*x5+rnorm(N,mean=1,sd=0.5)
  cor(x4, x1)
  cor(x4, x2)
  cor(x4, x5)
  #x4 = x2+0.5^x2^2+rnorm(N,mean=1,sd=0.5);cor(x2, x4)
  # impute x4
  x4.ii  = lm(x4~x2)$fitted.values
  #x4.i = lm(x4~x2+logx2)$fitted.values
  x4.i   = lm(x4~x2+x5)$fitted.values
  x4.iii = lm(x4~x2+x5_c)$fitted.values
  print(c(cor(x4, x4.i), cor(x4, x4.ii), cor(x4, x4.iii))) 
  # Generate the outcome 
  # population 
  pop = data.frame(id=1:N, x1, x2, x3, x4, x5, x5_c, x4.i, x4.ii, x4.iii, odds.s1, odds.s2)
  pop = pop[order(pop$odds.s1), ]
  pop$psu = rep(1:M.psu, each=psu.size)
  size.I = aggregate(pop$odds.s1, list(pop$psu), FUN=sum)[,2]
  # Generate the outcome
  if(!x1.mdl){
    beta = c(-3, 0.9, 0.5, 0.3)
    n_beta=length(beta)
    pop$y = get.y(betas = beta, design.x = cbind(1, pop$x3, pop$x4, pop$x3*pop$x4), ICCy = 0.1, 
                  M.psu = M.psu, psu.size = rep(psu.size, M.psu))
  }else{
    beta = c(-3, -0.7, 0.9, 0.5, 0.3)
    n_beta=length(beta)
    pop$y = get.y(betas = beta, design.x = cbind(1, pop$x1, pop$x3, pop$x4, pop$x3*pop$x4), ICCy = 0.1, 
                  M.psu = M.psu, psu.size = rep(psu.size, M.psu))
    
  }
  pop = pop[order(pop$id), ]
  y.mu = mean(pop$y)
  #ICCbare(pop$psu, pop$x1)
  #ICCbare(pop$psu, pop$x2)
  ##iccbin(pop$psu, pop$x3)
  #ICCbare(pop$psu, pop$x4)
  #ICCbare(pop$psu, pop$x5)
  #ICCbare(pop$psu, pop$x5_c)
  #ICCbare(pop$psu, pop$odds.s1)
  #ICCbare(pop$psu, pop$odds.s2)
  
  # Estimation methods
  method = c("naive.1", "wt.1", "calibX", 
             "calib.x4i", "calib.x4ii", "calib.x4iii", 
             "calib0.x4i", "calib0.x4ii", "calib0.x4iii",
             "naive.all", "naive.imp1", "naive.imp2", "naive.imp3", "wt.all", 
             "wt.imp1", "wt.imp2", "wt.imp3")
  length(method)
  # Fitted analysis model
  fit.y = c("y~x3*x4", "y~x1+x3*x4")[fm.i]
  beta.pop = rbind(beta.pop, glm(fit.y, family="binomial", pop)$coef); beta.pop
  n_beta.fit = ncol(beta.pop)
  #influence functions in FP
  imp.pop = glm(as.formula(gsub("x4", "x4.i", fit.y)), pop, family="binomial")
  p.tilde.pop = imp.pop$fitted.values
  x.tilde.pop = model.matrix(imp.pop)
  beta.tilde.pop.i = imp.pop$coefficients
  # influence functions
  U.pop_beta.pop=-t(c(p.tilde.pop*(1-p.tilde.pop))*x.tilde.pop)%*%x.tilde.pop
  v.mtx.pop.i = t(solve(U.pop_beta.pop)%*%
                    t(c(pop$y-p.tilde.pop)*x.tilde.pop))
  
  imp.pop = glm(as.formula(gsub("x4", "x4.ii", fit.y)), pop, family="binomial")
  p.tilde.pop = imp.pop$fitted.values
  x.tilde.pop = model.matrix(imp.pop)
  beta.tilde.pop.ii = imp.pop$coefficients
  # infleunce functions
  U.pop_beta.pop=-t(c(p.tilde.pop*(1-p.tilde.pop))*x.tilde.pop)%*%x.tilde.pop
  v.mtx.pop.ii = t(solve(U.pop_beta.pop)%*%
                     t(c(pop$y-p.tilde.pop)*x.tilde.pop))
  
  imp.pop = glm(as.formula(gsub("x4", "x4.i", fit.y)), pop, family="binomial")
  p.tilde.pop = imp.pop$fitted.values
  x.tilde.pop = model.matrix(imp.pop)
  beta.tilde.pop.i = imp.pop$coefficients
  # influence functions
  U.pop_beta.pop=-t(c(p.tilde.pop*(1-p.tilde.pop))*x.tilde.pop)%*%x.tilde.pop
  v.mtx.pop.i = t(solve(U.pop_beta.pop)%*%
                    t(c(pop$y-p.tilde.pop)*x.tilde.pop))
  
  imp.pop = glm(as.formula(gsub("x4", "x4.iii", fit.y)), pop, family="binomial")
  p.tilde.pop = imp.pop$fitted.values
  x.tilde.pop = model.matrix(imp.pop)
  beta.tilde.pop.iii = imp.pop$coefficients
  # infleunce functions
  U.pop_beta.pop=-t(c(p.tilde.pop*(1-p.tilde.pop))*x.tilde.pop)%*%x.tilde.pop
  v.mtx.pop.iii = t(solve(U.pop_beta.pop)%*%
                      t(c(pop$y-p.tilde.pop)*x.tilde.pop))
  beta.pop = rbind(beta.pop, beta.tilde.pop.i, beta.tilde.pop.ii, beta.tilde.pop.iii)

  for(c.j in 1:length(c.all)){
    c=c.all[c.j]
    n1 = round(n/c)
  beta.est       = array(0, c(NSIMU, n_beta.fit , length(method)))
  beta.var       = array(0, c(NSIMU, n_beta.fit , length(method)))
  #beta.tilde.est = array(0, c(NSIMU, n_beta.fit , 3))
  #beta.tilde.var = array(0, c(NSIMU, n_beta.fit , 3))
  #eta.est        = array(0, c(NSIMU, n_beta.fit , 6))
  #eta.var        = array(0, c(NSIMU, n_beta.fit , 6))
  print(date())
  simu=1
  for(simu in simu:NSIMU){
    set.seed(seed[simu,1])
    seed.sim = runif(c)*1e5
    samp.all = NULL
    i=1
    for (i in i:c){
      samp.i = samp.slct(seed     = seed.sim[i], 
                         fnt.pop  = pop, 
                         n        = n1, 
                         psu.name = "psu", 
                         m.psu    = m.psu, 
                         dsgn     = "pps-pps", 
                         size.I   = size.I, 
                         size     = "odds.s2")
      samp.i$cycle=i
      samp.all = rbind(samp.all, samp.i)
      print(i)
    }
    samp.all$wt.1 = samp.all$wt
    samp.all$wt = samp.all$wt/c
    v.mtx.i   = v.mtx.pop.i  [samp.all$id,]
    v.mtx.ii  = v.mtx.pop.ii [samp.all$id,]
    v.mtx.iii = v.mtx.pop.iii[samp.all$id,]
    samp.1 = samp.all[samp.all$cycle==1,]
    ds.1   = svydesign(ids=~psu, strata = NULL, weights=~wt.1, data=samp.1)
    ds.all = svydesign(ids=~psu, strata = ~cycle, weights=~wt, nest=T, data=samp.all)
    #naive.1
    beta.est[simu,,1] = glm(fit.y, samp.all[samp.all$cycle==1,], family="binomial")$coeff
    beta.var[simu,,1] = diag(vcov(glm(fit.y, samp.all[samp.all$cycle==1,], family="binomial")))
    #wt.1
    beta.est[simu,,2] = svyglm(fit.y, ds.1, family="binomial")$coeff
    beta.var[simu,,2] = diag(vcov(svyglm(fit.y, ds.1, family="binomial")))
    #naive.all
    beta.est[simu,,10] = glm(fit.y, samp.all, family="binomial")$coeff
    beta.var[simu,,10] = diag(vcov(glm(fit.y, samp.all, family="binomial")))
    #naive.imp1
    beta.est[simu,,11] = glm(as.formula(gsub("x4", "x4.i", fit.y)), samp.all, family="binomial")$coeff
    beta.var[simu,,11] = diag(vcov(glm(as.formula(gsub("x4", "x4.i", fit.y)), samp.all, family="binomial")))
    #naive.imp2
    beta.est[simu,,12] = glm(as.formula(gsub("x4", "x4.ii", fit.y)), samp.all, family="binomial")$coeff
    beta.var[simu,,12] = diag(vcov(glm(as.formula(gsub("x4", "x4.ii", fit.y)), samp.all, family="binomial")))
    #naive.imp2
    beta.est[simu,,13] = glm(as.formula(gsub("x4", "x4.ii", fit.y)), samp.all, family="binomial")$coeff
    beta.var[simu,,13] = diag(vcov(glm(as.formula(gsub("x4", "x4.ii", fit.y)), samp.all, family="binomial")))
    #wt.all
    beta.est[simu,,14] = svyglm(fit.y, ds.all, family="binomial")$coeff
    beta.var[simu,,14] = diag(vcov(svyglm(fit.y, ds.all, family="binomial")))
    #wt.imp1
    imp.mdl1 = svyglm(as.formula(gsub("x4", "x4.i", fit.y)), ds.all, family="binomial")
    beta.est[simu,,15] = imp.mdl1$coeff
    beta.var[simu,,15] = diag(vcov(imp.mdl1))
    #wt.imp2
    imp.mdl2 = svyglm(as.formula(gsub("x4", "x4.ii", fit.y)), ds.all, family="binomial")
    beta.est[simu,,16] = imp.mdl2$coeff
    beta.var[simu,,16] = diag(vcov(imp.mdl2))
    #wt.imp3
    imp.mdl3 = svyglm(as.formula(gsub("x4", "x4.iii", fit.y)), ds.all, family="binomial")
    beta.est[simu,,17] = imp.mdl3$coeff
    beta.var[simu,,17] = diag(vcov(imp.mdl3))
    # calibration 
    
    # calibration on x
    calibX = beta.est.varPop(ds = ds.1, fit = fit.y, 
                             v.mtx0 = as.matrix(samp.1[,c(2:4, 6)]), 
                             X.pop = c(colSums(pop[,c(2:4, 6)])))
    beta.est[simu,,3] = calibX$beta.est; beta.var[simu,,3] = calibX$var.all[1:n_beta.fit]
    
    calib.x4i  = beta.est.var(ds = ds.all, fit = fit.y, fit.i =gsub("x4", "x4.i", fit.y), 
                              sub = "cycle", sub.wt = "wt.1", structure="combine")
    calib.x4ii = beta.est.var(ds = ds.all, fit = fit.y, fit.i =gsub("x4", "x4.ii", fit.y), 
                              sub = "cycle", sub.wt = "wt.1", structure="combine")
    calib.x4iii= beta.est.var(ds = ds.all, fit = fit.y, fit.i =gsub("x4", "x4.iii", fit.y), 
                              sub = "cycle", sub.wt = "wt.1", structure="combine")
    beta.est[simu,,4] = calib.x4i  $beta.est; beta.var[simu,,4] = calib.x4i  $var.all[1:n_beta.fit]
    beta.est[simu,,5] = calib.x4ii $beta.est; beta.var[simu,,5] = calib.x4ii $var.all[1:n_beta.fit]
    beta.est[simu,,6] = calib.x4iii$beta.est; beta.var[simu,,6] = calib.x4iii$var.all[1:n_beta.fit]

    #eta.est [simu,,1] = calib.x4i  $eta     ; eta.var [simu,,1] = calib.x4i  $var.all[1:n_beta.fit+n_beta.fit]
    #eta.est [simu,,2] = calib.x4ii $eta     ; eta.var [simu,,2] = calib.x4ii $var.all[1:n_beta.fit+n_beta.fit]
    #eta.est [simu,,3] = calib.x4iii$eta     ; eta.var [simu,,3] = calib.x4iii$var.all[1:n_beta.fit+n_beta.fit]
    #
    #beta.tilde.est[simu,,1]  = calib.x4i  $beta.tilde; beta.tilde.var[simu,,1]  = calib.x4i  $var.all[1:n_beta.fit+n_beta.fit*2]
    #beta.tilde.est[simu,,2]  = calib.x4ii $beta.tilde; beta.tilde.var[simu,,2]  = calib.x4ii $var.all[1:n_beta.fit+n_beta.fit*2]
    #beta.tilde.est[simu,,3]  = calib.x4iii$beta.tilde; beta.tilde.var[simu,,3]  = calib.x4iii$var.all[1:n_beta.fit+n_beta.fit*2]

    calib0.x4i  = beta.est.var0(ds = ds.all, fit = fit.y, v.mtx =v.mtx.i,  
                                sub = "cycle", sub.wt = "wt.1", structure="combine")
    calib0.x4ii = beta.est.var0(ds = ds.all, fit = fit.y, v.mtx =v.mtx.ii, 
                                sub = "cycle", sub.wt = "wt.1", structure="combine")
    calib0.x4iii= beta.est.var0(ds = ds.all, fit = fit.y, v.mtx =v.mtx.iii, 
                                sub = "cycle", sub.wt = "wt.1", structure="combine")
    beta.est[simu,,7] = calib0.x4i  $beta.est; beta.var[simu,,7] = calib0.x4i  $var.all[1:n_beta.fit]
    beta.est[simu,,8] = calib0.x4ii $beta.est; beta.var[simu,,8] = calib0.x4ii $var.all[1:n_beta.fit]
    beta.est[simu,,9] = calib0.x4iii$beta.est; beta.var[simu,,9] = calib0.x4iii$var.all[1:n_beta.fit]
    #eta.est [simu,,4] = calib0.x4i  $eta     ; eta.var [simu,,4] = calib0.x4i  $var.all[1:n_beta.fit+n_beta.fit]
    #eta.est [simu,,5] = calib0.x4ii $eta     ; eta.var [simu,,5] = calib0.x4ii $var.all[1:n_beta.fit+n_beta.fit]
    #eta.est [simu,,6] = calib0.x4iii$eta     ; eta.var [simu,,6] = calib0.x4iii$var.all[1:n_beta.fit+n_beta.fit]
    print(paste0("simu", simu))
  }
  
  
  beta_est.mtx = cbind(simu_id = node.k, matrix(beta.est, NSIMU, n_beta.fit *length(method)))
  beta_var.mtx = cbind(simu_id = node.k, matrix(beta.var, NSIMU, n_beta.fit *length(method)))
  #beta.tilde_est.mtx = cbind(simu_id = node.k, matrix(beta.tilde.est, NSIMU, n_beta.fit *3))
  #beta.tilde_var.mtx = cbind(simu_id = node.k, matrix(beta.tilde.var, NSIMU, n_beta.fit *3))
  #eta_est.mtx = cbind(simu_id = node.k, matrix(eta.est, NSIMU, n_beta.fit *6))
  #eta_var.mtx = cbind(simu_id = node.k, matrix(eta.var, NSIMU, n_beta.fit *6))
  
  write.table(beta_est.mtx,       paste0(path, "beta.rho"          ,rho[j],"_","c",c.all[c.j],"_", node.k,".txt"), sep=",")
  write.table(beta_var.mtx,       paste0(path, "beta_var.rho"      ,rho[j],"_","c",c.all[c.j],"_", node.k,".txt"), sep=",")
  #write.table(beta.tilde_est.mtx, paste0(path, "beta.tilde.a"    ,a[j],"_",node.k,".txt"), sep=",")
  #write.table(beta.tilde_var.mtx, paste0(path, "beta.tilde_var.a",a[j],"_",node.k,".txt"), sep=",")
  #write.table(eta_est.mtx,        paste0(path, "eta.a"           ,a[j],"_",node.k,".txt"), sep=",")
  #write.table(eta_var.mtx,        paste0(path, "eta_var.a"       ,a[j],"_",node.k,".txt"), sep=",")
  print(paste0("c=",c, ";n1=", nrow(samp.1)))
}
#is.list(rlts)
#rlts
#write_xlsx(rlts,paste0("results", fm.i,".xlsx"))

pop_param = list(NSIMU_tot  = NSIMU_tot,
                 NSIMU      = NSIMU,
                 beta       = beta,
                 beta.pop   = beta.pop,
                 method     = method,
                 rho        = rho,
                 c.all      = c.all)
