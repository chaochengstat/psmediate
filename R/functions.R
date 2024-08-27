################all functions
# library(MASS)
# library(nleqslv)
# library(boot)
# library(caret)
# library(SuperLearner)
# library(ranger)


#' Expit function
#'
#' @param x A number vector.
#' @returns A numeric vector.
#' @examples
#' expit(0.1)
expit <- function(x){exp(x)/(1+exp(x))}


#' Expit function
#'
#' @param x A number vector.
#' @returns A numeric vector.
#' @examples
#' logit(0.5)
logit <- function(x){log(x/(1-x))}




#' Generate simulated dataset
#'
#' @param n sample size.
#' @param noncompliance "one-sided" or "two-sided" noncompliance type
#' @returns A simulated data with sample size n.
#' @examples
#' head(data_gen(n=100))
data_gen=function(n,noncompliance="two-sided") {
  # baseline covariates
  U <- MASS::mvrnorm(n,mu=cbind(0,0,0,0), Sigma=diag(c(1,1,1,1)))
  colnames(U)=c("U1","U2","U3","U4")
  # treatment
  pz <- expit(U %*% c(-1,0.5,-0.25,-0.1))
  Z <- rbinom(n,1,prob = pz)
  ps <- expit(-1 + 2*Z + c(U %*% c(1,-0.8,0.6,-1)))
  D <- rbinom(n,1,prob = ps)
  if (noncompliance=="one-sided") {
    D[which(Z==0 & S==1)] = 0
  }
  pm <- expit(-1.8 + 2*Z + 1.5*D + c(U %*% c(0.5,0.5,0.5,0.5)))
  M <- rbinom(n,1,prob = pm)
  eY = rnorm(n)
  Y = 10 + 3* Z - 2*D + 2*M + c(U %*% c(0.5,0.5,0.5,0.5)) + eY
  X = data.frame(X1=U[,1],X2=U[,2],X3=U[,3],X4=U[,4])
  dat=as.data.frame(cbind(Z,D,M,Y,X))

  return(dat)
}





#' Singly and multiply robust estimators for calculating the generalized mediation functional theta_(d1d0)^(zz')
#'
#' @param Xpi design matrix for the treatment probability model
#' @param Xe design matrix for the principal score model
#' @param Xm design matrix for the mediator model
#' @param Xo design matrix for the outcome model
#' @param Z  vector for the treatment
#' @param D vector for the post-treatment event
#' @param M vector for the mediator
#' @param Y vector for the outcome
#' @param group principal stratum ("compliers", "always takers", or "never takers")
#' @param type superscript for the mediation functional (i.e., zz'). The default value is "10".
#' @param noncompliance noncompliance type (one-sided or two-sided corresponding to strong and standard monotonicity assumption)
#'
#' @returns estimates for the generalized mediation functional based on the singly robust estimator as
#' well as the multiply robust estimator
GMF_par = function(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="10",noncompliance="two-sided") {
  # sample size
  n=length(Y)
  # outcome type
  if (range(Y)[1]==0 & range(Y)[2]==1) {
    type_y = "binary"
  } else {
    type_y = "continuous"
  }
  # identify d1 and d0
  if (group == "compliers") {d1=1;d0=0}
  if (group == "always takers") {d1=1;d0=1}
  if (group == "never takers") {d1=0;d0=0}
  # identify z and z'
  z  = as.integer(strsplit(type,split="")[[1]][1])
  zp = as.integer(strsplit(type,split="")[[1]][2])
  # identify dz and dz'
  dz = ifelse(z,d1,d0)
  dzp = ifelse(zp,d1,d0)
  # identify one-sided noncompliance and two-sided noncompliance
  #if (sum(D[Z==0])>=1) {noncompliance="two-sided"}
  #if (sum(D[Z==0])==0) {noncompliance="one-sided"}
  # identify k and z*d*
  k=abs(d1-d0)
  if (group == "compliers") {zs=1;ds=1}
  if (group == "always takers") {zs=0;ds=1}
  if (group == "never takers") {zs=1;ds=0}


  # propensity score model
  m_z=glm(Z~Xpi,family=binomial)
  pi_1x = predict(m_z,type = "response")
  epsilon = 0.01
  pi_1x[which(pi_1x<epsilon)]=epsilon;pi_1x[which(pi_1x>1-epsilon)]=1-epsilon
  pi_0x = 1-pi_1x
  if(z==1) {pi_zx=pi_1x}
  if(z==0) {pi_zx=pi_0x}
  if(zp==1) {pi_zpx=pi_1x}
  if(zp==0) {pi_zpx=pi_0x}
  if(zs==1) {pi_zsx=pi_1x}
  if(zs==0) {pi_zsx=pi_0x}


  # principal score model
  if (noncompliance=="two-sided") {
    m_p = glm(D~Z+Xe,family=binomial)
    Xe_mm = model.matrix(m_p)
    Xe_mm1 = Xe_mm0 = Xe_mm
    Xe_mm1[,"Z"]=1;Xe_mm0[,"Z"]=0
    p_11x = predict(m_p,as.data.frame(Xe_mm1),type = "response")
    p_01x = predict(m_p,as.data.frame(Xe_mm0),type = "response")
    p_11x[which(p_11x<epsilon)]=epsilon;p_11x[which(p_11x>1-epsilon)]=1-epsilon
    p_01x[which(p_01x<epsilon)]=epsilon;p_01x[which(p_01x>1-epsilon)]=1-epsilon
    p_00x = 1-p_01x
    p_10x = 1-p_11x
  } else {
    D_oneside = D[Z==1]
    Xe_oneside = Xe[Z==1,]
    m_p = glm(D_oneside~Xe_oneside,family=binomial)
    Xe_mm = model.matrix(D~Xe)
    p_11x = expit(c(Xe_mm %*% m_p$coefficients))
    p_01x = rep(0,n)
    p_11x[which(p_11x<epsilon)]=epsilon;p_11x[which(p_11x>1-epsilon)]=1-epsilon
    p_00x = 1-p_01x
    p_10x = 1-p_11x
  }
  if(z==1 & dz==1) {p_zdx=p_11x}
  if(z==0 & dz==1) {p_zdx=p_01x}
  if(z==1 & dz==0) {p_zdx=p_10x}
  if(z==0 & dz==0) {p_zdx=p_00x}
  if(zp==1 & dzp==1) {p_zpdpx=p_11x}
  if(zp==0 & dzp==1) {p_zpdpx=p_01x}
  if(zp==1 & dzp==0) {p_zpdpx=p_10x}
  if(zp==0 & dzp==0) {p_zpdpx=p_00x}
  if(zs==1 & ds==1) {p_zsdsx=p_11x}
  if(zs==0 & ds==1) {p_zsdsx=p_01x}
  if(zs==1 & ds==0) {p_zsdsx=p_10x}
  if(zs==0 & ds==0) {p_zsdsx=p_00x}

  # doubly robust estimator for p_zsds and p_01
  if (group == "compliers") {p_zsds = mean((Z==1)*((D==1)-p_11x)/pi_1x + p_11x)}
  if (group == "always takers") {p_zsds = mean((Z==0)*((D==1)-p_01x)/pi_0x + p_01x)}
  if (group == "never takers") {p_zsds = mean((Z==1)*((D==0)-p_10x)/pi_1x + p_10x)}
  if (noncompliance=="two-sided") {
    p_01 = mean((1-Z)*(D-p_01x)/pi_0x + p_01x)
  } else {
    p_01 = 0
  }
  # principal score weighting
  if (group == "compliers") {e_x = p_11x - k*p_01x}
  if (group == "always takers") {e_x = p_01x - k*p_01x}
  if (group == "never takers") {e_x = p_10x - k*p_01x}
  weight_ps = e_x/(p_zsds-k*p_01)

  # mediator model
  m_m = glm(M~Z+D+Xm,family=binomial)
  Xm_mm = model.matrix(m_m)
  Xm_mm11 = Xm_mm00 = Xm_mm10 = Xm_mm01 = Xm_mm
  Xm_mm10[,"Z"]=Xm_mm11[,"Z"]=1;Xm_mm00[,"Z"]=Xm_mm01[,"Z"]=0
  Xm_mm10[,"D"]=Xm_mm00[,"D"]=0;Xm_mm11[,"D"]=Xm_mm01[,"D"]=1
  r_001x = expit(c(Xm_mm00 %*% m_m$coefficients))
  r_001x[which(r_001x<epsilon)]=epsilon;r_001x[which(r_001x>1-epsilon)]=1-epsilon
  r_000x = 1-r_001x
  r_111x = expit(c(Xm_mm11 %*% m_m$coefficients))
  r_111x[which(r_111x<epsilon)]=epsilon;r_111x[which(r_111x>1-epsilon)]=1-epsilon
  r_110x = 1-r_111x
  r_101x = expit(c(Xm_mm10 %*% m_m$coefficients))
  r_101x[which(r_101x<epsilon)]=epsilon;r_101x[which(r_101x>1-epsilon)]=1-epsilon
  r_100x = 1-r_101x
  r_011x = expit(c(Xm_mm01 %*% m_m$coefficients))
  r_011x[which(r_011x<epsilon)]=epsilon;r_011x[which(r_011x>1-epsilon)]=1-epsilon
  r_010x = 1-r_011x
  if(z==1 & dz==1) {r_zd1x=r_111x;r_zd0x=r_110x}
  if(z==0 & dz==1) {r_zd1x=r_011x;r_zd0x=r_010x}
  if(z==1 & dz==0) {r_zd1x=r_101x;r_zd0x=r_100x}
  if(z==0 & dz==0) {r_zd1x=r_001x;r_zd0x=r_000x}
  if(zp==1 & dzp==1) {r_zpdp1x=r_111x;r_zpdp0x=r_110x}
  if(zp==0 & dzp==1) {r_zpdp1x=r_011x;r_zpdp0x=r_010x}
  if(zp==1 & dzp==0) {r_zpdp1x=r_101x;r_zpdp0x=r_100x}
  if(zp==0 & dzp==0) {r_zpdp1x=r_001x;r_zpdp0x=r_000x}


  # outcome model
  Y_sub = Y[Z==z & D==dz]
  Xo_sub =cbind(M,Xo)[which(Z==z & D==dz),]
  if (type_y=="continuous") {
    m_o = lm(Y_sub~Xo_sub)
    Xo_mm = model.matrix(Y~M+Xo)
    Xo_mm1 = Xo_mm0 = Xo_mm
    Xo_mm1[,"M"] = 1; Xo_mm0[,"M"] = 0
    mu_zd1x = c(Xo_mm1 %*% m_o$coefficients)
    mu_zd0x = c(Xo_mm0 %*% m_o$coefficients)
  } else {
    m_o = glm(Y_sub~Xo_sub,family=binomial)
    Xo_mm = model.matrix(Y~M+Xo)
    Xo_mm1 = Xo_mm0 = Xo_mm
    Xo_mm1[,"M"] = 1; Xo_mm0[,"M"] = 0
    mu_zd1x = expit(c(Xo_mm1 %*% m_o$coefficients))
    mu_zd0x = expit(c(Xo_mm0 %*% m_o$coefficients))
  }
  # moment-type estimators
  # M_{m+s+z} single robust estimator
  res1 = mean(weight_ps*(D==dz)*(Z==z)/(p_zdx*pi_zx)*(r_zpdp1x*M+r_zpdp0x*(1-M))/(r_zd1x*M+r_zd0x*(1-M))*Y)
  # M_{o+m+z} single robust estimator
  if (noncompliance=="two-sided") {
    res2 = mean(((D==ds)*(Z==zs)/(pi_zsx)-k*(1-Z)*D/pi_0x)*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))/(p_zsds-k*p_01)
  } else {
    res2 = mean(((D==ds)*(Z==zs)/(pi_zsx)-0)*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))/(p_zsds-k*p_01)
  }
  # M_{o+s+z} single robust estimator
  res3 = mean(weight_ps*(D==dzp)*(Z==zp)/(p_zpdpx*pi_zpx)*(mu_zd0x*(1-M)+mu_zd1x*M))
  # M_{o+m+s} single robust estimator
  res4 = mean(weight_ps*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))

  # multiply robust estimator
  if (noncompliance=="two-sided") {
    term1 = ((((D==ds)-p_zsdsx)*(Z==zs)/(pi_zsx)-k*(1-Z)*(D-p_01x)/pi_0x)*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))/(p_zsds-k*p_01)
  } else {
    term1 = ((((D==ds)-p_zsdsx)*(Z==zs)/(pi_zsx))*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))/(p_zsds-k*p_01)
  }
  term2 = weight_ps*(D==dz)*(Z==z)/(p_zdx*pi_zx)*(r_zpdp1x*M+r_zpdp0x*(1-M))/(r_zd1x*M+r_zd0x*(1-M))*(Y-mu_zd0x*(1-M)-mu_zd1x*M)
  term3 = weight_ps*(D==dzp)*(Z==zp)/(p_zpdpx*pi_zpx)*(mu_zd0x*(1-M)+mu_zd1x*M-(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))
  term4 = weight_ps*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x)
  res_mr = mean(term1+term2+term3+term4)
  out = c(res1,res2,res3,res4,res_mr)
  attr(out,"proportion") = p_zsds-k*p_01
  return(out)
}




#' Nonparametric efficient estimator for calculating the generalized mediation functional theta_(d1d0)^(zz')
#'
#' @param Xpi design matrix for the treatment probability model
#' @param Xe design matrix for the principal score model
#' @param Xm design matrix for the mediator model
#' @param Xo design matrix for the outcome model
#' @param Z  vector for the treatment
#' @param D vector for the post-treatment event
#' @param M vector for the mediator
#' @param Y vector for the outcome
#' @param group principal stratum ("compliers", "always takers", or "never takers")
#' @param type superscript for the mediation functional (i.e., zz'). The default value is "10".
#' @param noncompliance noncompliance type (one-sided or two-sided corresponding to strong and standard monotonicity assumption)
#' @param learners names of super learners (default GLM and random forest)
#' @param V number of folds of the cross-fitting procedure (default 5)
#' @param output.IF whether output the influence function (default false)
#'
#' @returns estimates for the generalized mediation functional based on the nonparametric efficient estimator
GMF_nonpar =function(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="10",
                     noncompliance="two-sided",
                     learners=c("SL.glm", "SL.ranger"),V=5,
                     output.IF=FALSE) {
  myseed=2024
  # sample size
  n=length(Y)
  # outcome type
  if (range(Y)[1]==0 & range(Y)[2]==1) {
    type_y = "binary"
  } else {
    type_y = "continuous"
  }
  # identify d1 and d0
  if (group == "compliers") {d1=1;d0=0}
  if (group == "always takers") {d1=1;d0=1}
  if (group == "never takers") {d1=0;d0=0}
  # identify z and z'
  z  = as.integer(strsplit(type,split="")[[1]][1])
  zp = as.integer(strsplit(type,split="")[[1]][2])
  # identify dz and dz'
  dz = ifelse(z,d1,d0)
  dzp = ifelse(zp,d1,d0)
  # identify one-sided noncompliance and two-sided noncompliance
  #if (sum(D[Z==0])>=1) {noncompliance="two-sided"}
  #if (sum(D[Z==0])==0) {noncompliance="one-sided"}
  # identify k and z*d*
  k=abs(d1-d0)
  if (group == "compliers") {zs=1;ds=1}
  if (group == "always takers") {zs=0;ds=1}
  if (group == "never takers") {zs=1;ds=0}

  IF1 = IF2=NULL # initialize influence function values
  set.seed(myseed)
  folds = createFolds(1:n, k = V)
  res=c()
  for (v in (1:V)) {
    #==========> Produce sample split <==========#
    Xpi_main=Xpi[-folds[[v]],]
    Xe_main=Xe[-folds[[v]],]
    Xm_main=Xm[-folds[[v]],]
    Xo_main=Xo[-folds[[v]],]
    Z_main=Z[-folds[[v]]]
    D_main=D[-folds[[v]]]
    M_main=M[-folds[[v]]]
    Y_main=Y[-folds[[v]]]

    Xpi_vali=Xpi[folds[[v]],]
    Xe_vali=Xe[folds[[v]],]
    Xm_vali=Xm[folds[[v]],]
    Xo_vali=Xo[folds[[v]],]
    Z_vali=Z[folds[[v]]]
    D_vali=D[folds[[v]]]
    M_vali=M[folds[[v]]]
    Y_vali=Y[folds[[v]]]


    # propensity score models
    m_z = SuperLearner::SuperLearner(
      Y          = Z_main,
      X          = as.data.frame(Xpi_main),
      family     = binomial(),
      SL.library = learners,
      control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
      cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
    )
    pi_1x    = predict(m_z, newdata = Xpi_vali, onlySL = TRUE)$pred[,1]
    epsilon = 0.0000001
    pi_1x[which(pi_1x<epsilon)]=epsilon;pi_1x[which(pi_1x>1-epsilon)]=1-epsilon
    pi_0x = 1-pi_1x
    if(z==1) {pi_zx=pi_1x}
    if(z==0) {pi_zx=pi_0x}
    if(zp==1) {pi_zpx=pi_1x}
    if(zp==0) {pi_zpx=pi_0x}
    if(zs==1) {pi_zsx=pi_1x}
    if(zs==0) {pi_zsx=pi_0x}


    # principal score models
    if (noncompliance=="two-sided") {
      m_p = SuperLearner::SuperLearner(
        Y          = D_main,
        X          = data.frame(Z=Z_main,Xe_main),
        family     = binomial(),
        SL.library = learners,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
        cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
      )
      Xe_mm1 = Xe_mm0 = data.frame(Z=Z_vali,Xe_vali)
      Xe_mm1[,"Z"]=1;Xe_mm0[,"Z"]=0
      p_11x = predict(m_p,newdata=Xe_mm1,onlySL = TRUE)$pred[,1]
      p_01x = predict(m_p,newdata=Xe_mm0,onlySL = TRUE)$pred[,1]
      p_11x[which(p_11x<epsilon)]=epsilon;p_11x[which(p_11x>1-epsilon)]=1-epsilon
      p_01x[which(p_01x<epsilon)]=epsilon;p_01x[which(p_01x>1-epsilon)]=1-epsilon
      p_00x = 1-p_01x
      p_10x = 1-p_11x
    } else {
      D_oneside = D_main[Z_main==1]
      Xe_oneside = Xe_main[Z_main==1,]
      m_p = SuperLearner::SuperLearner(
        Y          = D_oneside,
        X          = data.frame(Xe_oneside),
        family     = binomial(),
        SL.library = learners,
        control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
        cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
      )
      p_11x = predict(m_p,newdata=Xe_vali,onlySL = TRUE)$pred[,1]
      p_01x = rep(0,dim(Xe_vali)[1])
      p_11x[which(p_11x<epsilon)]=epsilon;p_11x[which(p_11x>1-epsilon)]=1-epsilon
      p_00x = 1-p_01x
      p_10x = 1-p_11x
    }
    if(z==1 & dz==1) {p_zdx=p_11x}
    if(z==0 & dz==1) {p_zdx=p_01x}
    if(z==1 & dz==0) {p_zdx=p_10x}
    if(z==0 & dz==0) {p_zdx=p_00x}
    if(zp==1 & dzp==1) {p_zpdpx=p_11x}
    if(zp==0 & dzp==1) {p_zpdpx=p_01x}
    if(zp==1 & dzp==0) {p_zpdpx=p_10x}
    if(zp==0 & dzp==0) {p_zpdpx=p_00x}
    if(zs==1 & ds==1) {p_zsdsx=p_11x}
    if(zs==0 & ds==1) {p_zsdsx=p_01x}
    if(zs==1 & ds==0) {p_zsdsx=p_10x}
    if(zs==0 & ds==0) {p_zsdsx=p_00x}

    # doubly robust estimator for p_zsds and p_01
    if (group == "compliers") {p_zsds = mean((Z_vali==1)*((D_vali==1)-p_11x)/pi_1x + p_11x)}
    if (group == "always takers") {p_zsds = mean((Z_vali==0)*((D_vali==1)-p_01x)/pi_0x + p_01x)}
    if (group == "never takers") {p_zsds = mean((Z_vali==1)*((D_vali==0)-p_10x)/pi_1x + p_10x)}
    if (noncompliance=="two-sided") {
      p_01 = mean((1-Z_vali)*(D_vali-p_01x)/pi_0x + p_01x)
    } else {
      p_01 = 0
    }
    # principal score weighting
    if (group == "compliers") {e_x = p_11x - k*p_01x}
    if (group == "always takers") {e_x = p_01x - k*p_01x}
    if (group == "never takers") {e_x = p_10x - k*p_01x}
    weight_ps = e_x/(p_zsds-k*p_01)

    # mediator model
    m_m = SuperLearner::SuperLearner(
      Y          = M_main,
      X          = data.frame(Z=Z_main,D=D_main,Xm_main),
      family     = binomial(),
      SL.library = learners,
      control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
      cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
    )
    Xm_mm = data.frame(Z=Z_vali,D=D_vali,Xm_vali)
    Xm_mm11 = Xm_mm00 = Xm_mm10 = Xm_mm01 = Xm_mm
    Xm_mm10[,"Z"]=Xm_mm11[,"Z"]=1;Xm_mm00[,"Z"]=Xm_mm01[,"Z"]=0
    Xm_mm10[,"D"]=Xm_mm00[,"D"]=0;Xm_mm11[,"D"]=Xm_mm01[,"D"]=1
    r_001x = predict(m_m,newdata=Xm_mm00,onlySL = TRUE)$pred[,1]
    r_001x[which(r_001x<epsilon)]=epsilon;r_001x[which(r_001x>1-epsilon)]=1-epsilon
    r_000x = 1-r_001x
    r_111x = predict(m_m,newdata=Xm_mm11,onlySL = TRUE)$pred[,1]
    r_111x[which(r_111x<epsilon)]=epsilon;r_111x[which(r_111x>1-epsilon)]=1-epsilon
    r_110x = 1-r_111x
    r_101x = predict(m_m,newdata=Xm_mm10,onlySL = TRUE)$pred[,1]
    r_101x[which(r_101x<epsilon)]=epsilon;r_101x[which(r_101x>1-epsilon)]=1-epsilon
    r_100x = 1-r_101x
    r_011x = predict(m_m,newdata=Xm_mm01,onlySL = TRUE)$pred[,1]
    r_011x[which(r_011x<epsilon)]=epsilon;r_011x[which(r_011x>1-epsilon)]=1-epsilon
    r_010x = 1-r_011x
    if(z==1 & dz==1) {r_zd1x=r_111x;r_zd0x=r_110x}
    if(z==0 & dz==1) {r_zd1x=r_011x;r_zd0x=r_010x}
    if(z==1 & dz==0) {r_zd1x=r_101x;r_zd0x=r_100x}
    if(z==0 & dz==0) {r_zd1x=r_001x;r_zd0x=r_000x}
    if(zp==1 & dzp==1) {r_zpdp1x=r_111x;r_zpdp0x=r_110x}
    if(zp==0 & dzp==1) {r_zpdp1x=r_011x;r_zpdp0x=r_010x}
    if(zp==1 & dzp==0) {r_zpdp1x=r_101x;r_zpdp0x=r_100x}
    if(zp==0 & dzp==0) {r_zpdp1x=r_001x;r_zpdp0x=r_000x}

    # outcome model
    Y_sub = Y_main[Z_main==z & D_main==dz]
    Xo_sub = data.frame(M=M_main,Xo_main)[which(Z_main==z & D_main==dz),]


    if(type_y=="continuous") {
      m_o = SuperLearner::SuperLearner(
        Y          = Y_sub,
        X          = Xo_sub,
        family     = gaussian(),
        SL.library = learners,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
      )
      Xo_mm = data.frame(M=M_vali,Xo_vali)
      Xo_mm1 = Xo_mm0 = Xo_mm
      Xo_mm1[,"M"] = 1; Xo_mm0[,"M"] = 0
      mu_zd1x = predict(m_o,newdata=Xo_mm1,onlySL = TRUE)$pred[,1]
      mu_zd0x = predict(m_o,newdata=Xo_mm0,onlySL = TRUE)$pred[,1]
    } else {
      m_o = SuperLearner::SuperLearner(
        Y          = Y_sub,
        X          = Xo_sub,
        family     = binomial(),
        SL.library = learners,
        control    = list(saveFitLibrary = TRUE),
        cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
      )
      Xo_mm = data.frame(M=M_vali,Xo_vali)
      Xo_mm1 = Xo_mm0 = Xo_mm
      Xo_mm1[,"M"] = 1; Xo_mm0[,"M"] = 0
      mu_zd1x = predict(m_o,newdata=Xo_mm1,onlySL = TRUE)$pred[,1]
      mu_zd0x = predict(m_o,newdata=Xo_mm0,onlySL = TRUE)$pred[,1]
    }


    # multiply robust estimator
    if (noncompliance=="two-sided") {
      term1 = ((((D_vali==ds)-p_zsdsx)*(Z_vali==zs)/(pi_zsx)-k*(1-Z_vali)*(D_vali-p_01x)/pi_0x)*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))
    } else {
      term1 = ((((D_vali==ds)-p_zsdsx)*(Z_vali==zs)/(pi_zsx))*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))
    }
    term2 = e_x*(D_vali==dz)*(Z_vali==z)/(p_zdx*pi_zx)*(r_zpdp1x*M_vali+r_zpdp0x*(1-M_vali))/(r_zd1x*M_vali+r_zd0x*(1-M_vali))*(Y_vali-mu_zd0x*(1-M_vali)-mu_zd1x*M_vali)
    term3 = e_x*(D_vali==dzp)*(Z_vali==zp)/(p_zpdpx*pi_zpx)*(mu_zd0x*(1-M_vali)+mu_zd1x*M_vali-(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x))
    term4 = e_x*(mu_zd0x*r_zpdp0x+mu_zd1x*r_zpdp1x)


    # doubly robust estimator for p_zsds and p_01
    if (group == "compliers") {term_psi1 = (Z_vali==1)*((D_vali==1)-p_11x)/pi_1x + p_11x}
    if (group == "always takers") {term_psi1 = (Z_vali==0)*((D_vali==1)-p_01x)/pi_0x + p_01x}
    if (group == "never takers") {term_psi1 = (Z_vali==1)*((D_vali==0)-p_10x)/pi_1x + p_10x}
    if (noncompliance=="two-sided") {
      term_psi2 = (1-Z_vali)*(D_vali-p_01x)/pi_0x + p_01x
    } else {
      term_psi2 = rep(0,length(Z_vali))
    }

    IF1=c(IF1,term1+term2+term3+term4)
    IF2=c(IF2,term_psi1-k*term_psi2)
    #res_mr = mean(term1+term2+term3+term4)/mean(p_11-p_01)
    #res=c(res,res_mr)
  }
  #print(mean(IF2))
  point = mean(IF1)/mean(IF2)
  SD = sqrt(var(1/mean(IF2)*(IF1-mean(IF1))-(mean(IF1)/mean(IF2)^2)*(IF2-mean(IF2)))/n)
  if (output.IF==FALSE) {
    return(c(point,point-1.96*SD,point+1.96*SD))
  } else {
    return(list(Estimation=c(point,point-1.96*SD,point+1.96*SD),
                IF1=IF1,
                IF2=IF2))
  }
}


#' Point estimates of the principal natural (in)direct effects based on the singly and multiply robust estimators
#'
#' @param data dataset
#' @param Xname names of the baseline covariates
#' @param Yname name of the outcome
#' @param Zname name of the treatment
#' @param Dname name of the post-treatment event
#' @param Mname name of the mediator
#' @param monotonicity "strong" or "standard" monotonicity assumption (default "standard")
#'
#' @returns Point estimate of the principal natural (in)direct effects
mediate_par_point = function(data,Xname,Yname,Zname,Dname,Mname,monotonicity="standard") {

  data = as.matrix(data)

  Xpi = data[,Xname,drop=FALSE]
  Xe  = data[,Xname,drop=FALSE]
  Xm  = data[,Xname,drop=FALSE]
  Xo  = data[,Xname,drop=FALSE]
  Z   = data[,Zname]
  D   = data[,Dname]
  M   = data[,Mname]
  Y   = data[,Yname]
  if (monotonicity=="strong") {noncompliance="one-sided"}
  if (monotonicity=="standard") {noncompliance="two-sided"}

  # calculate all GMFs
  M_c_10=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="10",noncompliance)
  M_c_11=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="11",noncompliance)
  M_c_00=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="00",noncompliance)
  if (noncompliance=="two-sided") {
    M_a_10=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="10",noncompliance)
    M_a_11=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="11",noncompliance)
    M_a_00=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="00",noncompliance)
  }
  M_n_10=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="10",noncompliance)
  M_n_11=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="11",noncompliance)
  M_n_00=GMF_par(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="00",noncompliance)
  # natural mediation effects
  NIE_c=M_c_11-M_c_10; NDE_c=M_c_10-M_c_00; TE_c=M_c_11-M_c_00
  NIE_n=M_n_11-M_n_10; NDE_n=M_n_10-M_n_00; TE_n=M_n_11-M_n_00
  if (noncompliance=="two-sided") {
    NIE_a=M_a_11-M_a_10; NDE_a=M_a_10-M_a_00; TE_a=M_a_11-M_a_00
  }
  if (noncompliance=="two-sided") {
    NIE = NIE_c*attr(NIE_c,"proportion")+NIE_n*attr(NIE_n,"proportion")+NIE_a*attr(NIE_a,"proportion")
    NDE = NDE_c*attr(NDE_c,"proportion")+NDE_n*attr(NDE_n,"proportion")+NDE_a*attr(NDE_a,"proportion")
    TE = TE_c*attr(TE_c,"proportion")+TE_n*attr(TE_n,"proportion")+TE_a*attr(TE_a,"proportion")
    out=matrix(c(NIE_c,NDE_c,TE_c,NIE_n,NDE_n,TE_n,
                 NIE_a,NDE_a,TE_a,NIE,NDE,TE),ncol=5,byrow=T)
    colnames(out)=c("a","b","c","d","mr")
    rownames(out)=c("NIE_c","NDE_c","TE_c","NIE_n","NDE_n","TE_n",
                    "NIE_a","NDE_a","TE_a","NIE","NDE","TE")
  } else {
    NIE = NIE_c*attr(NIE_c,"proportion")+NIE_n*attr(NIE_n,"proportion")
    NDE = NDE_c*attr(NDE_c,"proportion")+NDE_n*attr(NDE_n,"proportion")
    TE = TE_c*attr(TE_c,"proportion")+TE_n*attr(TE_n,"proportion")
    out=matrix(c(NIE_c,NDE_c,TE_c,NIE_n,NDE_n,TE_n,NIE,NDE,TE),ncol=5,byrow=T)
    colnames(out)=c("a","b","c","d","mr")
    rownames(out)=c("NIE_c","NDE_c","TE_c","NIE_n","NDE_n","TE_n","NIE","NDE","TE")
  }
  out
}


#' Singly and multiply robust estimation of the principal natural (in)direct effects
#'
#' @param data dataset
#' @param estimator For singly robust estimation, set it to "a", "b", "c", "d" for estimator a, b, c, d, respectively. For multiply robust estimation, set it to "mr"
#' @param Xname names of the baseline covariates
#' @param Yname name of the outcome
#' @param Zname name of the treatment
#' @param Dname name of the post-treatment event
#' @param Mname name of the mediator
#' @param monotonicity "strong" or "standard" monotonicity assumption (default "standard")
#' @param B number of iterations for the nonparametric bootstrap procedure
#'
#' @returns Point, SE, and 95% confidence interval estimates of the principal natural (in)direct effects
#' @examples
#' # set.seed(2024)
#' # data = data_gen(n=500)
#' # mediate_par(data,
#' #             estimator="mr",
#' #             Xname=c("X1","X2","X3","X4"),
#' #             Yname="Y",
#' #             Zname="Z",
#' #             Dname="D",
#' #             Mname="M",
#' #             monotonicity="standard",
#' #             B=200)
mediate_par = function(data,estimator="a",Xname,Yname,Zname,Dname,Mname,monotonicity="standard",B=500) {
  n=dim(data)[1]
  p=ifelse(monotonicity=="strong",5*9,5*12)
  out=matrix(NA,nrow=B,ncol=p)
  out[1,]=c(mediate_par_point(data,Xname,Yname,Zname,Dname,Mname,monotonicity))
  for (j in (2:B)) {
    ind = sample(1:n,n,replace=TRUE)
    out[j,]=c(mediate_par_point(data[ind,],Xname,Yname,Zname,Dname,Mname,monotonicity))
  }
  ci_low = matrix(apply(out,2,quantile, probs=0.025,na.rm=TRUE),ncol=5)
  ci_up = matrix(apply(out,2,quantile, probs=0.975,na.rm=TRUE),ncol=5)
  SE = matrix(apply(out,2,sd,na.rm=TRUE),ncol=5)
  res = list(point= matrix(out[1,],ncol=5),SE=SE, CI_low=ci_low,CI_up=ci_up)

  # estimators a, b, c, d
  res_a=cbind(res$point[,1],res$SE[,1],res$CI_low[,1],res$CI_up[,1])
  res_b=cbind(res$point[,2],res$SE[,2],res$CI_low[,2],res$CI_up[,2])
  res_c=cbind(res$point[,3],res$SE[,3],res$CI_low[,3],res$CI_up[,3])
  res_d=cbind(res$point[,4],res$SE[,4],res$CI_low[,4],res$CI_up[,4])
  # multiply robust estimator
  res_mr=cbind(res$point[,5],res$SE[,5],res$CI_low[,5],res$CI_up[,5])
  if (estimator=="a") out = res_a
  if (estimator=="b") out = res_b
  if (estimator=="c") out = res_c
  if (estimator=="d") out = res_d
  if (estimator=="mr") out = res_mr
  if (monotonicity=="strong") {
    rownames(out)=c("compliers: PNIE","compliers: PNDE","compliers: PCE",
                    "never-takers: PNIE","never-takers: PNDE","never-takers: PCE",
                    "ITT-NIE","ITT-NDE","ITT")
  } else {
    rownames(out)=c("compliers: PNIE","compliers: PNDE","compliers: PCE",
                    "never-takers: PNIE","never-takers: PNDE","never-takers: PCE",
                    "always-takers: PNIE","always-takers: PNDE","always-takers: PCE",
                    "ITT-NIE","ITT-NDE","ITT")
  }
  colnames(out) = c("Point","SE","CI_lower","CI_upper")
  out
}








#' Nonparametric efficient estimation of the principal natural (in)direct effects
#'
#' @param data dataset
#' @param Xname names of the baseline covariates
#' @param Yname name of the outcome
#' @param Zname name of the treatment
#' @param Dname name of the post-treatment event
#' @param Mname name of the mediator
#' @param monotonicity "strong" or "standard" monotonicity assumption (default "standard")
#' @param learners names of super learners (default GLM and random forest)
#' @param V number of folds for the cross-fitting procedure (default 5)
#'
#' @returns Point, SE, and 95% confidence interval estimates of the principal natural (in)direct effects
#' @examples
#' # set.seed(2024)
#' # data = data_gen(n=500)
#' # mediate_np(data,
#' #            Xname=c("X1","X2","X3","X4"),
#' #            Yname="Y",
#' #            Zname="Z",
#' #            Dname="D",
#' #            Mname="M",
#' #            monotonicity="standard",
#' #            learners=c("SL.glm", "SL.ranger"),V=3)
#'
mediate_np = function(data,Xname,Yname,Zname,Dname,Mname,monotonicity="standard",
                      learners=c("SL.glm", "SL.ranger"),V=5) {

  data = as.matrix(data)

  Xpi = data[,Xname,drop=FALSE]
  Xe  = data[,Xname,drop=FALSE]
  Xm  = data[,Xname,drop=FALSE]
  Xo  = data[,Xname,drop=FALSE]
  Z   = data[,Zname]
  D   = data[,Dname]
  M   = data[,Mname]
  Y   = data[,Yname]
  if (monotonicity=="strong") {noncompliance="one-sided"}
  if (monotonicity=="standard") {noncompliance="two-sided"}

  my_point_sd = function(m1,m2) {
    Point = m1$Estimation[1]-m2$Estimation[1]
    IF_m1 = 1/mean(m1$IF2)*(m1$IF1-mean(m1$IF1))-(mean(m1$IF1)/mean(m1$IF2)^2)*(m1$IF2-mean(m1$IF2))
    IF_m2 = 1/mean(m2$IF2)*(m2$IF1-mean(m2$IF1))-(mean(m2$IF1)/mean(m2$IF2)^2)*(m2$IF2-mean(m2$IF2))
    SD = sqrt(var(IF_m1-IF_m2)/length(IF_m1))
    return(c(Point,SD,Point-1.96*SD,Point+1.96*SD))
  }

  my_point_sd_itt_2sided = function(mc1,mc0,mn1,mn0,ma1,ma0) {
    Point = mean(mc1$IF1+mn1$IF1+ma1$IF1)-mean(mc0$IF1+mn0$IF1+ma0$IF1)
    SD = sqrt(var(mc1$IF1+mn1$IF1+ma1$IF1-mc0$IF1-mn0$IF1-ma0$IF1)/length(mc1$IF1))
    return(c(Point,SD,Point-1.96*SD,Point+1.96*SD))
  }

  my_point_sd_itt_1sided = function(mc1,mc0,mn1,mn0) {
    Point = mean(mc1$IF1+mn1$IF1)-mean(mc0$IF1+mn0$IF1)
    SD = sqrt(var(mc1$IF1+mn1$IF1-mc0$IF1-mn0$IF1)/length(mc1$IF1))
    return(c(Point,SD,Point-1.96*SD,Point+1.96*SD))
  }

  # calculate all GMFs
  M_c_10=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="10",noncompliance,learners,V=V,output.IF=TRUE)
  M_c_11=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="11",noncompliance,learners,V=V,output.IF=TRUE)
  M_c_00=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="compliers",type="00",noncompliance,learners,V=V,output.IF=TRUE)
  if (noncompliance=="two-sided") {
    M_a_10=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="10",noncompliance,learners,V=V,output.IF=TRUE)
    M_a_11=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="11",noncompliance,learners,V=V,output.IF=TRUE)
    M_a_00=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="always takers",type="00",noncompliance,learners,V=V,output.IF=TRUE)
  }
  M_n_10=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="10",noncompliance,learners,V=V,output.IF=TRUE)
  M_n_11=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="11",noncompliance,learners,V=V,output.IF=TRUE)
  M_n_00=GMF_nonpar(Xpi,Xe,Xm,Xo,Z,D,M,Y,group="never takers",type="00",noncompliance,learners,V=V,output.IF=TRUE)
  # natural mediation effects
  NIE_c=my_point_sd(M_c_11,M_c_10);
  NDE_c=my_point_sd(M_c_10,M_c_00);
  TE_c=my_point_sd(M_c_11,M_c_00);

  NIE_n=my_point_sd(M_n_11,M_n_10);
  NDE_n=my_point_sd(M_n_10,M_n_00);
  TE_n=my_point_sd(M_n_11,M_n_00);
  if (noncompliance=="two-sided") {
    NIE_a=my_point_sd(M_a_11,M_a_10);
    NDE_a=my_point_sd(M_a_10,M_a_00);
    TE_a=my_point_sd(M_a_11,M_a_00);
  }

  if (noncompliance=="two-sided") {
    NIE=my_point_sd_itt_2sided(M_c_11,M_c_10,M_n_11,M_n_10,M_a_11,M_a_10)
    NDE=my_point_sd_itt_2sided(M_c_10,M_c_00,M_n_10,M_n_00,M_a_10,M_a_00)
    TE =my_point_sd_itt_2sided(M_c_11,M_c_00,M_n_11,M_n_00,M_a_11,M_a_00)
    out=matrix(c(NIE_c,NDE_c,TE_c,NIE_n,NDE_n,TE_n,
                 NIE_a,NDE_a,TE_a,NIE , NDE, TE),ncol=4,byrow=T)
    colnames(out)=c("point","SE","CI_lower","CI_upper")
    rownames(out)=c("compliers: PNIE","compliers: PNDE","compliers: PCE",
                    "never-takers: PNIE","never-takers: PNDE","never-takers: PCE",
                    "always-takers: PNIE","always-takers: PNDE","always-takers: PCE",
                    "ITT-NIE","ITT-NDE","ITT")
  } else {
    NIE=my_point_sd_itt_1sided(M_c_11,M_c_10,M_n_11,M_n_10)
    NDE=my_point_sd_itt_1sided(M_c_10,M_c_00,M_n_10,M_n_00)
    TE =my_point_sd_itt_1sided(M_c_11,M_c_00,M_n_11,M_n_00)
    out=matrix(c(NIE_c,NDE_c,TE_c,NIE_n,NDE_n,TE_n,NIE,NDE,TE),ncol=4,byrow=T)
    colnames(out)=c("point","SE","CI_lower","CI_upper")
    rownames(out)=c("compliers: PNIE","compliers: PNDE","compliers: PCE",
                    "never-takers: PNIE","never-takers: PNDE","never-takers: PCE",
                    "ITT-NIE","ITT-NDE","ITT")
  }
  out
}









