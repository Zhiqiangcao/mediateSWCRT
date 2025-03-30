#' Mediation analysis function for continuous outcome and binary mediator in a stepped wedge design 
#' based on Hussey and Hughes model with a constant treatment effect
#' 
#' After obtaining parameter estimates from linear mixed effect model for outcome model and 
#' generalized linear mixed effect model for mediator model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param method Two approximate methods are provided, one is STA (second-order Taylor approximate), the other is GHQ (Gauss-Hermite Quadrature)
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design. When there is a implementation period, then corresponding treatment status 
#'       should be set to -1, which is mainly for convenience
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param exposure The exposure time with J levels (e.g., 0,1,2,...,J-1 when there is no implementation period). If dataset has no exposure, then set it to NULL
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression  (default is NULL)
#'@param a0 The reference treatment level in defining TE and NIE. The default value is 0
#'@param a1 The treatment level of interest in defining TE and NIE. The default value is 1
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP.
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 glmer  
#'@importFrom lme4 lmer
#'@importFrom stats as.formula  
#'@importFrom stats binomial
#'@importFrom stats qt 
#'@importFrom stats na.omit
#'@importFrom stats median
#'
#'@examples 
#' library(lme4)
#' I = 15 
#' J = 4
#' n = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' sigma_a = 0.334
#' sigma_tau = 0.605   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=1.2,eta=0.6,
#' sigma_a,sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=1)
#' # example 1: mediation analysis without covariates in outcome and mediator models 
#' res1 = mediate_contY_binaM_hhm(data=mydata1, method = "STA")
#' print(res1)
#' res1f = mediate_contY_binaM_hhm(data=mydata1, method = "GHQ")
#' print(res1f)
#' 
#' # example 2: mediation analysis with different covariates in outcome and mediator models
#' # generate two covariates
#' covdata = data.frame("X1" = numeric(), "X2" = numeric())
#' set.seed(100)
#' for(i in 1:I){
#'   for(j in 1:J){
#'      x1_ijk = rnorm(n,mean=0.5) 
#'      x2_ijk = rnorm(n,sd=0.5)
#'      covdata = rbind(covdata, data.frame(cbind(X1 = x1_ijk, X2 = x2_ijk)))
#'   }
#' }
#' 
#' mydata2 = data.frame(mydata1,covdata)
#' res2 = mediate_contY_binaM_hhm(data=mydata2, method = "STA",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' res2f = mediate_contY_binaM_hhm(data=mydata2, method = "GHQ",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2f)
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_contY_binaM_hhm(data=mydata2, method = "STA",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' res3f = mediate_contY_binaM_hhm(data=mydata2, method = "GHQ",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3f)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(123456)
#' mydata3 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=1,eta=0.8,
#' sigma_a,sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=1)
#' res4 = mediate_contY_binaM_hhm(data=mydata3,method = "STA")
#' print(res4)
#' res4f = mediate_contY_binaM_hhm(data=mydata3,method = "GHQ")
#' print(res4f)
#' 
#' # example 5: if there are no exposure variable
#' mydata4 = mydata3[,-4]
#' res5 = mediate_contY_binaM_hhm(data=mydata4,method = "STA",exposure=NULL)
#' print(res5)


mediate_contY_binaM_hhm = function(data, method = "STA", outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
                                    period = "period", exposure = "E", covariateY = NULL, covariateM = NULL, a0 = 0, a1 = 1){
  #first test whether there are "Y","M","A","cluster","period" in data
  data = as.data.frame(data)
  require_covariate = c(outcome,mediator,treatment,cluster,period,covariateY,covariateM)
  all_covariate = colnames(data)
  if(!all(require_covariate %in% all_covariate)){
    stop("Some specified required arguments are not in the data, please check!")
  }
  #check missing values in interested variables of data.frame
  if(!is.null(covariateY) | !is.null(covariateM)){
    #at least one covariate in covariateY or covariateM
    covariateYM = union(covariateM,covariateY) #extract common covariates
    data = data[,c(outcome,mediator,treatment,cluster,period,covariateYM)]
  }else{ #in this case, both covariateY = covariateM = NULL
    data = data[,c(outcome,mediator,treatment,cluster,period)]
  }
  #then remove all rows that contains at least one NA for those interested variables in dataset
  if(sum(is.na(data))>0){
    data = na.omit(data)
    warning("NA values detected in variables of data set, and we delete all rows that contains at NA before conducting mediation analysis!")
  }
  #if there is implementation period, i.e., treatment = -1, then we reset outcome and mediator to be NA
  if(sum(data[,3]==-1) > 0){ #for those implementation period
    index = data[,3]==-1
    data[index,1] = NA  #set outcome to be missing
    data[index,2] = NA  #set mediator to be missing
  }
  
  #outcome model and mediator model formulas
  if (is.null(covariateY)) {
    formula_Y = as.formula(paste(outcome, "~", period, "+", treatment, 
                                 "+", mediator, "+", "(1|cluster)", sep = ""))
  }
  else {
    formula_Y = as.formula(paste(outcome, "~", period, "+", treatment, "+", mediator,
                                 "+", paste(covariateY, collapse = "+"), 
                                 "+", "(1|cluster)",sep = ""))
  }
  if (is.null(covariateM)) {
    formula_M = as.formula(paste(mediator, "~", period, "+", treatment, 
                                 "+", "(1|cluster)", sep = ""))
  }
  else {
    formula_M = as.formula(paste(mediator, "~",period, "+", treatment, 
                                 "+", paste(covariateM, collapse = "+"),
                                 "+", "(1|cluster)", sep = ""))
  }
  #Gauss-Hermite Quadrature method
  HermiteCoefs = function(order) {
    x <- 1
    if (order > 0) 
      for (n in 1:order) x <- c(0, 2 * x) - c(((0:(n - 1)) * x)[-1L], 0, 0)
    return(x)
  }
  
  gauss.hermite = function(f, mu = 0, sd = 1, ..., order = 5) {
    stopifnot(is.function(f))
    stopifnot(length(mu) == 1)
    stopifnot(length(sd) == 1)
    Hn <- HermiteCoefs(order)
    Hn1 <- HermiteCoefs(order - 1)
    x <- sort(Re(polyroot(Hn)))
    Hn1x <- matrix(Hn1, nrow = 1) %*% t(outer(x, 0:(order - 1), "^"))
    w <- 2^(order - 1) * factorial(order) * sqrt(pi)/(order * Hn1x)^2
    ww <- w/sqrt(pi)
    xx <- mu + sd * sqrt(2) * x
    ans <- 0
    for (i in seq_along(x)) ans <- ans + ww[i] * f(xx[i], ...)
    return(ans)
  }
  
  #Gauss-Hermite Quadrature to compute kappa 
  kappa_GHQ = function(gammaj,eta,gammax,sigma_tau){
    f1 = function(x) {
      exp(gammaj + eta*a1 + sum(xm_j*gammax) +  x)/      #xm_j is a vector of median
        (1 + exp(gammaj + eta*a1 + sum(xm_j*gammax) + x))
    }
    f0 = function(x) {
      exp(gammaj + eta*a0 + sum(xm_j*gammax) +  x)/
        (1 + exp(gammaj + eta*a0 + sum(xm_j*gammax) + x))
    }
    p1s = gauss.hermite(f = f1, mu = 0, sd = sigma_tau, order = 40)
    p0s = gauss.hermite(f = f0, mu = 0, sd = sigma_tau, order = 40)
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  
  #Approximate likelihood estimation method to compute kappa
  Kappa_STA = function(gammaj,eta,gammax,sigma_tau){
    mu1 = exp(gammaj+eta*a1+sum(xm_j*gammax))/(1+exp(gammaj+eta*a1+sum(xm_j*gammax)))
    mu0 = exp(gammaj+eta*a0+sum(xm_j*gammax))/(1+exp(gammaj+eta*a0+sum(xm_j*gammax)))
    p1s = mu1+0.5*(mu1-3*mu1^2+2*mu1^3)*sigma_tau^2
    p0s = mu0+0.5*(mu0-3*mu0^2+2*mu0^3)*sigma_tau^2
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  
  #number of period
  J = length(unique(data[,5]))  #period
  clusterid = data[,4] #cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment)] = lapply(data[,c(period, cluster, treatment)], factor)
  #fit outcome model and mediator model
  model_Y = summary(lmer(formula_Y, data = data))
  model_M = summary(glmer(formula_M, data = data, family=binomial))
  theta_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  eta_hat = model_M$coefficients[J+1,1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  #standard error estimation for \tau_i in mediation model
  sigmatau = sqrt(sigmatau$cluster[1])
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    gamma_x = model_M$coefficients[-(1:(J+1)),1]
    xm = as.matrix(data[,names(gamma_x)],ncol=length(names(gamma_x)))
    nm = dim(xm)[2]
  }
  if(method == "GHQ"){
    #point estimate
    nie_j = numeric(J) 
    for(j in 1:J){
      gammaj = gamma_est[j]
      if(nm == 0) xm_j = 0 else {
        xm_j = as.matrix(xm[data$period==j,],ncol=nm)
        #modify, using median value instead
        xm_j = apply(xm_j,2,median)
      }
      niej = kappa_GHQ(gammaj=gammaj,eta=eta_hat,gammax=gamma_x,sigma_tau=sigmatau)
      nie_j[j] = beta_m_hat*(niej$trt-niej$ctrl)
    }
    nde_j = theta_hat
    te_j = nie_j+nde_j
    mp_m = mean(nie_j)/mean(te_j) 
    point_est = c(mean(nie_j),mean(nde_j),mean(te_j),mp_m)
    names(point_est) = c("NIE", "NDE", "TE", "MP")
    #Jackknife variance
    betaG = matrix(0,ng,4)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[J+2,1]
      theta_hat_g = model_Y_g$coefficients[J+1,1]
      eta_hat_g = model_M_g$coefficients[J+1,1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0 
      }else{
        gamma_x_g = model_M_g$coefficients[-(1:(J+1)),1]
        xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
        nm = dim(xm)[2]
      }
      nie_jg = numeric(J)
      for(j in 1:J){
        gammaj_g = gamma_est_g[j]
        if(nm == 0) xm_j = 0 else {
          xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
          #modify, using median value instead
          xm_j = apply(xm_j,2,median)
        }
        niej_g = kappa_GHQ(gammaj=gammaj_g,eta=eta_hat_g,gammax=gamma_x_g,sigma_tau=sigmatau_g)
        nie_jg[j] = beta_m_hat_g*(niej_g$trt-niej_g$ctrl)
      }
      nde_jg = theta_hat_g
      te_jg = nie_jg+nde_jg
      mp_mg = mean(nie_jg)/mean(te_jg) 
      point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
      betaG[g,] = point_est_g
    }
  }else{#STA method
    #point estimate
    nie_j = numeric(J) 
    for(j in 1:J){
      gammaj = gamma_est[j]
      if(nm == 0) xm_j = 0 else {
        xm_j = as.matrix(xm[data$period==j,],ncol=nm)
        #modify, using median value instead
        xm_j = apply(xm_j,2,median)
      }
      niej = Kappa_STA(gammaj=gammaj,eta=eta_hat,gammax=gamma_x,sigma_tau=sigmatau)
      nie_j[j] = beta_m_hat*(niej$trt-niej$ctrl)
    }
    nde_j = theta_hat
    te_j = nie_j+nde_j
    mp_m = mean(nie_j)/mean(te_j)  
    point_est = c(mean(nie_j),mean(nde_j),mean(te_j),mp_m)
    names(point_est) = c("NIE", "NDE", "TE", "MP")
    #Jackknife variance
    betaG = matrix(0,ng,4)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[J+2,1]
      theta_hat_g = model_Y_g$coefficients[J+1,1]
      eta_hat_g = model_M_g$coefficients[J+1,1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0
      }else{
        gamma_x_g = model_M_g$coefficients[-(1:(J+1)),1]
        xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
        nm = dim(xm)[2]
      }
      nie_jg = numeric(J)
      for(j in 1:J){
        gammaj_g = gamma_est_g[j]
        if(nm == 0) xm_j = 0 else {
          xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
          #modify, using median value instead
          xm_j = apply(xm_j,2,median)
        }
        niej_g = Kappa_STA(gammaj=gammaj_g,eta=eta_hat_g,gammax=gamma_x_g,sigma_tau=sigmatau_g)
        nie_jg[j] = beta_m_hat_g*(niej_g$trt-niej_g$ctrl)
      }
      nde_jg = theta_hat_g
      te_jg = nie_jg+nde_jg
      mp_mg = mean(nie_jg)/mean(te_jg)  
      point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
      betaG[g,] = point_est_g
    }
  }
  # warnings
  # check for NA values in mediation measures
  params = list(nie_j = nie_j,nde_j=nde_j)
  # identify parameters with NA values
  na_params = function(params){sapply(params,function(x) any(is.na(x)))}
  # if any NA values are found, print a warning
  if(sum(na_params(params)) > 0){
    warning("NA values detected in mediation measures!")
  }
  # check MP is out of range [0,1] or not
  if(mp_m < 0 | mp_m > 1){
    warning("MP is out of range [0,1]!")
  }
  #check NA in Jackknife estimation
  if(sum(is.na(betaG))>0){
    warning("NA values detected in Jackknife estimations for mediation measures!")
  }
  
  #Jackknife variance
  betag_tilde = colMeans(betaG)
  var_est = colSums(sweep(betaG,2,betag_tilde,FUN="-")^2)*((ng-1)/ng)
  names(var_est) = c("NIE", "NDE", "TE", "MP")
  sd_est = sqrt(var_est)
  names(sd_est) = c("NIE", "NDE", "TE", "MP")
  qta = qt(0.975,ng-1) 
  ci_est = rbind(point_est - qta * sd_est, point_est + qta * sd_est)
  rownames(ci_est) = c("ci_lower_confidence_limit", "ci_upper_confidence_limit")
  return(list(point_est = point_est, var_est = var_est, sd_est = sd_est, 
              ci_est = ci_est))
}
