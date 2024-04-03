#' Mediation analysis function for continuous outcome and binary mediator in a stepped wedge with a time-dependent treatment effect
#'
#' After obtaining parameter estimates from linear mixed effect model for outcome model and generalized linear mixed-effects for mediator model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP with a time-dependent treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param method two approximate methods are provided, one is GHQ (Gauss-Hermite Quadrature), the other is AMLE (approximate maximum likelihood estimation)
#'@param id The id of subjects in stepped wedge design
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design, here treatment is a factor variable with several levels
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression  (default is NULL)
#'@param a0 The baseline treatment level
#'@param a1 The new treatment level
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
#'
#'@examples
#' library(lme4)
#' n = 15 
#' J = 4
#' m = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' gammas = c(0.5,1.3,2.1) 
#' betas = c(0.6,1,1.4)
#' sigma_a = 0.334
#' sigma_tau = 0.605   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=1.2,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=1)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_contY_binaM_etm(data=mydata1, method = "GHQ")
#' print(res1)
#' res1f = mediate_contY_binaM_etm(data=mydata1, method = "AMLE")
#' print(res1f)
#' 
#' # example 2: mediation analysis with different covariates in both outcome and mediator models
#' # generate two covariates
#' covdata = data.frame("X1" = numeric(), "X2" = numeric())
#' set.seed(100)
#' for(i in 1:n){
#'   for(j in 1:J){
#'      x1_ijk = rnorm(m) 
#'      x2_ijk = rnorm(m,sd=0.5)
#'      covdata <- rbind(covdata, data.frame(cbind(X1 = x1_ijk, X2 = x2_ijk)))
#'   }
#' }
#' 
#' mydata2 = data.frame(mydata1,covdata)
#' res2 = mediate_contY_binaM_etm(data=mydata2, method = "GHQ",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' res2f = mediate_contY_binaM_etm(data=mydata2, method = "AMLE",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2f)
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_contY_binaM_etm(data=mydata2, method = "GHQ",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' res3f = mediate_contY_binaM_etm(data=mydata2, method = "AMLE",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3f)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' gammas = c(0.5,1.2,1.9,2.6) 
#' betas = c(0.3,0.6,0.9,1.2)
#' set.seed(1234)
#' mydata3 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=1.2,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=1)
#' res4 = mediate_contY_binaM_etm(data=mydata3,method = "GHQ")
#' print(res4)
#' res4f = mediate_contY_binaM_etm(data=mydata3,method = "AMLE")
#' print(res4f)
#'

mediate_contY_binaM_etm = function(data, method = "AMLE", id = "id", outcome = "Y", mediator = "M", treatment = "S", cluster = "cluster", 
                                   period = "period", covariateY = NULL, covariateM = NULL, a0 = 0, a1 = 1){
  
  data = as.data.frame(data)
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
  #gaussian hermite method
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
  kappa_GHQ = function(gammaj,gammaa,gammax,sigma_tau){
    f1 = function(x) {
      exp(gammaj + gammaa*a1 + as.numeric(xm_j%*%gammax) +  x)/
        (1 + exp(gammaj + gammaa*a1 + as.numeric(xm_j%*%gammax) + x))
    }
    f0 = function(x) {
      exp(gammaj + gammaa*a0 + as.numeric(xm_j%*%gammax) +  x)/
        (1 + exp(gammaj + gammaa*a0 + as.numeric(xm_j%*%gammax) + x))
    }
    p1s = gauss.hermite(f = f1, mu = 0, sd = sigma_tau, order = 40)
    p0s = gauss.hermite(f = f0, mu = 0, sd = sigma_tau, order = 40)
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  
  #Approximate likelihood estimation method to compute kappa
  Kappa_AMLE = function(gammaj,gammaa,gammax,sigma_tau){
    mu1 = exp(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax)))
    mu0 = exp(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax)))
    p1s = mu1+0.5*(mu1-3*mu1^2+2*mu1^3)*sigma_tau^2
    p0s = mu0+0.5*(mu0-3*mu0^2+2*mu0^3)*sigma_tau^2
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  #number of period
  J = max(data$period)
  clusterid = data$cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment, id)] = lapply(data[,c(period, cluster, treatment, id)], factor)
  #fit outcome model and mediator model
  model_Y = summary(lmer(formula_Y, data = data))
  model_M = summary(glmer(formula_M, data = data, family=binomial))
  beta_m_hat = model_Y$coefficients[2*J,1]
  beta_s_hat = model_Y$coefficients[(J+1):(2*J-1),1]
  gamma_s_hat = model_M$coefficients[(J+1):(2*J-1),1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    nm = length(covariateM)
    gamma_x = model_M$coefficients[(2*J):(2*J+nm-1),1]
    xm = as.matrix(data[,covariateM],ncol=nm)
  }
  if(method == "GHQ"){
    nie_s = nde_s = numeric(J-1)
    for(s in 1:(J-1)){
      gammas = gamma_s_hat[s]
      gammajs = gamma_est[(s+1):J] #corresponding s
      nie_js = numeric(J-s)
      for(j in 1:(J-s)){
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==(j+s),ncol=nm])
        niej = kappa_GHQ(gammaj=gammajs[j],gammaa=gammas,gammax = gamma_x, sigma_tau=sigmatau)
        nie_js[j] = mean(beta_m_hat*(niej$trt-niej$ctrl)) #average all i and k 
      } #average all j
      nie_s[s] = mean(nie_js)
      nde_s[s] = beta_s_hat[s]
    }
    te_s = nie_s + nde_s
    mp_s = nie_s/te_s
    mp_m = mean(nie_s)/mean(te_s)
    #point estimation
    point_est = data.frame(NIE = c(mean(nie_s),nie_s),NDE = c(mean(nde_s),nde_s),
                           TE = c(mean(te_s),te_s), MP = c(mp_m,mp_s))
    row.names(point_est) = c("overall",paste("extime",1:(J-1),sep=""))
    #Jacknife variance
    betaG = array(0,dim=c(ng,4,J))
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[2*J,1]
      beta_s_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
      gamma_s_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0; xm = 0 
      }else{
        gamma_x_g = model_M_g$coefficients[(2*J):(2*J+nm-1),1]
        xm = as.matrix(data_g[,covariateM],ncol=nm)
      }
      nie_s_g = nde_s_g = numeric(J-1)
      for(s in 1:(J-1)){
        gammas_g = gamma_s_hat_g[s]
        gammajs_g = gamma_est_g[(s+1):J] 
        nie_js_g = numeric(J-s)
        for(j in 1:(J-s)){
          if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==(j+s),ncol=nm])
          niej_g = kappa_GHQ(gammaj=gammajs_g[j],gammaa=gammas_g,gammax = gamma_x_g, sigma_tau=sigmatau_g)
          nie_js_g[j] = mean(beta_m_hat_g*(niej_g$trt-niej_g$ctrl))
        }
        nie_s_g[s] = mean(nie_js_g)
        nde_s_g[s] = beta_s_hat_g[s]
      }
      te_s_g = nie_s_g + nde_s_g
      mp_s_g = nie_s_g/te_s_g
      est_res_g = cbind(nie_s_g,nde_s_g,te_s_g,mp_s_g)
      for(s in 1:J){
        if(s == 1){ #overall NIE, NDE, TE and MP
          overall_m = colMeans(est_res_g)
          overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
          betaG[g,,s] = overall_m
        }else{
          betaG[g,,s] = est_res_g[s-1,]
        }
      }
    }
  }else{#AMLE method
    nie_s = nde_s = numeric(J-1)
    for(s in 1:(J-1)){
      gammas = gamma_s_hat[s]
      gammajs = gamma_est[(s+1):J] #corresponding s
      nie_js = numeric(J-s)
      for(j in 1:(J-s)){
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==(j+s),ncol=nm])
        niej = Kappa_AMLE(gammaj=gammajs[j],gammaa=gammas,gammax = gamma_x, sigma_tau=sigmatau)
        nie_js[j] = mean(beta_m_hat*(niej$trt-niej$ctrl)) #average all i and k
      } #average all j
      nie_s[s] = mean(nie_js)
      nde_s[s] = beta_s_hat[s]
    }
    te_s = nie_s + nde_s
    mp_s = nie_s/te_s
    mp_m = mean(nie_s)/mean(te_s)
    #point estimation
    point_est = data.frame(NIE = c(mean(nie_s),nie_s),NDE = c(mean(nde_s),nde_s),
                           TE = c(mean(te_s),te_s), MP = c(mp_m,mp_s))
    row.names(point_est) = c("overall", paste("extime",1:(J-1),sep=""))
    #Jacknife variance
    betaG = array(0,dim=c(ng,4,J))
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[2*J,1]
      beta_s_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
      gamma_s_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0 
      }else{
        gamma_x_g = model_M_g$coefficients[(2*J):(2*J+nm-1),1]
        xm = as.matrix(data_g[,covariateM],ncol=nm)
      }
      nie_s_g = nde_s_g = numeric(J-1)
      for(s in 1:(J-1)){
        gammas_g = gamma_s_hat_g[s]
        gammajs_g = gamma_est_g[(s+1):J] 
        nie_js_g = numeric(J-s)
        for(j in 1:(J-s)){
          if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==(j+s),ncol=nm])
          niej_g = Kappa_AMLE(gammaj=gammajs_g[j],gammaa=gammas_g,gammax=gamma_x_g, sigma_tau=sigmatau_g)
          nie_js_g[j] = mean(beta_m_hat_g*(niej_g$trt-niej_g$ctrl))
        }
        nie_s_g[s] = mean(nie_js_g)
        nde_s_g[s] = beta_s_hat_g[s]
      }
      te_s_g = nie_s_g + nde_s_g
      mp_s_g = nie_s_g/te_s_g
      est_res_g = cbind(nie_s_g,nde_s_g,te_s_g,mp_s_g)
      for(s in 1:J){
        if(s == 1){ #overall NIE, NDE, TE and MP
          overall_m = colMeans(est_res_g)
          overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
          betaG[g,,s] = overall_m
        }else{
          betaG[g,,s] = est_res_g[s-1,]
        }
      }
    }
  }
  #Contional jackknife varaince
  var_est = matrix(0,J,4)
  for(j in 1:J){
    betaG_j = betaG[,,j]
    betag_tilde_j = colMeans(betaG_j)
    var_est_j = colSums(sweep(betaG_j,2,betag_tilde_j,FUN="-")^2)*((ng-1)/ng)
    var_est[j,] = var_est_j
  }
  colnames(var_est) = c("NIE", "NDE", "TE", "MP")
  row.names(var_est) = row.names(point_est)
  var_est = as.data.frame(var_est)
  sd_est = sqrt(var_est)
  qta = qt(0.975,ng-1) 
  ci_est_low = point_est - qta * sd_est
  ci_est_high = point_est + qta * sd_est
  return(list(point_est = point_est, var_est = var_est, sd_est = sd_est, 
              ci_low_boundary = ci_est_low,ci_high_boundary = ci_est_high))
}
