#' Mediation analysis function for binary outcome and binary mediator in a stepped wedge design based on Hussey and Hughes model
#'
#' After obtaining parameter estimates from generalized linear mixed-effects model for both mediator model and outcome model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param method two approximate methods are provided, one is GHQ (Gauss-Hermite Quadrature), the other is AMLE (approximate maximum likelihood estimation)
#'@param id The id of subjects in stepped wedge design
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design
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
#' sigma_a = 0.605
#' sigma_tau = 0.605   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_hhm(n,J,m,beta,gamma,betaA=1,betaM=1,gammaA=1,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=1)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_binaY_binaM_hhm(data=mydata1, method = "GHQ")
#' print(res1)
#' res1f = mediate_binaY_binaM_hhm(data=mydata1, method = "AMLE")
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
#' res2 = mediate_binaY_binaM_hhm(data=mydata2, method = "GHQ",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' res2f = mediate_binaY_binaM_hhm(data=mydata2, method = "AMLE",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2f)
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_binaY_binaM_hhm(data=mydata2, method = "GHQ",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' res3f = mediate_binaY_binaM_hhm(data=mydata2, method = "AMLE",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3f)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 16  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(1234)
#' mydata3 = gen_data_hhm(n,J,m,beta,gamma,betaA=1,betaM=1,gammaA=1,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=1)
#' res4 = mediate_binaY_binaM_hhm(data=mydata3,method = "GHQ")
#' print(res4)
#' res4f = mediate_binaY_binaM_hhm(data=mydata3,method = "AMLE")
#' print(res4f)

mediate_binaY_binaM_hhm = function(data, method = "AMLE", id = "id", outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
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
  
  #Gauss-Hermite Quadrature to compute eta 
  eta_GHQ = function(gammaj,gammaa,gammax,sigma_tau){
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
    res = data.frame(p1s,p0s)
    return(res)
  }
  ##Gauss-Hermite Quadrature to compute kappa
  kappa_GHQ = function(betaj,betaa,betam,betax,sigma_a){
    f11 = function(x){
      exp(betaj + betaa*a1 + betam*a1 + as.numeric(xy_j%*%betax) + x)/(1+exp(betaj + betaa*a1 + betam*a1 + as.numeric(xy_j%*%betax) + x))
    }
    f10 = function(x){
      exp(betaj + betaa*a1 + betam*a0 + as.numeric(xy_j%*%betax) + x)/(1+exp(betaj + betaa*a1 + betam*a0 + as.numeric(xy_j%*%betax) + x))
    }
    f01 = function(x){
      exp(betaj + betaa*a0 + betam*a1 + as.numeric(xy_j%*%betax) + x)/(1+exp(betaj + betaa*a0 + betam*a1 + as.numeric(xy_j%*%betax) + x))
    }
    f00 = function(x){
      exp(betaj + betaa*a0 + betam*a0 + as.numeric(xy_j%*%betax) + x)/(1+exp(betaj + betaa*a0 + betam*a0 + as.numeric(xy_j%*%betax) + x))
    }
    p11s = gauss.hermite(f = f11, mu = 0, sd = sigma_a, order = 40)
    p10s = gauss.hermite(f = f10, mu = 0, sd = sigma_a, order = 40)
    p01s = gauss.hermite(f = f01, mu = 0, sd = sigma_a, order = 40)
    p00s = gauss.hermite(f = f00, mu = 0, sd = sigma_a, order = 40)
    res = data.frame(p11s,p10s,p01s,p00s)
    return(res)
  }
  #Approximate likelihood estimation method to compute probability
  cal_prob_AMLE = function(gammaj,gammaa,gammax,sigma_tau,betaj,betaa,betam,betax,sigma_a){
    mu1 = exp(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax)))
    mu0 = exp(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax)))
    eta1 = mu1+0.5*(mu1-3*mu1^2+2*mu1^3)*sigma_tau^2
    eta0 = mu0+0.5*(mu0-3*mu0^2+2*mu0^3)*sigma_tau^2
    mu11 = exp(betaj+betaa*a1+betam*a1+as.numeric(xy_j%*%betax))/(1+exp(betaj+betaa*a1+betam*a1+as.numeric(xy_j%*%betax)))
    mu10 = exp(betaj+betaa*a1+betam*a0+as.numeric(xy_j%*%betax))/(1+exp(betaj+betaa*a1+betam*a0+as.numeric(xy_j%*%betax)))
    mu01 = exp(betaj+betaa*a0+betam*a1+as.numeric(xy_j%*%betax))/(1+exp(betaj+betaa*a0+betam*a1+as.numeric(xy_j%*%betax)))
    mu00 = exp(betaj+betaa*a0+betam*a0+as.numeric(xy_j%*%betax))/(1+exp(betaj+betaa*a0+betam*a0+as.numeric(xy_j%*%betax)))
    kappa11 = mu11+0.5*(mu11-3*mu11^2+2*mu11^3)*sigma_a^2
    kappa10 = mu10+0.5*(mu10-3*mu10^2+2*mu10^3)*sigma_a^2
    kappa01 = mu01+0.5*(mu01-3*mu01^2+2*mu01^3)*sigma_a^2
    kappa00 = mu00+0.5*(mu00-3*mu00^2+2*mu00^3)*sigma_a^2
    res = data.frame(eta1,eta0,kappa11,kappa10,kappa01,kappa00)
    return(res)
  }
  
  #number of period
  J = max(data$period)
  clusterid = data$cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment, id)] = lapply(data[,c(period, cluster, treatment, id)], factor)
  #utcome model 
  model_Y = summary(glmer(formula_Y, data = data, family=binomial))
  beta_a_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  beta_est = model_Y$coefficients[1:J,1]
  sigmaa = model_Y$varcor
  sigmaa = sqrt(sigmaa$cluster[1])
  #mediator model
  model_M = summary(glmer(formula_M, data = data, family=binomial))
  gamma_a_hat = model_M$coefficients[J+1,1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  
  #check whether covariates are included in outcome and mediator models
  if(is.null(covariateY)){
    beta_x = 0; ny = 0
  }else{
    ny = length(covariateY)
    beta_x = model_Y$coefficients[(J+3):(J+2+ny),1]
    xy = as.matrix(data[,covariateY],ncol=ny)
  }
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    nm = length(covariateM)
    gamma_x = model_M$coefficients[(J+2):(J+1+nm),1]
    xm = as.matrix(data[,covariateM],ncol=nm)
  }
  
  if(method == "GHQ"){
    #point estimate
    nie_j = nde_j = numeric(J) 
    for(j in 1:J){
      gammaj = gamma_est[j]
      betaj = beta_est[j]
      if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==j,],ncol=nm)
      if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data$period==j,],ncol=ny)
      #Gaussian-Hermite Quadrature approach
      eta_res = eta_GHQ(gammaj=gammaj,gammaa=gamma_a_hat,gammax=gamma_x,sigma_tau=sigmatau)
      eta1 = eta_res$p1s; eta0 = eta_res$p0s
      kappa_res = kappa_GHQ(betaj=betaj,betaa=beta_a_hat,betam=beta_m_hat,betax=beta_x,sigma_a=sigmaa)
      kappa11 = kappa_res$p11s; kappa10 = kappa_res$p10s
      kappa01 = kappa_res$p01s; kappa00 = kappa_res$p00s
      P111 = kappa10*(1-eta1) + kappa11*eta1
      P110 = kappa10*(1-eta0) + kappa11*eta0
      P010 = kappa00*(1-eta0) + kappa01*eta0
      niej = log(P111/(1-P111))-log(P110/(1-P110))
      ndej = log(P110/(1-P110))-log(P010/(1-P010))
      nie_j[j] = mean(niej)
      nde_j[j] = mean(ndej)
    }
    te_j = nie_j+nde_j
    mp_m = mean(nie_j)/mean(te_j)
    point_est = c(mean(nie_j),mean(nde_j),mean(te_j),mp_m)
    names(point_est) = c("NIE", "NDE", "TE", "MP")
    #Jacknife variance
    betaG = matrix(0,ng,4)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      #outcome model
      model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[J+2,1]
      beta_a_hat_g = model_Y_g$coefficients[J+1,1]
      beta_est_g = model_Y_g$coefficients[1:J,1]
      sigmaa_g = model_Y_g$varcor
      sigmaa_g = sqrt(sigmaa_g$cluster[1])
      #mediate model
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      gamma_a_hat_g = model_M_g$coefficients[J+1,1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      #check whether covariates are included in outcome and mediator models
      if(ny == 0){
        beta_x_g = 0
      }else{
        beta_x_g = model_Y_g$coefficients[(J+3):(J+2+ny),1]
        xy = as.matrix(data_g[,covariateY],ncol=ny)
      }
      if(nm == 0){
        gamma_x_g = 0
      }else{
        gamma_x_g = model_M_g$coefficients[(J+2):(J+1+nm),1]
        xm = as.matrix(data_g[,covariateM],ncol=nm)
      }
      nie_jg = nde_jg = numeric(J)
      for(j in 1:J){
        gammaj_g = gamma_est_g[j]
        betaj_g = beta_est_g[j]
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
        if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data_g$period==j,],ncol=ny)
        #Gaussian-Hermite Quadrature approach
        eta_res_g = eta_GHQ(gammaj=gammaj_g,gammaa=gamma_a_hat_g,gammax=gamma_x_g,sigma_tau=sigmatau_g)
        eta1_g = eta_res_g$p1s; eta0_g = eta_res_g$p0s
        kappa_res_g = kappa_GHQ(betaj=betaj_g,betaa=beta_a_hat_g,betam=beta_m_hat_g,betax=beta_x_g,sigma_a=sigmaa_g)
        kappa11_g = kappa_res_g$p11s; kappa10_g = kappa_res_g$p10s
        kappa01_g = kappa_res_g$p01s; kappa00_g = kappa_res_g$p00s
        P111_g = kappa10_g*(1-eta1_g) + kappa11_g*eta1_g
        P110_g = kappa10_g*(1-eta0_g) + kappa11_g*eta0_g
        P010_g = kappa00_g*(1-eta0_g) + kappa01_g*eta0_g
        niej_g = log(P111_g/(1-P111_g))-log(P110_g/(1-P110_g))
        ndej_g = log(P110_g/(1-P110_g))-log(P010_g/(1-P010_g))
        nie_jg[j] = mean(niej_g)
        nde_jg[j] = mean(ndej_g)
      }
      te_jg = nie_jg+nde_jg
      mp_mg = mean(nie_jg)/mean(te_jg)
      point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
      betaG[g,] = point_est_g
    }
  }else{#ALE method
    #point estimate
    nie_j = nde_j = numeric(J) 
    for(j in 1:J){
      gammaj = gamma_est[j]
      betaj = beta_est[j]
      if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==j,],ncol=nm)
      if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data$period==j,],ncol=ny)
      #approximate likelihood estimation method
      approx_res = cal_prob_AMLE(gammaj=gammaj,gammaa=gamma_a_hat,gammax=gamma_x,sigma_tau=sigmatau,
                                 betaj=betaj,betaa=beta_a_hat,betam=beta_m_hat,betax=beta_x,sigma_a=sigmaa)
      eta1 = approx_res$eta1;  eta0 = approx_res$eta0;
      kappa11 = approx_res$kappa11; kappa10 = approx_res$kappa10
      kappa01 = approx_res$kappa01; kappa00 = approx_res$kappa00
      P111 = kappa10*(1-eta1) + kappa11*eta1
      P110 = kappa10*(1-eta0) + kappa11*eta0
      P010 = kappa00*(1-eta0) + kappa01*eta0
      niej = log(P111/(1-P111))-log(P110/(1-P110))
      ndej = log(P110/(1-P110))-log(P010/(1-P010))
      nie_j[j] = mean(niej)
      nde_j[j] = mean(ndej)
    }
    te_j = nie_j+nde_j
    mp_m = mean(nie_j)/mean(te_j)
    point_est = c(mean(nie_j),mean(nde_j),mean(te_j),mp_m)
    names(point_est) = c("NIE", "NDE", "TE", "MP")
    #Jacknife variance
    betaG = matrix(0,ng,4)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      #outcome model
      model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[J+2,1]
      beta_a_hat_g = model_Y_g$coefficients[J+1,1]
      beta_est_g = model_Y_g$coefficients[1:J,1]
      sigmaa_g = model_Y_g$varcor
      sigmaa_g = sqrt(sigmaa_g$cluster[1])
      #mediate model
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      gamma_a_hat_g = model_M_g$coefficients[J+1,1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      #check whether covariates are included in outcome and mediator models
      if(ny == 0){
        beta_x_g = 0
      }else{
        beta_x_g = model_Y_g$coefficients[(J+3):(J+2+ny),1]
        xy = as.matrix(data_g[,covariateY],ncol=ny)
      }
      if(nm == 0){
        gamma_x_g = 0
      }else{
        gamma_x_g = model_M_g$coefficients[(J+2):(J+1+nm),1]
        xm = as.matrix(data_g[,covariateM],ncol=nm)
      }
      nie_jg = nde_jg = numeric(J)
      for(j in 1:J){
        gammaj_g = gamma_est_g[j]
        betaj_g = beta_est_g[j]
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
        if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data_g$period==j,],ncol=ny)
        #approximate likelihood estimation
        approx_res_g = cal_prob_AMLE(gammaj=gammaj_g,gammaa=gamma_a_hat_g,gammax=gamma_x_g,sigma_tau=sigmatau_g,
                                     betaj=betaj_g,betaa=beta_a_hat_g,betam=beta_m_hat_g,betax=beta_x_g,sigma_a=sigmaa_g)
        eta1_g = approx_res_g$eta1;  eta0_g = approx_res_g$eta0;
        kappa11_g = approx_res_g$kappa11; kappa10_g = approx_res_g$kappa10
        kappa01_g = approx_res_g$kappa01; kappa00_g = approx_res_g$kappa00
        P111_g = kappa10_g*(1-eta1_g) + kappa11_g*eta1_g
        P110_g = kappa10_g*(1-eta0_g) + kappa11_g*eta0_g
        P010_g = kappa00_g*(1-eta0_g) + kappa01_g*eta0_g
        niej_g = log(P111_g/(1-P111_g))-log(P110_g/(1-P110_g))
        ndej_g = log(P110_g/(1-P110_g))-log(P010_g/(1-P010_g))
        nie_jg[j] = mean(niej_g)
        nde_jg[j] = mean(ndej_g)
      }
      te_jg = nie_jg+nde_jg
      mp_mg = mean(nie_jg)/mean(te_jg)
      point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
      betaG[g,] = point_est_g
    }
  }
  #Contional jackknife varaince
  betag_tilde = colMeans(betaG)
  var_est = colSums(sweep(betaG,2,betag_tilde,FUN="-")^2)*((ng-1)/ng)
  names(var_est) = c("NIE", "NDE", "TE", "MP")
  sd_est = sqrt(var_est)
  names(sd_est) = c("NIE", "NDE", "TE", "MP")
  qta = qt(0.975,ng-1) 
  ci_est = rbind(point_est - qta * sd_est, point_est + qta * sd_est)
  rownames(ci_est) = c("Lower boundary", "Upper boundary")
  return(list(point_est = point_est, var_est = var_est, sd_est = sd_est, 
              ci_est = ci_est))
}
