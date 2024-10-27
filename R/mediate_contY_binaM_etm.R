#' Mediation analysis function for continuous outcome and binary mediator in a stepped wedge 
#' with an exposure-time dependent treatment effect
#'
#' After obtaining parameter estimates from linear mixed effect model for outcome model and 
#' generalized linear mixed effect model for mediator model, this function can obtain mediation 
#' measures including NIE, NDE, TE and MP with an exposure-time dependent treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param method Two approximate methods are provided, one is STA (second-order Taylor approximate), the other is GHQ (Gauss-Hermite Quadrature)
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design, here treatment is exposure time E with J levels (e.g., 0,1,2,...,J-1)
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression  (default is NULL)
#'@param a0 The baseline treatment level
#'@param a1 The new treatment level
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP at each exposure time 
#'         and overall summary mediation measures, as well as Chi-square test result of TE
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 glmer  
#'@importFrom lme4 lmer
#'@importFrom stats as.formula  
#'@importFrom stats binomial
#'@importFrom stats qt 
#'@importFrom stats pchisq 
#'@importFrom stats na.omit
#'
#'@examples
#' library(lme4)
#' I = 15 
#' J = 4
#' n = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' eta_e = c(0.5,1.3,2.1) 
#' theta_e = c(0.6,1,1.4)
#' sigma_a = 0.334
#' sigma_tau = 0.605   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=1.2,eta_e,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=1)
#' # example 1: mediation analysis without covariates in outcome and mediator models
#' res1 = mediate_contY_binaM_etm(data=mydata1, method = "STA")
#' print(res1)
#' res1f = mediate_contY_binaM_etm(data=mydata1, method = "GHQ")
#' print(res1f)
#' 
#' # example 2: mediation analysis with different covariates in outcome and mediator models
#' # generate two covariates
#' covdata = data.frame("X1" = numeric(), "X2" = numeric())
#' set.seed(100)
#' for(i in 1:I){
#'   for(j in 1:J){
#'      x1_ijk = rnorm(n) 
#'      x2_ijk = rnorm(n,sd=0.5)
#'      covdata = rbind(covdata, data.frame(cbind(X1 = x1_ijk, X2 = x2_ijk)))
#'   }
#' }
#' 
#' mydata2 = data.frame(mydata1,covdata)
#' res2 = mediate_contY_binaM_etm(data=mydata2, method = "STA",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' res2f = mediate_contY_binaM_etm(data=mydata2, method = "GHQ",
#' covariateY = c("X1"), covariateM = c("X2"))
#' print(res2f)
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_contY_binaM_etm(data=mydata2, method = "STA",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' res3f = mediate_contY_binaM_etm(data=mydata2, method = "GHQ",
#' covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3f)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' eta_e = c(0.5,1.2,1.9,2.6) 
#' theta_e = c(0.3,0.6,0.9,1.2)
#' set.seed(123456)
#' mydata3 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=1.2,eta_e,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=1)
#' res4 = mediate_contY_binaM_etm(data=mydata3,method = "STA")
#' print(res4)
#' res4f = mediate_contY_binaM_etm(data=mydata3,method = "GHQ")
#' print(res4f)
#'

mediate_contY_binaM_etm = function(data, method = "STA", outcome = "Y", mediator = "M", treatment = "E", cluster = "cluster", 
                                    period = "period", covariateY = NULL, covariateM = NULL, a0 = 0, a1 = 1){
  #first test whether there are "Y","M","E","cluster","period" in data
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
      exp(gammaj + eta*a1 + as.numeric(xm_j%*%gammax) +  x)/
        (1 + exp(gammaj + eta*a1 + as.numeric(xm_j%*%gammax) + x))
    }
    f0 = function(x) {
      exp(gammaj + eta*a0 + as.numeric(xm_j%*%gammax) +  x)/
        (1 + exp(gammaj + eta*a0 + as.numeric(xm_j%*%gammax) + x))
    }
    p1s = gauss.hermite(f = f1, mu = 0, sd = sigma_tau, order = 40)
    p0s = gauss.hermite(f = f0, mu = 0, sd = sigma_tau, order = 40)
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  
  #Approximate likelihood estimation method to compute kappa
  Kappa_STA = function(gammaj,eta,gammax,sigma_tau){
    mu1 = exp(gammaj+eta*a1+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+eta*a1+as.numeric(xm_j%*%gammax)))
    mu0 = exp(gammaj+eta*a0+as.numeric(xm_j%*%gammax))/(1+exp(gammaj+eta*a0+as.numeric(xm_j%*%gammax)))
    p1s = mu1+0.5*(mu1-3*mu1^2+2*mu1^3)*sigma_tau^2
    p0s = mu0+0.5*(mu0-3*mu0^2+2*mu0^3)*sigma_tau^2
    res = data.frame(trt=p1s,ctrl=p0s)
    return(res)
  }
  
  #number of period
  J = length(unique(data$period))
  clusterid = data$cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment)] = lapply(data[,c(period, cluster, treatment)], factor)
  #fit outcome model and mediator model
  model_Y = summary(lmer(formula_Y, data = data))
  model_M = summary(glmer(formula_M, data = data, family=binomial))
  beta_m_hat = model_Y$coefficients[2*J,1]
  theta_e_hat = model_Y$coefficients[(J+1):(2*J-1),1]
  eta_e_hat = model_M$coefficients[(J+1):(2*J-1),1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    gamma_x = model_M$coefficients[-(1:(2*J-1)),1]
    xm = as.matrix(data[,names(gamma_x)],ncol=length(names(gamma_x)))
    nm = dim(xm)[2]
  }
  if(method == "GHQ"){
    nie_e = nde_e = numeric(J-1)
    for(e in 1:(J-1)){
      etae = eta_e_hat[e]
      gammaje = gamma_est[(e+1):J] #corresponding e
      nie_je = numeric(J-e)
      for(j in 1:(J-e)){
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==(j+e),],ncol=nm)
        niej = kappa_GHQ(gammaj=gammaje[j],eta=etae,gammax = gamma_x, sigma_tau=sigmatau)
        nie_je[j] = mean(beta_m_hat*(niej$trt-niej$ctrl)) #average all i and k 
      } #average all j
      nie_e[e] = mean(nie_je)
      nde_e[e] = theta_e_hat[e]
    }
    te_e = nie_e + nde_e
    mp_e = nie_e/te_e
    mp_m = mean(nie_e)/mean(te_e)
    #point estimation
    point_est = data.frame(NIE = c(nie_e,mean(nie_e)),NDE = c(nde_e,mean(nde_e)),
                           TE = c(te_e,mean(te_e)), MP = c(mp_e,mp_m))
    row.names(point_est) = c(paste("e=",1:(J-1),sep=""),"overall")
    #Jackknife variance
    betaG = array(0,dim=c(ng,4,J))
    #Test=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
    te_jack = matrix(0,ng,J-2)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[2*J,1]
      theta_e_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
      eta_e_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0; xm = 0 
      }else{
        gamma_x_g = model_M_g$coefficients[-(1:(2*J-1)),1]
        xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
        nm = dim(xm)[2]
      }
      nie_e_g = nde_e_g = numeric(J-1)
      for(e in 1:(J-1)){
        etae_g = eta_e_hat_g[e]
        gammaje_g = gamma_est_g[(e+1):J] 
        nie_je_g = numeric(J-e)
        for(j in 1:(J-e)){
          if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==(j+e),],ncol=nm)
          niej_g = kappa_GHQ(gammaj=gammaje_g[j],eta=etae_g,gammax = gamma_x_g, sigma_tau=sigmatau_g)
          nie_je_g[j] = mean(beta_m_hat_g*(niej_g$trt-niej_g$ctrl))
        }
        nie_e_g[e] = mean(nie_je_g)
        nde_e_g[e] = theta_e_hat_g[e]
      }
      te_e_g = nie_e_g + nde_e_g
      mp_e_g = nie_e_g/te_e_g
      est_res_g = cbind(nie_e_g,nde_e_g,te_e_g,mp_e_g)
      #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
      te_jack[g,] = te_e_g[1]-te_e_g[-1]
      for(e in 1:J){
        if(e < J){
          betaG[g,,e] = est_res_g[e,]
        }else{
          #overall NIE, NDE, TE and MP
          overall_m = colMeans(est_res_g)
          overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
          betaG[g,,e] = overall_m
        }
      }
    }
  }else{#STA method
    nie_e = nde_e = numeric(J-1)
    for(e in 1:(J-1)){
      etae = eta_e_hat[e]
      gammaje = gamma_est[(e+1):J] #corresponding e
      nie_je = numeric(J-e)
      for(j in 1:(J-e)){
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==(j+e),],ncol=nm)
        niej = Kappa_STA(gammaj=gammaje[j],eta=etae,gammax = gamma_x, sigma_tau=sigmatau)
        nie_je[j] = mean(beta_m_hat*(niej$trt-niej$ctrl)) #average all i and k
      } #average all j
      nie_e[e] = mean(nie_je)
      nde_e[e] = theta_e_hat[e]
    }
    te_e = nie_e + nde_e
    mp_e = nie_e/te_e
    mp_m = mean(nie_e)/mean(te_e)
    #point estimation
    point_est = data.frame(NIE = c(nie_e,mean(nie_e)),NDE = c(nde_e,mean(nde_e)),
                           TE = c(te_e,mean(te_e)), MP = c(mp_e,mp_m))
    row.names(point_est) = c(paste("e=",1:(J-1),sep=""),"overall")
    #Jackknife variance
    betaG = array(0,dim=c(ng,4,J))
    #Test=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
    te_jack = matrix(0,ng,J-2)
    for(g in 1:ng){
      g_exc = (clusterid != clusters[g])
      data_g = data[g_exc,]
      model_Y_g = summary(lmer(formula_Y, data = data_g))
      model_M_g = summary(glmer(formula_M, data = data_g, family=binomial))
      beta_m_hat_g = model_Y_g$coefficients[2*J,1]
      theta_e_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
      eta_e_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
      gamma_est_g = model_M_g$coefficients[1:J,1]
      sigmatau_g = model_M_g$varcor
      sigmatau_g = sqrt(sigmatau_g$cluster[1])
      if(nm == 0){
        gamma_x_g = 0 
      }else{
        gamma_x_g = model_M_g$coefficients[-(1:(2*J-1)),1]
        xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
        nm = dim(xm)[2]
      }
      nie_e_g = nde_e_g = numeric(J-1)
      for(e in 1:(J-1)){
        etae_g = eta_e_hat_g[e]
        gammaje_g = gamma_est_g[(e+1):J] 
        nie_je_g = numeric(J-e)
        for(j in 1:(J-e)){
          if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==(j+e),],ncol=nm)
          niej_g = Kappa_STA(gammaj=gammaje_g[j],eta=etae_g,gammax=gamma_x_g, sigma_tau=sigmatau_g)
          nie_je_g[j] = mean(beta_m_hat_g*(niej_g$trt-niej_g$ctrl))
        }
        nie_e_g[e] = mean(nie_je_g)
        nde_e_g[e] = theta_e_hat_g[e]
      }
      te_e_g = nie_e_g + nde_e_g
      mp_e_g = nie_e_g/te_e_g
      est_res_g = cbind(nie_e_g,nde_e_g,te_e_g,mp_e_g)
      #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
      te_jack[g,] = te_e_g[1]-te_e_g[-1]
      for(e in 1:J){
        if(e < J){
          betaG[g,,e] = est_res_g[e,]
        }else{  #overall NIE, NDE, TE and MP
          overall_m = colMeans(est_res_g)
          overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
          betaG[g,,e] = overall_m
        }
      }
    }
  }
  # warnings
  # check for NA values in mediation measures
  params = list(nie_e=nie_e,nde_e=nde_e,te_e=te_e,mp_e=mp_e)
  # identify parameters with NA values
  na_params = function(params){sapply(params,function(x) any(is.na(x)))}
  # if any NA values are found, print a warning
  if(sum(na_params(params)) > 0){
    warning("NA values detected in mediation measures!")
  }
  # check MP is out of range [0,1] or not
  MP = c(mp_e,mp_m)
  if(any(MP < 0) | any(MP > 1)){
    warning("MP is out of range [0,1]!")
  }
  #check NA in Jackknife estimation
  if(sum(is.na(betaG))>0){
    warning("NA values detected in Jackknife estimations for mediation measures!")
  }
  
  #Jackknife variance
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
  
  #construct chi-square test statistics
  chi_square_sta_cal = function(est,jack_est,ng,J){
    #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
    #est is point estimation of theta(e) at e=1,2,...,J-1
    #jack_est is a ng*(J-2) matrix for T
    point.est = matrix(est[1]-est[-1],ncol=1)  #(J-1)*1 matrix
    jack_est_m = colMeans(jack_est)
    jack_est_matrix = jack_est-matrix(rep(jack_est_m,each=ng),ncol=J-2)
    jack_var_est = ((ng-1)/ng)*(t(jack_est_matrix)%*%jack_est_matrix)
    T.stat = as.numeric(t(point.est)%*%solve(jack_var_est)%*%point.est)
    P.value = 1-pchisq(T.stat, df=J-2)
    resu = c(T.stat,P.value)
    return(resu)
  }
  #test TE(e) can be treated as constant
  chisq_test = chi_square_sta_cal(point_est[-J,3],te_jack,ng,J)
  names(chisq_test) = c("Test.stat","P.value")
  
  return(list(point_est = point_est, var_est = var_est, sd_est = sd_est, 
              ci_low_boundary = ci_est_low,ci_high_boundary = ci_est_high,
              chisq_test = chisq_test))
}
