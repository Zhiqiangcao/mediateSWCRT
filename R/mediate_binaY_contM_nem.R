#' Mediation analysis function for binary outcome and continuous mediator in a stepped wedge design 
#' based on nested exchangeable correlation model with a constant treatment effect
#' 
#' After obtaining parameter estimates from linear mixed effect model for mediator model and
#' generalized linear mixed effect model for outcome model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design. When there is a implementation period, then corresponding treatment status 
#'       should be set to -1, which is mainly for convenience  
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param exposure The exposure time with J levels (e.g., 0,1,2,...,J-1 when there is no implementation period). If dataset has no exposure variable, then set it to NULL
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression (default is NULL)
#'@param a0 The reference treatment level in defining TE and NIE. The default value is 0
#'@param a1 The treatment level of interest in defining TE and NIE. The default value is 1
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP.
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 lmer   
#'@importFrom lme4 glmer
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
#' sigma_tau = sigma_a = 0.334
#' sigma_phi = sigma_psi = 0.6
#' sigma_em = sigma_ey = 0.8
#' set.seed(123456)
#' mydata1 = gen_data_nem(I,J,n,beta,gamma,theta=0.75,beta_M=1,eta=0.4,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=1,binary.m=0)
#' # example 1: mediation analysis without covariates in outcome and mediator models
#' res1 = mediate_binaY_contM_nem(data=mydata1)
#' print(res1)
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
#' res2 = mediate_binaY_contM_nem(data = mydata2,covariateY = c("X1"), 
#' covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_binaY_contM_nem(data = mydata2,covariateY = c("X1","X2"), 
#' covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12 
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(123456)
#' mydata3 = gen_data_nem(I,J,n,beta,gamma,theta=0.75,beta_M=1,eta=0.4,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=1,binary.m=0)
#' res4 = mediate_binaY_contM_nem(data=mydata3)
#' print(res4)
#' 

mediate_binaY_contM_nem = function(data, outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
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
    formula_Y = as.formula(paste(outcome, "~", period, "+", treatment, "+", mediator, 
                                 "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  else {
    formula_Y = as.formula(paste(outcome, "~", period, "+", treatment, "+", mediator,
                                 "+", paste(covariateY, collapse = "+"), 
                                 "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  if (is.null(covariateM)) {
    formula_M = as.formula(paste(mediator, "~", period, "+", treatment, 
                                 "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  else {
    formula_M = as.formula(paste(mediator, "~",period, "+", treatment, 
                                 "+", paste(covariateM, collapse = "+"),
                                 "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  #compute triple integral \mu (STA method)
  cal_prob = function(betaj,theta,betam,betax,gammaj,eta,gammax,sigma_em,sigma_tau,sigma_a,sigma_phi,sigma_psi){
    q_11 = exp(betaj+theta*a1+betam*(gammaj+eta*a1+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a1+betam*(gammaj+eta*a1+sum(xm_j*gammax))+sum(xy_j*betax)))
    A1_11 = q_11+(q_11-3*q_11^2+2*q_11^3)*0.5*sigma_phi^2
    A2_11 = q_11^2+(4*q_11^2-10*q_11^3+6*q_11^4)*0.5*sigma_phi^2
    A3_11 = q_11^3+(9*q_11^3-21*q_11^4+12*q_11^5)*0.5*sigma_phi^2
    A4_11 = q_11^4+(16*q_11^4-36*q_11^5+20*q_11^6)*0.5*sigma_phi^2
    A5_11 = q_11^5+(25*q_11^5-55*q_11^6+30*q_11^7)*0.5*sigma_phi^2
    
    p1_11 = A1_11
    p2_11 = (A1_11-3*A2_11+2*A3_11)*0.5*(sigma_a^2+betam^2*(sigma_em^2+sigma_tau^2+sigma_psi^2))
    p3_11 = (A1_11-15*A2_11+50*A3_11-60*A4_11+24*A5_11)*1/4*betam^2*sigma_a^2*(sigma_em^2+sigma_tau^2+sigma_psi^2)
    pr_11 = p1_11+p2_11+p3_11
    #
    q_10 = exp(betaj+theta*a1+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a1+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax)))
    A1_10 = q_10+(q_10-3*q_10^2+2*q_10^3)*0.5*sigma_phi^2
    A2_10 = q_10^2+(4*q_10^2-10*q_10^3+6*q_10^4)*0.5*sigma_phi^2
    A3_10 = q_10^3+(9*q_10^3-21*q_10^4+12*q_10^5)*0.5*sigma_phi^2
    A4_10 = q_10^4+(16*q_10^4-36*q_10^5+20*q_10^6)*0.5*sigma_phi^2
    A5_10 = q_10^5+(25*q_10^5-55*q_10^6+30*q_10^7)*0.5*sigma_phi^2
    
    p1_10 = A1_10
    p2_10 = (A1_10-3*A2_10+2*A3_10)*0.5*(sigma_a^2+betam^2*(sigma_em^2+sigma_tau^2+sigma_psi^2))
    p3_10 = (A1_10-15*A2_10+50*A3_10-60*A4_10+24*A5_10)*1/4*betam^2*sigma_a^2*(sigma_em^2+sigma_tau^2+sigma_psi^2)
    pr_10 = p1_10+p2_10+p3_10
    #
    q_00 = exp(betaj+theta*a0+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a0+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax)))
    A1_00 = q_00+(q_00-3*q_00^2+2*q_00^3)*0.5*sigma_phi^2
    A2_00 = q_00^2+(4*q_00^2-10*q_00^3+6*q_00^4)*0.5*sigma_phi^2
    A3_00 = q_00^3+(9*q_00^3-21*q_00^4+12*q_00^5)*0.5*sigma_phi^2
    A4_00 = q_00^4+(16*q_00^4-36*q_00^5+20*q_00^6)*0.5*sigma_phi^2
    A5_00 = q_00^5+(25*q_00^5-55*q_00^6+30*q_00^7)*0.5*sigma_phi^2
    
    p1_00 = A1_00
    p2_00 = (A1_00-3*A2_00+2*A3_00)*0.5*(sigma_a^2+betam^2*(sigma_em^2+sigma_tau^2+sigma_psi^2))
    p3_00 = (A1_00-15*A2_00+50*A3_00-60*A4_00+24*A5_00)*1/4*betam^2*sigma_a^2*(sigma_em^2+sigma_tau^2+sigma_psi^2)
    pr_00 = p1_00+p2_00+p3_00
    
    res = data.frame(pr_11,pr_10,pr_00)
    return(res)
  }
  
  #number of period
  J = length(unique(data[,5]))  #period
  clusterid = data[,4] #cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment)] = lapply(data[,c(period, cluster, treatment)], factor)
  #fit outcome model
  model_Y = summary(glmer(formula_Y, data = data, family=binomial))
  theta_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  beta_est = model_Y$coefficients[1:J,1]  
  sigma_yy = model_Y$varcor
  #standard error estimation for \alpha_i and \phi_{ij} in mediation model
  sigmaa = sqrt(sigma_yy$cluster[1])
  sigmaphi = sqrt(sigma_yy$`period:cluster`[1])
  #fit mediator model
  model_M = summary(lmer(formula_M, data = data))
  eta_hat = model_M$coefficients[J+1,1]
  gamma_est = model_M$coefficients[1:J,1] #\gamma_{0j} for j
  sigma_mm = model_M$varcor
  #standard error estimation for \tau_i and \psi_{ij} in mediation model
  sigmatau = sqrt(sigma_mm$cluster[1])
  sigmapsi = sqrt(sigma_mm$`period:cluster`[1])
  #standard error for error term in mediator model
  sigmaem = model_M$sigma
  if(is.null(covariateY)){
    beta_x = 0; ny = 0
  }else{
    #in case fixed-effect model matrix is rank deficient
    beta_x = model_Y$coefficients[-(1:(J+2)),1]
    xy = as.matrix(data[,names(beta_x)],ncol=length(names(beta_x)))
    ny = dim(xy)[2]
  }
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    gamma_x = model_M$coefficients[-(1:(J+1)),1]
    xm = as.matrix(data[,names(gamma_x)],ncol=length(names(gamma_x)))
    nm = dim(xm)[2]
  }
  #point estimate
  nie_j = nde_j = numeric(J) 
  for(j in 1:J){
    gammaj = gamma_est[j];
    betaj = beta_est[j];
    if(nm == 0) xm_j = 0 else {
      xm_j = as.matrix(xm[data$period==j,],ncol=nm)
      #modify, using median
      xm_j = apply(xm_j,2,median)
    }
    if(ny == 0) xy_j = 0 else {
      xy_j = as.matrix(xy[data$period==j,],ncol=ny)
      #modify, using median
      xy_j = apply(xy_j,2,median)
    }
    #STA method
    prob_res = cal_prob(betaj=betaj,theta=theta_hat,betam=beta_m_hat,betax=beta_x,gammaj=gammaj,
                        eta=eta_hat, gammax = gamma_x,sigma_em=sigmaem,sigma_tau=sigmatau,
                        sigma_a=sigmaa, sigma_phi = sigmaphi, sigma_psi = sigmaphi)
    pr11 = prob_res$pr_11
    pr10 = prob_res$pr_10
    pr00 = prob_res$pr_00
    nie_j[j] = log(pr11/(1-pr11))-log(pr10/(1-pr10))
    nde_j[j] = log(pr10/(1-pr10))-log(pr00/(1-pr00))
  }
  te_j = nie_j+nde_j
  mp_m = mean(nie_j)/mean(te_j)
  point_est = c(mean(nie_j),mean(nde_j),mean(te_j),mp_m)
  names(point_est) = c("NIE", "NDE", "TE", "MP")
  #Jackknife variance
  betaG = matrix(0,ng,4)
  for(g in 1:ng){
    g_exc = (clusterid != clusters[g])
    data_g = data[g_exc,]
    #fit outcome model
    model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
    theta_hat_g = model_Y_g$coefficients[J+1,1]
    beta_m_hat_g = model_Y_g$coefficients[J+2,1]
    beta_est_g = model_Y_g$coefficients[1:J,1]  
    sigma_yy_g = model_Y_g$varcor
    #standard error estimation for \alpha_i and \phi_{ij} in mediation model
    sigmaa_g = sqrt(sigma_yy_g$cluster[1])
    sigmaphi_g = sqrt(sigma_yy_g$`period:cluster`[1])
    #fit mediator model
    model_M_g = summary(lmer(formula_M, data = data_g))
    eta_hat_g = model_M_g$coefficients[J+1,1]
    gamma_est_g = model_M_g$coefficients[1:J,1] #\gamma_{0j} for j
    sigma_mm_g = model_M_g$varcor
    #standard error estimation for \tau_i and \psi_{ij} in mediation model
    sigmatau_g = sqrt(sigma_mm_g$cluster[1])
    sigmapsi_g = sqrt(sigma_mm_g$`period:cluster`[1])
    #standard error for error term in mediator model
    sigmaem_g = model_M_g$sigma
    #check whether covariates are included in outcome and mediator models
    if(ny == 0){
      beta_x_g = 0;
    }else{
      beta_x_g = model_Y_g$coefficients[-(1:(J+2)),1]
      xy = as.matrix(data_g[,names(beta_x_g)],ncol=length(names(beta_x_g)))
      ny = dim(xy)[2]
      xy = as.matrix(data_g[,covariateY],ncol=ny)
    }
    if(nm == 0){
      gamma_x_g = 0
    }else{
      gamma_x_g = model_M_g$coefficients[-(1:(J+1)),1]
      xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
      nm = dim(xm)[2]
    }
    nie_jg = nde_jg = numeric(J)
    for(j in 1:J){
      gammaj_g = gamma_est_g[j];
      betaj_g = beta_est_g[j];
      if(nm == 0) xm_j = 0 else {
        xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
        #modify, using median
        xm_j = apply(xm_j,2,median)
      }
      if(ny == 0) xy_j = 0 else {
        xy_j = as.matrix(xy[data_g$period==j,],ncol=ny)
        #modify, using median
        xy_j = apply(xy_j,2,median)
      }
      ##STA method
      prob_res_g = cal_prob(betaj=betaj_g,theta=theta_hat_g,betam=beta_m_hat_g,betax=beta_x_g,gammaj=gammaj_g,
                            eta=eta_hat_g,gammax = gamma_x_g,sigma_em=sigmaem_g,sigma_tau=sigmatau_g,
                            sigma_a=sigmaa_g,sigma_phi = sigmaphi_g,sigma_psi = sigmapsi_g)
      pr11_g = prob_res_g$pr_11
      pr10_g = prob_res_g$pr_10
      pr00_g = prob_res_g$pr_00
      nie_jg[j] = log(pr11_g/(1-pr11_g))-log(pr10_g/(1-pr10_g))
      nde_jg[j] = log(pr10_g/(1-pr10_g))-log(pr00_g/(1-pr00_g))
    }
    te_jg = nie_jg+nde_jg
    mp_mg = mean(nie_jg)/mean(te_jg)
    point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
    betaG[g,] = point_est_g
  }
  # warnings
  # check for NA values in mediation measures
  params = list(nie_j=nie_j,nde_j=nde_j,te_j=te_j)
  # identify parameters with NA values
  na_params = function(params){sapply(params,function(x) any(is.na(x)))}
  # if any NA values are found, print a warning
  if(sum(na_params(params)) > 0){
    warning("NA values detected in mediation measures!")
  }
  # identify mp is out of range [0,1]
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
