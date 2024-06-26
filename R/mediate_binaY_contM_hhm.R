#' Mediation analysis function for binary outcome and continuous mediator in a stepped wedge design based on Hussey and Hughes model
#'
#' After obtaining parameter estimates from linear mixed effect model for mediator model and generalized linear mixed-effects for outcome model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
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
#'@importFrom lme4 lmer   
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
#' sigma_tau = 0.334   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=1,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=0)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_binaY_contM_hhm(data=mydata1)
#' print(res1)
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
#' res2 = mediate_binaY_contM_hhm(data=mydata2,covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_binaY_contM_hhm(data=mydata2,covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 12 
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(1234)
#' mydata3 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=1,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=0)
#' res4 = mediate_binaY_contM_hhm(data=mydata3)
#' print(res4)

mediate_binaY_contM_hhm = function(data, id = "id", outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
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
  #compute double integral s (AMLE method)
  cal_prob = function(betaj,betaa,betam,betax,gammaj,gammaa,gammax,sigma_em,sigma_tau,sigma_a){
    nu_11 = exp(betaj+betaa*a1+betam*(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax))/(1+
                                                                                                               exp(betaj+betaa*a1+betam*(gammaj+gammaa*a1+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax)))
    p1_11 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_11+(2*nu_11^3-3*nu_11^2+nu_11)*0.5*sigma_a^2)
    p2_11 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_11^2+(6*nu_11^4-10*nu_11^3+4*nu_11^2)*0.5*sigma_a^2)
    p3_11 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_11^3+(12*nu_11^5-21*nu_11^4+9*nu_11^3)*0.5*sigma_a^2)
    pr_11 = p1_11+p2_11+p3_11
    #
    nu_01 = exp(betaj+betaa*a1+betam*(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax))/(1+
                                                                                                               exp(betaj+betaa*a1+betam*(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax)))
    p1_01 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_01+(2*nu_01^3-3*nu_01^2+nu_01)*0.5*sigma_a^2)
    p2_01 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_01^2+(6*nu_01^4-10*nu_01^3+4*nu_01^2)*0.5*sigma_a^2)
    p3_01 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_01^3+(12*nu_01^5-21*nu_01^4+9*nu_01^3)*0.5*sigma_a^2)
    pr_01 = p1_01+p2_01+p3_01
    #
    nu_00 = exp(betaj+betaa*a0+betam*(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax))/(1+
                                                                                                               exp(betaj+betaa*a0+betam*(gammaj+gammaa*a0+as.numeric(xm_j%*%gammax))+as.numeric(xy_j%*%betax)))
    p1_00 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_00+(2*nu_00^3-3*nu_00^2+nu_00)*0.5*sigma_a^2)
    p2_00 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_00^2+(6*nu_00^4-10*nu_00^3+4*nu_00^2)*0.5*sigma_a^2)
    p3_00 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_00^3+(12*nu_00^5-21*nu_00^4+9*nu_00^3)*0.5*sigma_a^2)
    pr_00 = p1_00+p2_00+p3_00
    res = data.frame(pr_11,pr_01,pr_00)
    return(res)
  }
  #number of period
  J = max(data$period)
  clusterid = data$cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment, id)] = lapply(data[,c(period, cluster, treatment, id)], factor)
  #fit outcome model
  model_Y = summary(glmer(formula_Y, data = data, family=binomial))
  beta_a_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  beta_est = model_Y$coefficients[1:J,1]  
  sigmaa = model_Y$varcor
  sigmaa = sqrt(sigmaa$cluster[1])
  #fit mediator model
  model_M = summary(lmer(formula_M, data = data))
  gamma_a_hat = model_M$coefficients[J+1,1]
  gamma_est = model_M$coefficients[1:J,1] #\gamma_{0j} for j=1,2,3,4,5
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  sigmaem = model_M$sigma
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
  #point estimate
  nie_j = nde_j = numeric(J) 
  for(j in 1:J){
    gammaj = gamma_est[j];
    betaj = beta_est[j];
    if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==j,],ncol=nm)
    if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data$period==j,],ncol=ny)
    ##approximate method
    prob_res = cal_prob(betaj=betaj,betaa=beta_a_hat,betam=beta_m_hat,betax=beta_x,gammaj=gammaj,gammaa=gamma_a_hat,
                        gammax = gamma_x,sigma_em=sigmaem,sigma_tau=sigmatau,sigma_a=sigmaa)
    pr11 = prob_res$pr_11
    pr01 = prob_res$pr_01
    pr00 = prob_res$pr_00
    niej_res = log(pr11/(1-pr11))-log(pr01/(1-pr01))
    ndej_res = log(pr01/(1-pr01))-log(pr00/(1-pr00))
    nie_j[j] = mean(niej_res)
    nde_j[j] = mean(ndej_res)
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
    model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
    beta_a_hat_g = model_Y_g$coefficients[J+1,1]
    beta_m_hat_g = model_Y_g$coefficients[J+2,1]
    beta_est_g = model_Y_g$coefficients[1:J,1]  
    sigmaa_g = model_Y_g$varcor
    sigmaa_g = sqrt(sigmaa_g$cluster[1])
    #fit mediator model
    model_M_g = summary(lmer(formula_M, data = data_g))
    gamma_a_hat_g = model_M_g$coefficients[J+1,1]
    gamma_est_g = model_M_g$coefficients[1:J,1] #\gamma_{0j} for j=1,2,3,4,5
    sigmatau_g = model_M_g$varcor
    sigmatau_g = sqrt(sigmatau_g$cluster[1])
    sigmaem_g = model_M_g$sigma
    #check whether covariates are included in outcome and mediator models
    if(ny == 0){
      beta_x_g = 0;
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
      gammaj_g = gamma_est_g[j];
      betaj_g = beta_est_g[j];
      if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==j,],ncol=nm)
      if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data_g$period==j,],ncol=ny)
      ##approximate method
      prob_res_g = cal_prob(betaj=betaj_g,betaa=beta_a_hat_g,betam=beta_m_hat_g,betax=beta_x_g,gammaj=gammaj_g,
                            gammaa=gamma_a_hat_g,gammax = gamma_x_g,sigma_em=sigmaem_g,sigma_tau=sigmatau_g,sigma_a=sigmaa_g)
      pr11_g = prob_res_g$pr_11
      pr01_g = prob_res_g$pr_01
      pr00_g = prob_res_g$pr_00
      niej_res_g = log(pr11_g/(1-pr11_g))-log(pr01_g/(1-pr01_g))
      ndej_res_g = log(pr01_g/(1-pr01_g))-log(pr00_g/(1-pr00_g))
      nie_jg[j] = mean(niej_res_g)
      nde_jg[j] = mean(ndej_res_g)
    }
    te_jg = nie_jg+nde_jg
    mp_mg = mean(nie_jg)/mean(te_jg)
    point_est_g = c(mean(nie_jg),mean(nde_jg),mean(te_jg),mp_mg)
    betaG[g,] = point_est_g
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
