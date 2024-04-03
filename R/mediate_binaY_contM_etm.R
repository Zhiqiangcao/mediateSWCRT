#' Mediation analysis function for binary outcome and continuous mediator in a stepped wedge design with a time-dependent treatment effect
#'
#' After obtaining parameter estimates from linear mixed effect model for mediator model and generalized linear mixed-effects for outcome model, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP with a time-dependent treatment effect
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
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
#' gammas = c(0.8,1.2,1.6)
#' betas = c(0.60,0.75,0.90)
#' sigma_a = 0.605
#' sigma_tau = 0.334   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=0.625,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=0)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_binaY_contM_etm(data=mydata1)
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
#' res2 = mediate_binaY_contM_etm(data=mydata2,covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_binaY_contM_etm(data=mydata2,covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 12 
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' gammas = c(0.4,0.8,1.2,1.6)
#' betas = c(0.60,0.75,0.90,1.05)
#' set.seed(1234)
#' mydata3 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=0.625,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=0)
#' res4 = mediate_binaY_contM_etm(data=mydata3)
#' print(res4)

mediate_binaY_contM_etm = function(data, id = "id", outcome = "Y", mediator = "M", treatment = "S", cluster = "cluster", 
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
  #fit outcome model and mediator model
  model_Y = summary(glmer(formula_Y, data = data, family=binomial))
  model_M = summary(lmer(formula_M, data = data))
  beta_m_hat = model_Y$coefficients[2*J,1]
  beta_s_hat = model_Y$coefficients[(J+1):(2*J-1),1]
  beta_est = model_Y$coefficients[1:J,1]
  sigmaa = model_Y$varcor
  sigmaa = sqrt(sigmaa$cluster[1])
  gamma_s_hat = model_M$coefficients[(J+1):(2*J-1),1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  sigmaem = model_M$sigma
  #check whether covariates are included in outcome and mediator models
  if(is.null(covariateY)){
    beta_x = 0; ny = 0
  }else{
    ny = length(covariateY)
    beta_x = model_Y$coefficients[(2*J+1):(2*J+ny),1]
    xy = as.matrix(data[,covariateY],ncol=ny)
  }
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    nm = length(covariateM)
    gamma_x = model_M$coefficients[(2*J):(2*J+nm-1),1]
    xm = as.matrix(data[,covariateM],ncol=nm)
  }
  #point estimate
  nie_s = nde_s = numeric(J-1)
  for(s in 1:(J-1)){
    betaas =  beta_s_hat[s]
    gammas = gamma_s_hat[s]
    betajs = beta_est[(s+1):J]
    gammajs = gamma_est[(s+1):J] #corresponding s
    nie_js = nde_js = numeric(J-s)
    for(j in 1:(J-s)){
      if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data$period==(j+s),ncol=ny])
      if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data$period==(j+s),ncol=nm])
      prob_res = cal_prob(betaj=betajs[j],betaa=betaas,betam=beta_m_hat,betax=beta_x,gammaj=gammajs[j],
                          gammaa=gammas,gammax=gamma_x,sigma_em=sigmaem,sigma_tau=sigmatau,sigma_a=sigmaa)
      pr11 = prob_res$pr_11
      pr01 = prob_res$pr_01
      pr00 = prob_res$pr_00
      niej_res = log(pr11/(1-pr11))-log(pr01/(1-pr01))
      ndej_res = log(pr01/(1-pr01))-log(pr00/(1-pr00))
      nie_js[j] = mean(niej_res)  #average all i and k 
      nde_js[j] = mean(ndej_res)  #average all i and k 
    } 
    nie_s[s] = mean(nie_js) #average all j
    nde_s[s] = mean(nde_js) #average all j
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
    model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
    model_M_g = summary(lmer(formula_M, data = data_g))
    beta_m_hat_g = model_Y_g$coefficients[2*J,1]
    beta_s_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
    beta_est_g = model_Y_g$coefficients[1:J,1]
    sigmaa_g = model_Y_g$varcor
    sigmaa_g = sqrt(sigmaa_g$cluster[1])
    
    gamma_s_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
    gamma_est_g = model_M_g$coefficients[1:J,1]
    sigmatau_g = model_M_g$varcor
    sigmatau_g = sqrt(sigmatau_g$cluster[1])
    sigmaem_g = model_M_g$sigma
    if(ny == 0){
      beta_x_g = 0
    }else{
      beta_x_g = model_Y_g$coefficients[(2*J+1):(2*J+ny),1]
      xy = as.matrix(data_g[,covariateY],ncol=ny)
    }
    if(nm == 0){
      gamma_x_g = 0
    }else{
      gamma_x_g = model_M_g$coefficients[(2*J):(2*J+nm-1),1]
      xm = as.matrix(data_g[,covariateM],ncol=nm)
    }
    nie_s_g = nde_s_g = numeric(J-1)
    for(s in 1:(J-1)){
      betaas_g =  beta_s_hat_g[s]
      gammas_g = gamma_s_hat_g[s]
      betajs_g = beta_est_g[(s+1):J]
      gammajs_g = gamma_est_g[(s+1):J] 
      nie_js_g = nde_js_g = numeric(J-s)
      for(j in 1:(J-s)){
        if(ny == 0) xy_j = 0 else xy_j = as.matrix(xy[data_g$period==(j+s),ncol=ny])
        if(nm == 0) xm_j = 0 else xm_j = as.matrix(xm[data_g$period==(j+s),ncol=nm])
        prob_res_g = cal_prob(betaj=betajs_g[j],betaa=betaas_g,betam=beta_m_hat_g,betax=beta_x_g,gammaj=gammajs_g[j],
                              gammaa=gammas_g,gammax=gamma_x_g,sigma_em=sigmaem_g,sigma_tau=sigmatau_g,sigma_a=sigmaa_g)
        pr11_g = prob_res_g$pr_11
        pr01_g = prob_res_g$pr_01
        pr00_g = prob_res_g$pr_00
        niej_res_g = log(pr11_g/(1-pr11_g))-log(pr01_g/(1-pr01_g))
        ndej_res_g = log(pr01_g/(1-pr01_g))-log(pr00_g/(1-pr00_g))
        nie_js_g[j] = mean(niej_res_g)  #average all i and k 
        nde_js_g[j] = mean(ndej_res_g)  #average all i and k 
      }
      nie_s_g[s] = mean(nie_js_g)
      nde_s_g[s] = mean(nde_js_g)
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
