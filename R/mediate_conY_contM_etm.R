#' Mediation analysis function for continuous outcome and continuous mediator in a stepped wedge design with a time-dependent treatment effect
#'
#' After obtaining parameter estimates from linear mixed-effects models for both outcome and mediator models, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP with a time-dependent treatment effect 
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param id The id of subjects in stepped wedge design
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design, here treatment is a factor variable with J levels
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression  (default is NULL)
#'@param a0 The baseline treatment level
#'@param a1 The new treatment level
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP at each exposure time, as well as the average mediation measures
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 lmer 
#'@importFrom stats as.formula  
#'@importFrom stats qt 
#'
#'@examples 
#' library(lme4)
#' n = 15 
#' J = 4
#' m = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' gammas = c(0.5,0.8,1.1) 
#' betas = c(0.5,1.0,1.5)
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=0.8,gammas,sigma_a,sigma_ey,
#' sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_contY_contM_etm(data=mydata1)
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
#' mydata2 = data.frame(mydata1,covdata)
#' res2 = mediate_contY_contM_etm(data=mydata2,covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_contY_contM_etm(data=mydata2,covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' gammas = c(0.4,0.6,0.8,1.0) 
#' betas = c(0.6,0.75,0.9,1.05)
#' set.seed(1234)
#' mydata3 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=0.8,gammas,sigma_a,sigma_ey,
#' sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' res4 = mediate_contY_contM_etm(data=mydata3)
#' print(res4)
#' 


mediate_contY_contM_etm = function(data, id = "id", outcome = "Y", mediator = "M", treatment = "S", cluster = "cluster", 
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
  #number of period
  J = max(data$period)
  clusterid = data$cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, treatment, id)] = lapply(data[,c(period, cluster, treatment, id)], factor)
  #fit outcome model and mediator model
  model_Y = summary(lmer(formula_Y, data = data))
  model_M = summary(lmer(formula_M, data = data))
  beta_m_hat = model_Y$coefficients[2*J,1]
  beta_s_hat = model_Y$coefficients[(J+1):(2*J-1),1]
  gamma_s_hat = model_M$coefficients[(J+1):(2*J-1),1]
  nie_s = beta_m_hat*gamma_s_hat*(a1-a0) 
  nde_s = beta_s_hat*(a1-a0) 
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
    g_exc <- (clusterid != clusters[g])
    data_g = data[g_exc,]
    model_Y_g = summary(lmer(formula_Y, data = data_g))
    model_M_g = summary(lmer(formula_M, data = data_g))
    beta_m_hat_g = model_Y_g$coefficients[2*J,1]
    beta_s_hat_g = model_Y_g$coefficients[(J+1):(2*J-1),1]
    gamma_s_hat_g = model_M_g$coefficients[(J+1):(2*J-1),1]
    nie_s_g = beta_m_hat_g*gamma_s_hat_g*(a1-a0) 
    nde_s_g = beta_s_hat_g*(a1-a0) 
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