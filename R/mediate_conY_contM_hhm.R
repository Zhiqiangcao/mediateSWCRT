#' Mediation analysis function for continuous outcome and continuous mediator in a stepped wedge design based on Hussey and Hughes model
#'
#' After obtaining parameter estimates from linear mixed-effects models for both outcome and mediator models, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect 
#' 
#'@param data A dataset, which should include id, outcome, mediator, treatment, cluster, period variables
#'@param id The id of subjects in stepped wedge design
#'@param outcome The outcome variable in stepped wedge design (continuous)
#'@param mediator The mediator variable in stepped wedge design (continuous)
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
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=0.625,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' # example 1: mediation analysis without covariates in both outcome and mediator models
#' res1 = mediate_contY_contM_hhm(data=mydata1,covariateY = NULL, covariateM = NULL)
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
#' res2 = mediate_contY_contM_hhm(data=mydata2,covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)  
#' 
#' # example 3: mediation analysis with the same covariates in both outcome and mediator models
#' res3 = mediate_contY_contM_hhm(data=mydata2,covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' n = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(1234)
#' mydata3 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=0.625,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' res4 = mediate_contY_contM_hhm(data=mydata3,covariateY = NULL, covariateM = NULL)
#' print(res4)

mediate_contY_contM_hhm = function(data, id = "id", outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
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
  beta_a_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  gamma_a_hat = model_M$coefficients[J+1,1]
  #point estimation
  nie = beta_m_hat*gamma_a_hat*(a1-a0) 
  nde = beta_a_hat*(a1-a0) 
  te = nie + nde
  mp = nie/te
  point_est = c(nie,nde,te,mp)
  names(point_est) = c("NIE", "NDE", "TE", "MP")
  #Jacknife variance
  betaG = matrix(0,ng,4)
  for(g in 1:ng){
    g_exc <- (clusterid != clusters[g])
    data_g = data[g_exc,]
    model_Y_g = summary(lmer(formula_Y, data = data_g))
    model_M_g = summary(lmer(formula_M, data = data_g))
    beta_m_hat_g = model_Y_g$coefficients[J+2,1]
    beta_a_hat_g = model_Y_g$coefficients[J+1,1]
    gamma_a_hat_g = model_M_g$coefficients[J+1,1]
    betaG[g,1] = beta_m_hat_g*gamma_a_hat_g*(a1-a0)  #NIE
    betaG[g,2] = beta_a_hat_g*(a1-a0)      #NDE
    betaG[g,3] = betaG[g,1]+betaG[g,2]   #TE
    betaG[g,4] = betaG[g,1]/betaG[g,3]   #MP
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


