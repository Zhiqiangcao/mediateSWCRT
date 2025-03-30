#' Mediation analysis function for continuous outcome and continuous mediator in a stepped wedge design 
#' based on nested exchangeable correlation model with a constant treatment effect 
#' 
#' After obtaining parameter estimates from linear mixed-effects models for both outcome and mediator models, 
#' this function can obtain mediation measures including NIE, NDE, TE and MP assuming a constant treatment effect 
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param outcome The outcome variable in stepped wedge design (continuous)
#'@param mediator The mediator variable in stepped wedge design (continuous)
#'@param treatment The treatment variable in stepped wedge design. When there is a implementation period, then corresponding treatment status 
#'       should be set to -1, which is mainly for convenience 
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param exposure The exposure time with J levels (e.g., 0,1,2,...,J-1 when there is no implementation period). If dataset has no exposure variable, then set it to NULL
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
#'@importFrom lme4 lmer  
#'@importFrom stats as.formula  
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
#' mydata1 = gen_data_nem(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=0)
#' # example 1: mediation analysis without covariates in outcome and mediator models
#' res1 = mediate_contY_contM_nem(data=mydata1)
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
#' mydata2 = data.frame(mydata1,covdata)
#' res2 = mediate_contY_contM_nem(data = mydata2,covariateY = c("X1"), 
#' covariateM = c("X2"))
#' print(res2)  
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_contY_contM_nem(data = mydata2,covariateY = c("X1","X2"), 
#' covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' set.seed(123456)
#' mydata3 = gen_data_nem(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=0)
#' res4 = mediate_contY_contM_nem(data=mydata3)
#' print(res4)

mediate_contY_contM_nem = function(data, outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
                                   period = "period", exposure = "E", covariateY = NULL, covariateM = NULL, a0 = 0, a1 = 1){
  #first test whether there are "Y","M","A","cluster","period" in data
  data = as.data.frame(data)
  require_covariate = c(outcome,mediator,treatment,cluster,period,covariateY,covariateM)
  all_covariate = colnames(data)
  if(!all(require_covariate %in% all_covariate)){
    stop("Some specified required arguments are not in the data, please check!")
  }
  #check missing variables in required variables of data.frame
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
                                 "+", mediator, "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  else {
    formula_Y = as.formula(paste(outcome, "~", period, "+", treatment, "+", mediator,
                                 "+", paste(covariateY, collapse = "+"), 
                                 "+", "(1|cluster)", "+", "(1|period:cluster)",sep = ""))
  }
  if (is.null(covariateM)) {
    formula_M = as.formula(paste(mediator, "~", period, "+", treatment, 
                                 "+", "(1|cluster)", "+", "(1|period:cluster)", sep = ""))
  }
  else {
    formula_M = as.formula(paste(mediator, "~",period, "+", treatment, 
                                 "+", paste(covariateM, collapse = "+"),
                                 "+", "(1|cluster)", "+", "(1|period:cluster)",  sep = ""))
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
  model_M = summary(lmer(formula_M, data = data))
  theta_hat = model_Y$coefficients[J+1,1]
  beta_m_hat = model_Y$coefficients[J+2,1]
  eta_hat = model_M$coefficients[J+1,1]
  #point estimation
  nie = beta_m_hat*eta_hat*(a1-a0) 
  nde = theta_hat*(a1-a0) 
  te = nie + nde
  mp = nie/te
  point_est = c(nie,nde,te,mp)
  names(point_est) = c("NIE", "NDE", "TE", "MP")
  # warnings 
  # check for NA values in mediation measures
  params = list(nie=nie,nde=nde,te=te,mp=mp)
  # identify parameters with NA values
  na_params = function(params){sapply(params,function(x) any(is.na(x)))}
  # if any NA values are found, print a warning
  if(sum(na_params(params)) > 0){
    warning("NA values detected in mediation measures!")
  }
  # check mp is out of range [0,1] or not
  if(mp < 0 | mp > 1){
    warning("MP is out of range [0,1]!")
  }
  #Jackknife variance
  betaG = matrix(0,ng,4)
  for(g in 1:ng){
    g_exc = (clusterid != clusters[g])
    data_g = data[g_exc,]
    model_Y_g = summary(lmer(formula_Y, data = data_g))
    model_M_g = summary(lmer(formula_M, data = data_g))
    beta_m_hat_g = model_Y_g$coefficients[J+2,1]
    theta_hat_g = model_Y_g$coefficients[J+1,1]
    eta_hat_g = model_M_g$coefficients[J+1,1]
    betaG[g,1] = beta_m_hat_g*eta_hat_g*(a1-a0)  #NIE
    betaG[g,2] = theta_hat_g*(a1-a0)      #NDE
    betaG[g,3] = betaG[g,1]+betaG[g,2]   #TE
    betaG[g,4] = betaG[g,1]/betaG[g,3]   #MP
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


