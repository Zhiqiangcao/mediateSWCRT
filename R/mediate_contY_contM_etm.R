#' Mediation analysis function for continuous outcome and continuous mediator in a stepped wedge design with an exposure-time dependent treatment effect
#' 
#' After obtaining parameter estimates from linear mixed-effects models for both outcome and mediator models, this function can obtain 
#' mediation measures including NIE, NDE, TE and MP with an exposure-time dependent treatment effect 
#' 
#'@param data A dataset, which should include outcome, mediator, treatment, cluster, period variables
#'@param outcome The outcome variable in stepped wedge design
#'@param mediator The mediator variable in stepped wedge design
#'@param treatment The treatment variable in stepped wedge design. When there is a implementation period, then corresponding treatment status 
#'       should be set to -1, which is mainly for convenience  
#'@param cluster The cluster variable in stepped wedge design
#'@param period The period variable in stepped wedge design
#'@param exposure The exposure time with J levels (e.g., 0,1,2,...,J-1 when there is no implementation period)
#'@param covariateY A vector of confounders in the outcome regression (default is NULL)
#'@param covariateM A vector of confounders in the mediator regression  (default is NULL)
#'@param a0 The reference treatment level in defining TE and NIE. The default value is 0
#'@param a1 The treatment level of interest in defining TE and NIE. The default value is 1
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP at each exposure time and
#'        overall summary mediation measures, as well as Chi-square test result of TE
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 lmer 
#'@importFrom stats as.formula  
#'@importFrom stats qt 
#'@importFrom stats pchisq 
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
#' eta_e = c(0.5,0.8,1.1) 
#' theta_e = c(0.5,1.0,1.5)
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=0.8,eta_e,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=0)
#' # example 1: mediation analysis without covariates in outcome and mediator models
#' res1 = mediate_contY_contM_etm(data = mydata1)
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
#' res2 = mediate_contY_contM_etm(data = mydata2, covariateY = c("X1"), 
#' covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_contY_contM_etm(data = mydata2, covariateY = c("X1","X2"), 
#' covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' eta_e = c(0.4,0.6,0.8,1.0) 
#' theta_e = c(0.6,0.75,0.9,1.05)
#' set.seed(123456)
#' mydata3 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=0.8,eta_e,sigma_a,sigma_ey,
#' sigma_tau,sigma_em,binary.o=0,binary.m=0)
#' res4 = mediate_contY_contM_etm(data=mydata3)
#' print(res4)
#' 

mediate_contY_contM_etm = function(data, outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
                                    period = "period", exposure = "E", covariateY = NULL, covariateM = NULL, a0 = 0, a1 = 1){
  #first test whether there are "Y","M","A","E","cluster","period" in data
  data = as.data.frame(data)
  require_covariate = c(outcome,mediator,treatment,cluster,period,exposure,covariateY,covariateM)
  all_covariate = colnames(data)
  if(!all(require_covariate %in% all_covariate)){
    stop("Some specified required arguments are not in the data, please check!")
  }
  #check missing variables in interested variables of data.frame
  if(!is.null(covariateY) | !is.null(covariateM)){
    #at least one covariate in covariateY or covariateM
    covariateYM = union(covariateM,covariateY) #extract common covariates
    data = data[,c(outcome,mediator,treatment,cluster,period,exposure,covariateYM)]
  }else{ #in this case, both covariateY = covariateM = NULL
    data = data[,c(outcome,mediator,treatment,cluster,period,exposure)]
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
    formula_Y = as.formula(paste(outcome, "~", period, "+", exposure, 
                                 "+", mediator, "+", "(1|cluster)", sep = ""))
  }
  else {
    formula_Y = as.formula(paste(outcome, "~", period, "+", exposure, "+", mediator,
                                 "+", paste(covariateY, collapse = "+"), 
                                 "+", "(1|cluster)",sep = ""))
  }
  if (is.null(covariateM)) {
    formula_M = as.formula(paste(mediator, "~", period, "+", exposure, 
                                 "+", "(1|cluster)", sep = ""))
  }
  else {
    formula_M = as.formula(paste(mediator, "~",period, "+", exposure, 
                                 "+", paste(covariateM, collapse = "+"),
                                 "+", "(1|cluster)", sep = ""))
  }
  #number of period
  J = length(unique(data[,5]))  #period
  El = length(unique(data[,6])) #maximum value of exposure time
  clusterid = data[,4] #cluster
  clusters = unique(clusterid)
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, exposure)] = lapply(data[,c(period, cluster, exposure)], factor)
  #fit outcome model and mediator model
  model_Y = summary(lmer(formula_Y, data = data))
  model_M = summary(lmer(formula_M, data = data))
  beta_m_hat = model_Y$coefficients[J+El,1]
  theta_e_hat = model_Y$coefficients[(J+1):(J+El-1),1]
  eta_e_hat = model_M$coefficients[(J+1):(J+El-1),1]
  nie_e = beta_m_hat*eta_e_hat*(a1-a0) 
  nde_e = theta_e_hat*(a1-a0) 
  te_e = nie_e + nde_e
  mp_e = nie_e/te_e
  mp_m = mean(nie_e)/mean(te_e)
  #point estimation
  point_est = data.frame(NIE = c(nie_e,mean(nie_e)),NDE = c(nde_e,mean(nde_e)),
                         TE = c(te_e,mean(te_e)), MP = c(mp_e,mp_m))
  row.names(point_est) = c(paste("e=",1:(El-1),sep=""),"overall")
  # warnings 
  # check for NA values in mediation measures
  params = list(nie_e=nie_e,nde_e=nde_e,te_e=te_e,mp_e=mp_e)
  # identify parameters with NA values
  na_params = function(params){sapply(params,function(x) any(is.na(x)))}
  # if any NA values are found, print a warning
  if(sum(na_params(params)) > 0){
    warning("NA values detected in mediation measures!")
  }
  # check mp is out of range [0,1]
  MP = c(mp_e,mp_m)
  if(any(MP < 0) | any(MP > 1)){
    warning("MP is out of range [0,1]!")
  }
  
  #Jackknife variance
  betaG = array(0,dim=c(ng,4,El))
  #Test=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
  te_jack = matrix(0,ng,El-2)
  for(g in 1:ng){
    g_exc = (clusterid != clusters[g])
    data_g = data[g_exc,]
    model_Y_g = summary(lmer(formula_Y, data = data_g))
    model_M_g = summary(lmer(formula_M, data = data_g))
    beta_m_hat_g = model_Y_g$coefficients[J+El,1]
    theta_e_hat_g = model_Y_g$coefficients[(J+1):(J+El-1),1]
    eta_e_hat_g = model_M_g$coefficients[(J+1):(J+El-1),1]
    nie_e_g = beta_m_hat_g*eta_e_hat_g*(a1-a0) 
    nde_e_g = theta_e_hat_g*(a1-a0) 
    te_e_g = nie_e_g + nde_e_g
    mp_e_g = nie_e_g/te_e_g
    est_res_g = cbind(nie_e_g,nde_e_g,te_e_g,mp_e_g)
    #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(El-1))
    te_jack[g,] = te_e_g[1]-te_e_g[-1]
    for(e in 1:El){
      if(e < El){ 
        betaG[g,,e] = est_res_g[e,]
      }else{ #overall NIE, NDE, TE and MP
        overall_m = colMeans(est_res_g)
        overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
        betaG[g,,e] = overall_m
      }
    }
  }
  #check NA in Jackknife estimation
  if(sum(is.na(betaG))>0){
    warning("NA values detected in Jackknife estimations for mediation measures!")
  }
  #Jackknife variance
  var_est = matrix(0,El,4)
  for(j in 1:El){
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
  chi_square_sta_cal = function(est,jack_est,ng,El){
    #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(E-1))
    #est is point estimation of theta(e) at e=1,2,...,E-1
    #jack_est is a ng*(J-2) matrix for T
    point.est = matrix(est[1]-est[-1],ncol=1)  #(E-1)*1 matrix
    jack_est_m = colMeans(jack_est)
    jack_est_matrix = jack_est-matrix(rep(jack_est_m,each=ng),ncol=El-2)
    jack_var_est = ((ng-1)/ng)*(t(jack_est_matrix)%*%jack_est_matrix)
    T.stat = as.numeric(t(point.est)%*%solve(jack_var_est)%*%point.est)
    P.value = 1-pchisq(T.stat, df=El-2)
    resu = c(T.stat,P.value)
    return(resu)
  }
  #test TE(e) can be treated as constant
  chisq_test = chi_square_sta_cal(point_est[-El,3],te_jack,ng,El)
  names(chisq_test) = c("Test.stat","P.value")
  
  return(list(point_est = point_est, var_est = var_est, sd_est = sd_est, 
              ci_lower_confidence_limit = ci_est_low,ci_upper_confidence_limit = ci_est_high,
              chisq_test = chisq_test))
}
