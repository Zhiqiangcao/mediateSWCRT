#' Mediation analysis function for binary outcome and continuous mediator in a stepped wedge design
#' with an exposure-time dependent treatment effect
#' 
#' After obtaining parameter estimates from linear mixed effect model for mediator model and 
#' generalized linear mixed effect model for outcome model, this function can obtain 
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
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP at each exposure time and overall summary mediation measures, as well as Chi-square test result of TE
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'
#'@importFrom lme4 lmer   
#'@importFrom lme4 glmer
#'@importFrom stats as.formula  
#'@importFrom stats binomial
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
#' eta_e = c(0.8,1.2,1.6)
#' theta_e = c(0.60,0.75,0.90)
#' sigma_a = 0.605
#' sigma_tau = 0.334   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=1,binary.m=0)
#' # example 1: mediation analysis without covariates in outcome and mediator models
#' res1 = mediate_binaY_contM_etm(data=mydata1)
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
#' res2 = mediate_binaY_contM_etm(data=mydata2,covariateY = c("X1"), covariateM = c("X2"))
#' print(res2)
#' 
#' # example 3: mediation analysis with the same covariates in outcome and mediator models
#' res3 = mediate_binaY_contM_etm(data=mydata2,covariateY = c("X1","X2"), covariateM = c("X1","X2"))
#' print(res3)
#' 
#' # example 4: mediation analysis with different clusters and periods
#' I = 12 
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' eta_e = c(0.4,0.8,1.2,1.6)
#' theta_e = c(0.60,0.75,0.90,1.05)
#' set.seed(123456)
#' mydata3 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=1,binary.m=0)
#' res4 = mediate_binaY_contM_etm(data=mydata3)
#' print(res4)

mediate_binaY_contM_etm = function(data, outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
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
  #compute double integral \mu (STA method)
  cal_prob = function(betaj,theta,betam,betax,gammaj,eta,gammax,sigma_em,sigma_tau,sigma_a){
    nu_11 = exp(betaj+theta*a1+betam*(gammaj+eta*a1+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a1+betam*(gammaj+eta*a1+sum(xm_j*gammax))+sum(xy_j*betax)))
    p1_11 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_11+(2*nu_11^3-3*nu_11^2+nu_11)*0.5*sigma_a^2)
    p2_11 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_11^2+(6*nu_11^4-10*nu_11^3+4*nu_11^2)*0.5*sigma_a^2)
    p3_11 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_11^3+(12*nu_11^5-21*nu_11^4+9*nu_11^3)*0.5*sigma_a^2)
    pr_11 = p1_11+p2_11+p3_11
    #
    nu_01 = exp(betaj+theta*a1+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a1+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax)))
    p1_01 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_01+(2*nu_01^3-3*nu_01^2+nu_01)*0.5*sigma_a^2)
    p2_01 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_01^2+(6*nu_01^4-10*nu_01^3+4*nu_01^2)*0.5*sigma_a^2)
    p3_01 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_01^3+(12*nu_01^5-21*nu_01^4+9*nu_01^3)*0.5*sigma_a^2)
    pr_01 = p1_01+p2_01+p3_01
    #
    nu_00 = exp(betaj+theta*a0+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax))/(1+exp(betaj+theta*a0+betam*(gammaj+eta*a0+sum(xm_j*gammax))+sum(xy_j*betax)))
    p1_00 = (1+0.5*betam^2*(sigma_em^2+sigma_tau^2))*(nu_00+(2*nu_00^3-3*nu_00^2+nu_00)*0.5*sigma_a^2)
    p2_00 = -3/2*betam^2*(sigma_em^2+sigma_tau^2)*(nu_00^2+(6*nu_00^4-10*nu_00^3+4*nu_00^2)*0.5*sigma_a^2)
    p3_00 = betam^2*(sigma_em^2+sigma_tau^2)*(nu_00^3+(12*nu_00^5-21*nu_00^4+9*nu_00^3)*0.5*sigma_a^2)
    pr_00 = p1_00+p2_00+p3_00
    res = data.frame(pr_11,pr_01,pr_00)
    return(res)
  }
  
  #number of period
  J = length(unique(data[,5]))  #period
  El = length(unique(data[,6])) #maximum value of exposure time
  clusterid = data[,4] #cluster
  clusters = unique(clusterid)
  dje = J-El  #if there are no buffer period, then J=El and dje=0
  #number of clusters
  ng = length(clusters)
  data[,c(period, cluster, exposure)] = lapply(data[,c(period, cluster, exposure)], factor)
  #fit outcome model 
  model_Y = summary(glmer(formula_Y, data = data, family=binomial))
  beta_m_hat = model_Y$coefficients[J+El,1]
  theta_e_hat = model_Y$coefficients[(J+1):(J+El-1),1]
  beta_est = model_Y$coefficients[1:J,1]
  sigmaa = model_Y$varcor
  sigmaa = sqrt(sigmaa$cluster[1])
  #fit mediator model
  model_M = summary(lmer(formula_M, data = data))
  eta_e_hat = model_M$coefficients[(J+1):(J+El-1),1]
  gamma_est = model_M$coefficients[1:J,1]
  sigmatau = model_M$varcor
  sigmatau = sqrt(sigmatau$cluster[1])
  sigmaem = model_M$sigma
  if(is.null(covariateY)){
    beta_x = 0; ny = 0
  }else{
    #in case fixed-effect model matrix is rank deficient
    beta_x = model_Y$coefficients[-(1:(J+El)),1]
    xy = as.matrix(data[,names(beta_x)],ncol=length(names(beta_x)))
    ny = dim(xy)[2]
  }
  if(is.null(covariateM)){
    gamma_x = 0; nm = 0
  }else{
    gamma_x = model_M$coefficients[-(1:(J+El-1)),1]
    xm = as.matrix(data[,names(gamma_x)],ncol=length(names(gamma_x)))
    nm = dim(xm)[2]
  }
  #point estimate
  nie_e = nde_e = numeric(El-1)
  for(e in 1:(El-1)){
    thetae =  theta_e_hat[e]
    etae = eta_e_hat[e]
    betaje = beta_est[(e+1+dje):J]
    gammaje = gamma_est[(e+1+dje):J] #corresponding to e
    nie_je = nde_je = numeric(J-e-dje)
    for(j in 1:(J-e-dje)){
      if(ny == 0) xy_j = 0 else {
        xy_j = as.matrix(xy[data$period==(j+e+dje),],ncol=ny)
        #modify, using median
        xy_j = apply(xy_j,2,median)
      }
      if(nm == 0) xm_j = 0 else {
        xm_j = as.matrix(xm[data$period==(j+e+dje),],ncol=nm)
        #modify, using median
        xm_j = apply(xm_j,2,median)
      }
      prob_res = cal_prob(betaj=betaje[j],theta=thetae,betam=beta_m_hat,betax=beta_x,gammaj=gammaje[j],
                          eta=etae,gammax=gamma_x,sigma_em=sigmaem,sigma_tau=sigmatau,sigma_a=sigmaa)
      pr11 = prob_res$pr_11
      pr01 = prob_res$pr_01
      pr00 = prob_res$pr_00
      nie_je[j] = log(pr11/(1-pr11))-log(pr01/(1-pr01))
      nde_je[j] = log(pr01/(1-pr01))-log(pr00/(1-pr00))
    } 
    nie_e[e] = mean(nie_je) #average all j
    nde_e[e] = mean(nde_je) #average all j
  }
  te_e = nie_e + nde_e
  mp_e = nie_e/te_e
  mp_m = mean(nie_e)/mean(te_e)
  #point estimation
  point_est = data.frame(NIE = c(nie_e,mean(nie_e)),NDE = c(nde_e,mean(nde_e)),
                         TE = c(te_e,mean(te_e)), MP = c(mp_e,mp_m))
  row.names(point_est) = c(paste("e=",1:(El-1),sep=""),"overall")
  #Jackknife variance
  betaG = array(0,dim=c(ng,4,El))
  #Test=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(J-1))
  te_jack = matrix(0,ng,El-2)
  for(g in 1:ng){
    g_exc = (clusterid != clusters[g])
    data_g = data[g_exc,]
    #fit outcome model
    model_Y_g = summary(glmer(formula_Y, data = data_g, family=binomial))
    beta_m_hat_g = model_Y_g$coefficients[J+El,1]
    theta_e_hat_g = model_Y_g$coefficients[(J+1):(J+El-1),1]
    beta_est_g = model_Y_g$coefficients[1:J,1]
    sigmaa_g = model_Y_g$varcor
    sigmaa_g = sqrt(sigmaa_g$cluster[1])
    #fit mediator model
    model_M_g = summary(lmer(formula_M, data = data_g))
    eta_e_hat_g = model_M_g$coefficients[(J+1):(J+El-1),1]
    gamma_est_g = model_M_g$coefficients[1:J,1]
    sigmatau_g = model_M_g$varcor
    sigmatau_g = sqrt(sigmatau_g$cluster[1])
    sigmaem_g = model_M_g$sigma
    if(ny == 0){
      beta_x_g = 0
    }else{
      beta_x_g = model_Y_g$coefficients[-(1:(J+El)),1]
      xy = as.matrix(data_g[,names(beta_x_g)],ncol=length(names(beta_x_g)))
      ny = dim(xy)[2]
    }
    if(nm == 0){
      gamma_x_g = 0
    }else{
      gamma_x_g = model_M_g$coefficients[-(1:(J+El-1)),1]
      xm = as.matrix(data_g[,names(gamma_x_g)],ncol=length(names(gamma_x_g)))
      nm = dim(xm)[2]
    }
    nie_e_g = nde_e_g = numeric(El-1)
    for(e in 1:(El-1)){
      thetae_g =  theta_e_hat_g[e]
      etae_g = eta_e_hat_g[e]
      betaje_g = beta_est_g[(e+1+dje):J]
      gammaje_g = gamma_est_g[(e+1+dje):J] 
      nie_je_g = nde_je_g = numeric(J-e-dje)
      for(j in 1:(J-e-dje)){
        if(ny == 0) xy_j = 0 else {
          xy_j = as.matrix(xy[data_g$period==(j+e+dje),],ncol=ny)
          #modify, using median
          xy_j = apply(xy_j,2,median)
        }
        if(nm == 0) xm_j = 0 else {
          xm_j = as.matrix(xm[data_g$period==(j+e+dje),],ncol=nm)
          #modify, using median
          xm_j = apply(xm_j,2,median)
        }
        prob_res_g = cal_prob(betaj=betaje_g[j],theta=thetae_g,betam=beta_m_hat_g,betax=beta_x_g,gammaj=gammaje_g[j],
                              eta=etae_g,gammax=gamma_x_g,sigma_em=sigmaem_g,sigma_tau=sigmatau_g,sigma_a=sigmaa_g)
        pr11_g = prob_res_g$pr_11
        pr01_g = prob_res_g$pr_01
        pr00_g = prob_res_g$pr_00
        nie_je_g[j] = log(pr11_g/(1-pr11_g))-log(pr01_g/(1-pr01_g))
        nde_je_g[j] = log(pr01_g/(1-pr01_g))-log(pr00_g/(1-pr00_g))
      }
      nie_e_g[e] = mean(nie_je_g)
      nde_e_g[e] = mean(nde_je_g)
    }
    te_e_g = nie_e_g + nde_e_g
    mp_e_g = nie_e_g/te_e_g
    est_res_g = cbind(nie_e_g,nde_e_g,te_e_g,mp_e_g)
    #T=(theta(1)-theta(2),theta(1)-theta(3),...,theta(1)-theta(E-1))
    te_jack[g,] = te_e_g[1]-te_e_g[-1]
    for(e in 1:El){
      if(e < El){
        betaG[g,,e] = est_res_g[e,]
      }else{
        #overall NIE, NDE, TE and MP
        overall_m = colMeans(est_res_g)
        overall_m[4] = overall_m[1]/overall_m[3]  #using mean(nie)/mean(te) as mean of mp
        betaG[g,,e] = overall_m
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
