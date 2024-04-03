#' Mediation analysis in a stepped wedge cluster randomized trials
#' 
#' Function for calculating the point and interval estimates for the NIE, NDE, TE and MP in a stepped 
#' wedge cluster randomized trials with constant treatment effect or time-dependent treatment effect
#' 
#'@param data (Required) The name of the dataset.
#'@param method two approximate methods are provided, one is GHQ (Gauss-Hermite Quadrature), the other is AMLE (approximate maximum likelihood estimation)
#'       when both the outcome and mediator are continuous, the method choice is useless; when the outcome is binary and mediator is continuous, we default apply AMLE method.
#'@param id (Required) The id of subjects in dataset
#'@param outcome (Required) Name of the outcome variable, which should be either a continuous or binary datatype.
#'@param mediator (Required) Name of the mediator variable, which should be either a continuous or binary datatype.
#'@param treatment (Required) Name of the treatment variable, which should be either a continuous or binary datatype. For constant treatment effect, default notation is A; for time-dependent treatment effect, default notation is S.  
#'@param cluster (Required) Name of the cluster variable, which should be a factor
#'@param period (Required) Name of the period variable, which should be a factor
#'@param covariate.outcome A vector of names showing the confounding variables used in the outcome model. The default value 
#'       is NULL, which represents no confounding variables. We only accepted continuous and binary confounding variables, 
#'       if one confounding variable is categorical, please set it to a series of binary variables in advance.
#'@param covariate.mediator A vector of names showing the confounding variables used in the mediator model. The default value
#'       is NULL, which represents no confounding variables. We only accepted continuous and binary confounding variables, if one
#'       confounding variable is categorical, please set it to a series of binary variables in advance.
#'@param a0 The baseline treatment level (i.e., $a^*$). The default value is 0.
#'@param a1 The new treatment level (i.e., $a$). The default value is 1.
#'@param binary.outcome (Required) If the outcome is binary, set to 1. If the outcome is continuous, set to 0. 
#'@param binary.mediator (Required) If the mediator is binary, set to 1. If the mediator is continuous, set to 0.
#'@param time.dependent (Required) If the treatment effect is time-dependent, set to 1. If the treatment effect is constant, set to 0.
#'
#'@return A list containing point and interval estimates of NIE, NDE, TE and MP.
#'@export
#'
#'@author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#'@references Cao Z. and Fan Li. Analysis of stepped wedge cluster randomized trials in the. Statistics in Medicine. 2024;0(0):1-25.
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
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata1 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=0.625,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' # example 1: constant treatment effect, mediation analysis of continuous outcome 
#' # and continuous mediator without covariates
#' res1 = mediate_swcrt(data=mydata1)
#' print(res1)
#' 
#' # example 2: constant treatment effect, mediation analysis of binary outcome 
#' # and continuous mediator without covariates
#' sigma_a = 0.605
#' sigma_tau = 0.334
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata2 = gen_data_hhm(n,J,m,beta,gamma,betaA=0.75,betaM=1,gammaA=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=1,binary.mediator=0)
#' res2 = mediate_swcrt(data=mydata2,binary.outcome = 1)
#' print(res2)
#' 
#' # example 3: time-dependent treatment effect, mediation analysis of continuous outcome 
#' # and continuous mediator without covariates
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' gammas = c(0.5,0.8,1.1) 
#' betas = c(0.5,1.0,1.5)
#' set.seed(123456)
#' mydata3 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=0.8,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=0)
#' res3 = mediate_swcrt(data=mydata3,treatment = "S", time.dependent = 1)
#' print(res3)
#' 
#' # example 4: time-dependent treatment effect, mediation analysis of continuous outcome
#' # and binary mediator without covariates
#' gammas = c(0.5,1.3,2.1) 
#' betas = c(0.6,1,1.4)
#' sigma_a = 0.334
#' sigma_tau = 0.605   
#' sigma_em = sigma_ey = 1
#' set.seed(123456)
#' mydata4 = gen_data_etm(n,J,m,beta,gamma,betas,betaM=1.2,gammas,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.outcome=0,binary.mediator=1)
#' # using default AMLE method
#' res4 = mediate_swcrt(data=mydata4,treatment = "S", 
#' binary.mediator = 1, time.dependent = 1)
#' print(res4)
#' # using GHQ method
#' res5 = mediate_swcrt(data=mydata4,method = "GHQ",treatment = "S", 
#' binary.mediator = 1, time.dependent = 1)
#' print(res5)
#' 
mediate_swcrt = function(data, method = "AMLE", id = "id", outcome = "Y", mediator = "M", treatment = "A", cluster = "cluster", 
                         period = "period", covariate.outcome = NULL, covariate.mediator = NULL, a0 = 0, a1 = 1,
                         binary.outcome = 0, binary.mediator = 0, time.dependent = 0){
  data = as.data.frame(data)
  covariateY = covariate.outcome
  covariateM = covariate.mediator
  #constant treatment effect
  if(time.dependent == 0){ #constant treatment effect
    if (binary.outcome == 0 & binary.mediator == 0){ #type 1: continuous outcome and continuous mediator
      delta_res = mediate_contY_contM_hhm(data = data, id = id, outcome = outcome, mediator = mediator, treatment = treatment, cluster = cluster, 
                                          period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else if(binary.outcome == 0 & binary.mediator == 1){ #type 2: continuous outcome and binary mediator
      delta_res = mediate_contY_binaM_hhm(data = data, method = method, id = id, outcome = outcome, mediator = mediator, treatment = treatment, 
                                          cluster = cluster, period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else if(binary.outcome == 1 & binary.mediator == 0){ #type 3: binary outcome and continuous mediator
      delta_res = mediate_binaY_contM_hhm(data = data, id = id, outcome = outcome, mediator = mediator, treatment = treatment, cluster = cluster, 
                                         period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else{#type 4: binary outcome and binary mediator
      delta_res = mediate_binaY_binaM_hhm(data = data, method = method, id = id, outcome = outcome, mediator = mediator, treatment = treatment, 
                                          cluster = cluster, period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }
  }
  if(time.dependent == 1){ #time-dependent treatment effect
    if (binary.outcome == 0 & binary.mediator == 0){ #type 1: continuous outcome and continuous mediator
      delta_res = mediate_contY_contM_etm(data = data, id = id, outcome = outcome, mediator = mediator, treatment = treatment, cluster = cluster, 
                                          period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else if(binary.outcome == 0 & binary.mediator == 1){ #type 2: continuous outcome and binary mediator
      delta_res = mediate_contY_binaM_etm(data = data, method = method, id = id, outcome = outcome, mediator = mediator, treatment = treatment, 
                                          cluster = cluster, period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else if(binary.outcome == 1 & binary.mediator == 0){ #type 3: binary outcome and continuous mediator
      delta_res = mediate_binaY_contM_etm(data = data, id = id, outcome = outcome, mediator = mediator, treatment = treatment, cluster = cluster, 
                                          period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }else{#type 4: binary outcome and binary mediator
      delta_res = mediate_binaY_binaM_etm(data = data, method = method, id = id, outcome = outcome, mediator = mediator, treatment = treatment, 
                                          cluster = cluster, period = period, covariateY = covariateY, covariateM = covariateM, a0 = a0, a1 = a1)
    }
  }
  return(delta_res)
}
