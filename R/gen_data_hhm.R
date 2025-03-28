#' Generate data in a stepped wedge design based on Hussey and Hughes model with a constant treatment effect
#' 
#' This function can generate longitudinal data in a stepped wedge design based on Hussey and Hughes model
#' with a constant treatment effect.
#' 
#' @param I The number of cluster
#' @param J The number of time periods
#' @param n The number of individuals per each cluster
#' @param beta The underling time trend effect in a outcome model
#' @param gamma The underling time trend effect in a mediator model
#' @param theta The exposure effect on the outcome conditional on mediator and other possible confounders
#' @param beta_M The mediator effect on the outcome conditional on other possible confounders
#' @param eta The exposure effect on the mediator conditional on other possible confounders
#' @param sigma_a The standard error of random effect in outcome model (assuming it follows normal distribution with mean zero)
#' @param sigma_ey The standard error of error term in outcome model (assuming it follows normal distribution with mean zero)
#' @param sigma_tau The standard error of random effect in mediator model (assuming it follows normal distribution with mean zero)
#' @param sigma_em The standard error of error term in mediator model (assuming it follows normal distribution with mean zero)
#' @param binary.o (Required) If the outcome is binary, set to 1. If the outcome is continuous, set to 0
#' @param binary.m (Required) If the mediator is binary, set to 1. If the mediator is continuous, set to 0
#' 
#' @return A data frame including cluster, period, id, E (time since intervention, i.e., exposure time), A (treatment status indicator),
#'         c (the start time of the treatment), alpha (random effects of observations within a cluster), tau (random effects of mediator within a cluster),
#'         M (mediator), Y (outcome)
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#' @references Hussey M.A. and Hughes J.P. Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials.
#' 2007;28(2):182-191.
#' 
#' Kenny A., Voldal E.C., Xia F., Heagerty P.J. and Hughes J.P. Analysis of stepped wedge cluster randomized trials in the
#' presence of a time-varying treatment effect. Statistics in Medicine. 2022;41(22):4311-4339.
#' 
#' @importFrom stats rnorm 
#' @importFrom stats rbinom
#' 
#' @examples
#' I = 9  
#' J = 4
#' n = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' # generate continuous outcome and continuous mediator
#' mydata1 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=0)
#' 
#' # generate continuous outcome and binary mediator
#' sigma_tau = 0.605    #using this value can produce ICC of mediator model is about 0.1 since  
#' # ICC of binary variable is calculated as sigma_tau^2/(sigma_tau^2+pi^2/3)
#' mydata2 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=1)
#' 
#' # generate binary outcome and continuous mediator
#' sigma_a = 0.605
#' mydata3 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=1,binary.m=0)
#' 
#' # generate binary outcome and binary mediator
#' sigma_tau = sigma_a = 0.605
#' mydata4 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=1,binary.m=1)
#' 
#' # generate continuous outcome and continuous mediator with different clusters and period
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' sigma_tau = sigma_a = 0.334
#' sigma_em = sigma_ey = 1
#' mydata5 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
#' sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=0)


gen_data_hhm = function(I,J,n,beta,gamma,theta,beta_M,eta,sigma_a,sigma_ey,
                        sigma_tau,sigma_em,binary.o=0,binary.m=0){
  # Generate data frame
  data = data.frame(
    "cluster" = integer(), # cluster
    "period" = integer(), # time; 1=baseline, J=endline
    "id" = integer(), # individual
    "E" = integer(), # time since intervention
    "A" = integer(), # treatment state indicator
    "c" = integer(), # the start time of the treatment
    "alpha" = numeric(), # random effects of observations within a cluster
    "tau" = numeric(), # random effects of mediation within a cluster
    "M" = numeric(), # mediator
    "Y" = numeric()  # outcome
  )
  # Generate crossover times (assumes a "balanced and complete" design)
  n_clust_per_time = I/(J-1)
  if (n_clust_per_time %% 1 != 0) {
    stop("n_clusters must be divisible by period-1")
  }
  crossover_times = rep(2:J, each=n_clust_per_time)
  
  # Loop through clusters, time, and individuals
  for(i in 1:I){ #I: number of cluster
    alpha_i = rnorm(1,0,sigma_a)
    tau_i = rnorm(1,0,sigma_tau)
    c_i = crossover_times[i]-1 # the start time of the treatment
    for(j in 1:J){
      A_ij = ifelse(j < crossover_times[i], 0, 1)
      e = ifelse(j < crossover_times[i], 0, (j-crossover_times[i]) + 1) # time since treatment
      e_ijk = rnorm(n,0,sigma_em)
      epsilon_ijk = rnorm(n,0,sigma_ey)
      if(binary.o == 0 & binary.m == 0){ #type 1: continuous outcome and continuous mediator
        m_ijk = gamma[j] + eta*A_ij + tau_i + e_ijk
        y_ijk = beta[j] + theta*A_ij + beta_M*m_ijk + alpha_i + epsilon_ijk
      }else if(binary.o == 0 & binary.m == 1){ #type 2: continuous outcome and binary mediator
        mu_ij = gamma[j] + eta*A_ij + tau_i
        mu_ij = exp(mu_ij)/(1+exp(mu_ij))  #to make mu_ij between 0 and 1
        m_ijk = rbinom(n=n, size=1, prob=mu_ij)
        y_ijk = beta[j] + theta*A_ij + beta_M*m_ijk + alpha_i + epsilon_ijk
      }else if(binary.o == 1 & binary.m == 0){ #type 3: binary outcome and continuous mediator
        m_ijk = gamma[j] + eta*A_ij + tau_i + e_ijk
        mu_ijy = beta[j] + theta*A_ij + beta_M*m_ijk + alpha_i
        mu_ijy = exp(mu_ijy)/(1+exp(mu_ijy))
        y_ijk = rbinom(n=n, size=1, prob=mu_ijy)
      }else{  #type 4: binary outcome and binary mediator
        mu_ij = gamma[j] + eta*A_ij + tau_i
        mu_ij = exp(mu_ij)/(1+exp(mu_ij))  #to make mu_ij between 0 and 1
        m_ijk = rbinom(n=n, size=1, prob=mu_ij)
        mu_ijy =  beta[j] + theta*A_ij + beta_M*m_ijk + alpha_i
        mu_ijy = exp(mu_ijy)/(1+exp(mu_ijy))
        y_ijk = rbinom(n=n, size=1, prob=mu_ijy)
      } 
      data = rbind(data, data.frame(cbind(
        cluster = rep(i,n), period = rep(j,n), id = c(1:n), E = rep(e,n),
        A = rep(A_ij,n), c = rep(c_i,n), alpha = rep(alpha_i,n),
        tau = rep(tau_i,n), M = m_ijk, Y=y_ijk))
      )
    }
  }
  dataset = as.data.frame(data)
  return(dataset)
}
