#' Generate data in a stepped wedge design based on nested exchangeable correlation model with an exposure-time dependent treatment effect
#' 
#' This function can generate longitudinal data in a stepped wedge design nested exchangeable correlation model in the presence of
#' an exposure-time dependent treatment effect
#' 
#' @param I The number of cluster
#' @param J The number of time periods
#' @param n The number of individuals per each cluster
#' @param beta The underling time trend effect in a outcome model
#' @param gamma The underling time trend effect in a mediator model
#' @param theta_e The exposure effect (a J-1 vector) on the outcome conditional on mediator and other possible confounders
#' @param beta_M The mediator effect on the outcome conditional on other possible confounders
#' @param eta_e The exposure effect (a J-1 vector) on the mediator conditional on other possible confounders
#' @param sigma_a The standard error of random effect in outcome model (assuming it follows normal distribution with mean zero)
#' @param sigma_phi The standard error of random effect of cluster-by-time interaction in outcome model (assuming it follows normal distribution with mean zero)
#' @param sigma_ey The standard error of error term in outcome model (assuming it follows normal distribution with mean zero)
#' @param sigma_tau The standard error of random effect in mediator model (assuming it follows normal distribution with mean zero)
#' @param sigma_psi The standard error of random effect of cluster-by-time interaction in mediator model (assuming it follows normal distribution with mean zero)
#' @param sigma_em The standard error of error term in mediator model (assuming it follows normal distribution with mean zero)
#' @param binary.o (Required) If the outcome is binary, set to 1. If the outcome is continuous, set to 0
#' @param binary.m (Required) If the mediator is binary, set to 1. If the mediator is continuous, set to 0
#' 
#' @return A data frame including cluster, period, id, E (time since intervention, i.e., exposure time), A (treatment status indicator),
#'         c (the start time of the treatment), alpha (random effects of outcome within a cluster), phi (random effects of cluster-by-time interaction in outcome model),
#'         tau (random effects of mediator within a cluster), psi (random effects of cluster-by-time interaction in mediator model), M (mediator), Y (outcome)
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk} and Fan Li \email{fan.f.li@@yale.edu} 
#' 
#' @references Li F., Hughes J.P., Hemming K., Taljaard M., Melnick E.R. and Heagerty P.J. Mixed-effects models for the design
#' and analysis of stepped wedge cluster randomized trials: An overview. Statistical Methods in Medical Research. 2021;30(2):612-639.
#' 
#' Kenny A., Voldal E.C., Xia F., Heagerty P.J. and Hughes J.P. Analysis of stepped wedge cluster randomized trials in the
#' presence of a time-varying treatment effect. Statistics in Medicine. 2022;41(22):4311-4339.
#' 
#' @importFrom stats rnorm 
#' @importFrom stats rbinom
#' @examples
#' I = 9  
#' J = 4
#' n = 20
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
#' sigma_tau = sigma_a = 0.334
#' sigma_phi = sigma_psi = 0.6
#' sigma_em = sigma_ey = 0.8
#' eta_e = c(0.32,0.4,0.48) 
#' theta_e = c(0.6,0.75,0.9)
#' # generate continuous outcome and continuous mediator
#' mydata1 = gen_data_nem_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=0)
#' 
#' # generate continuous outcome and binary mediator
#' mydata2 = gen_data_nem_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=1)
#' 
#' # generate binary outcome and continuous mediator
#' mydata3 = gen_data_nem_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=1,binary.m=0)
#' 
#' # generate binary outcome and binary mediator
#' mydata4 = gen_data_nem_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=1,binary.m=1)
#' 
#' # generate continuous outcome and continuous mediator with different clusters and period
#' I = 12  
#' J = 5
#' beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
#' gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2),0.3/(2^3)))
#' eta_e = c(0.32,0.4,0.48,0.56) 
#' theta_e = c(0.6,0.75,0.9,1.05)
#' mydata5 = gen_data_nem_etm(I,J,n,beta,gamma,theta_e,beta_M=0.625,eta_e,sigma_a,
#' sigma_phi,sigma_ey,sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=0)
#' 

gen_data_nem_etm = function(I,J,n,beta,gamma,theta_e,beta_M,eta_e,sigma_a,sigma_phi,sigma_ey,
                        sigma_tau,sigma_psi,sigma_em,binary.o=0,binary.m=0){
  # Generate data frame
  data <- data.frame(
    "cluster" = integer(), # cluster
    "period" = integer(), # time; 1=baseline, J=endline
    "id" = integer(), # individual
    "E" = integer(), # time since intervention
    "A" = integer(), # treatment state indicator
    "c" = integer(), # the start time of the treatment
    "alpha" = numeric(), # random effects of observations within a cluster
    "phi" = numeric(), # random effects of cluster-by-time interaction in outcome model
    "tau" = numeric(), # random effects of mediation within a cluster
    "psi" = numeric(), # random effects of cluster-by-time interaction in mediator model
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
      thetae = ifelse(e > 0, theta_e[e], 0)
      etae = ifelse(e > 0, eta_e[e], 0) 
      phi_ij = rnorm(1,0,sigma_phi)
      psi_ij = rnorm(1,0,sigma_psi)
      e_ijk = rnorm(n,0,sigma_em)
      epsilon_ijk = rnorm(n,0,sigma_ey)
      if(binary.o == 0 & binary.m == 0){ #type 1: continuous outcome and continuous mediator
        m_ijk = gamma[j] + etae*A_ij + tau_i + psi_ij + e_ijk
        y_ijk = beta[j] + thetae*A_ij + beta_M*m_ijk + alpha_i + phi_ij + epsilon_ijk
      }else if(binary.o == 0 & binary.m == 1){ #type 2: continuous outcome and binary mediator
        mu_ij = gamma[j] + etae*A_ij + tau_i + psi_ij
        mu_ij = exp(mu_ij)/(1+exp(mu_ij))  #to make mu_ij between 0 and 1
        m_ijk = rbinom(n=n, size=1, prob=mu_ij)
        y_ijk = beta[j] + thetae*A_ij + beta_M*m_ijk + alpha_i + phi_ij + epsilon_ijk
      }else if(binary.o == 1 & binary.m == 0){ #type 3: binary outcome and continuous mediator
        m_ijk = gamma[j] + etae*A_ij + tau_i + psi_ij + e_ijk
        mu_ijy = beta[j] + thetae*A_ij + beta_M*m_ijk + alpha_i + phi_ij
        mu_ijy = exp(mu_ijy)/(1+exp(mu_ijy))
        y_ijk = rbinom(n=n, size=1, prob=mu_ijy)
      }else{  #type 4: binary outcome and binary mediator
        mu_ij = gamma[j] + etae*A_ij + tau_i + psi_ij
        mu_ij = exp(mu_ij)/(1+exp(mu_ij))  #to make mu_ij between 0 and 1
        m_ijk = rbinom(n=n, size=1, prob=mu_ij)
        mu_ijy =  beta[j] + thetae*A_ij + beta_M*m_ijk + alpha_i + phi_ij
        mu_ijy = exp(mu_ijy)/(1+exp(mu_ijy))
        y_ijk = rbinom(n=n, size=1, prob=mu_ijy)
      }
      data = rbind(data, data.frame(cbind(
        cluster = rep(i,n), period = rep(j,n), id = c(1:n), E = rep(e,n),
        A = rep(A_ij,n), c = rep(c_i,n), alphi = rep(alpha_i,n), phi = rep(phi_ij,n),
        tau = rep(tau_i,n), psi = rep(psi_ij,n), M = m_ijk, Y=y_ijk))
      )
    }
  }
  dataset <- as.data.frame(data)
  return(dataset)
}

