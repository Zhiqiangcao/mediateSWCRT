% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediate_swcrt.R
\name{mediate_swcrt}
\alias{mediate_swcrt}
\title{Mediation analysis in a stepped wedge cluster randomized trials based on Hussey and Hughes model}
\usage{
mediate_swcrt(
  data,
  method = "STA",
  outcome = "Y",
  mediator = "M",
  treatment = "A",
  cluster = "cluster",
  period = "period",
  exposure = "E",
  covariate.outcome = NULL,
  covariate.mediator = NULL,
  a0 = 0,
  a1 = 1,
  binary.o = 0,
  binary.m = 0,
  time.dependent = FALSE
)
}
\arguments{
\item{data}{(Required) The name of the dataset}

\item{method}{Two approximate methods are provided, one is STA (second-order Taylor approximate), the other is GHQ (Gauss-Hermite Quadrature)
when both the outcome and mediator are continuous, the method choice is irrelevant; when the outcome is binary and mediator is continuous, we default to the STA method}

\item{outcome}{(Required) Name of the outcome variable, which should be either a continuous or binary datatype}

\item{mediator}{(Required) Name of the mediator variable, which should be either a continuous or binary datatype}

\item{treatment}{(Required) Name of the treatment variable, which should be either a continuous or binary datatype. When there is a implementation period, then corresponding treatment status
should be set to -1, which is mainly for convenience}

\item{cluster}{(Required) Name of the cluster variable, which should be a factor}

\item{period}{(Required) Name of the period variable, which should be a factor}

\item{exposure}{Name of the exposure time variable, which should be a factor. For time-dependent model, this argument is required. For constant treatment effect model,
if dataset has no exposure variable, then set it to NULL}

\item{covariate.outcome}{A vector of names showing the confounding variables used in the outcome model. The default value
is NULL, which represents no confounding variables. We only accepted continuous and binary confounding variables,
if one confounding variable is categorical, please set it to a series of binary variables in advance}

\item{covariate.mediator}{A vector of names showing the confounding variables used in the mediator model. The default value
is NULL, which represents no confounding variables. We only accepted continuous and binary confounding variables, if one
confounding variable is categorical, please set it to a series of binary variables in advance}

\item{a0}{The reference treatment level in defining TE and NIE. The default value is 0}

\item{a1}{The treatment level of interest in defining TE and NIE. The default value is 1}

\item{binary.o}{(Required) If the outcome is binary, set to 1. If the outcome is continuous, set to 0}

\item{binary.m}{(Required) If the mediator is binary, set to 1. If the mediator is continuous, set to 0}

\item{time.dependent}{(Required) If the treatment effect is time-dependent, set to TRUE. If the treatment effect is constant, set to FALSE}
}
\value{
A list containing point and interval estimates of NIE, NDE, TE and MP under constant treatment effect structure; or a list
containing point and interval estimates of NIE, NDE, TE and MP at each exposure time and overall summary mediation measures,
as well as Chi-square test result of TE under exposure-time treatment effect structure
}
\description{
Function for calculating the point and interval estimates for the NIE, NDE, TE and MP in a stepped
wedge cluster randomized trials with constant treatment effect or exposure-time dependent treatment effect
}
\examples{
library(lme4)
I = 15 
J = 4
n = 20
beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
sigma_tau = sigma_a = 0.334
sigma_em = sigma_ey = 1
set.seed(123456)
mydata1 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=0.625,eta=0.4,sigma_a,
sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=0)
covdata = data.frame("X1" = numeric(), "X2" = numeric())
set.seed(100)
for(i in 1:I){
  for(j in 1:J){
     x1_ijk = rnorm(n) 
     x2_ijk = rnorm(n,sd=0.1)
     covdata = rbind(covdata, data.frame(cbind(X1 = x1_ijk, X2 = x2_ijk)))
  }
}
mydata1 = data.frame(mydata1,covdata)
# example 1: constant treatment effect, mediation analysis of continuous outcome 
# and continuous mediator without covariates, with different covariates in outcome 
# and mediator models, as well as with the same covariates in both models
res1 = mediate_swcrt(data=mydata1)
print(res1)
res1f = mediate_swcrt(data=mydata1,covariate.outcome = c("X1"),
covariate.mediator = c("X2"))
print(res1f)
res1ff = mediate_swcrt(data=mydata1,covariate.outcome = c("X1","X2"),
covariate.mediator = c("X1","X2"))
print(res1ff)

# example 2: constant treatment effect, mediation analysis of binary outcome 
# and continuous mediator without covariates
I = 15 
J = 4
n = 20
beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
sigma_a = 0.605
sigma_tau = 0.334
sigma_em = sigma_ey = 1
set.seed(123456)
mydata2 = gen_data_hhm(I,J,n,beta,gamma,theta=0.75,beta_M=1,eta=0.4,sigma_a,
sigma_ey,sigma_tau,sigma_em,binary.o=1,binary.m=0)
res2 = mediate_swcrt(data=mydata2,binary.o = 1)
print(res2)

# example 3: exposure-time dependent treatment effect, mediation analysis of
# continuous outcome and continuous mediator without covariates
I = 15 
J = 4
n = 20
beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
sigma_tau = sigma_a = 0.334
sigma_em = sigma_ey = 1
eta_e = c(0.5,0.8,1.1) 
theta_e = c(0.5,1.0,1.5)
set.seed(123456)
mydata3 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=0.8,eta_e,sigma_a,
sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=0)
res3 = mediate_swcrt(data=mydata3,time.dependent = TRUE)
print(res3)

# example 4: exposure-time dependent treatment effect, mediation analysis of 
# continuous outcome and binary mediator with the same covariates in both models
I = 15 
J = 4
n = 20
beta = cumsum(c(0,0.1,0.1/2,0.1/(2^2)))
gamma = cumsum(c(0,0.3,0.3/2,0.3/(2^2)))
eta_e = c(0.5,1.3,2.1) 
theta_e = c(0.6,1,1.4)
sigma_a = 0.334
sigma_tau = 0.605   
sigma_em = sigma_ey = 1
set.seed(123456)
mydata4 = gen_data_etm(I,J,n,beta,gamma,theta_e,beta_M=1.2,eta_e,sigma_a,
sigma_ey,sigma_tau,sigma_em,binary.o=0,binary.m=1)
set.seed(100)
for(i in 1:I){
  for(j in 1:J){
     x1_ijk = rnorm(n) 
     x2_ijk = rnorm(n,sd=0.5)
     covdata = rbind(covdata, data.frame(cbind(X1 = x1_ijk, X2 = x2_ijk)))
  }
}
mydata4 = data.frame(mydata4,covdata)
# using default STA method
res4 = mediate_swcrt(data=mydata4,covariate.outcome = c("X1","X2"), 
covariate.mediator = c("X1","X2"),binary.m = 1, time.dependent = TRUE)
print(res4)
# using GHQ method
res4f = mediate_swcrt(data=mydata4,method = "GHQ",covariate.outcome = c("X1","X2"), 
covariate.mediator = c("X1","X2"),binary.m = 1, time.dependent = TRUE)
print(res4f)


}
\references{
Cao Z. and Fan Li. Assessing mediation in cross-sectional stepped wedge cluster randomized trials. under review. 2024;0(0):1-24.
}
\author{
Zhiqiang Cao \email{zcaoae@connect.ust.hk} and Fan Li \email{fan.f.li@yale.edu}
}
