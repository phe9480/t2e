% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexPower.R
\name{flexPower}
\alias{flexPower}
\title{Power Calculation Based on Flexible Distributions in Group Sequential Design}
\usage{
flexPower(
  T = c(24, 36),
  n = 450,
  r = 1/2,
  sf = gsDesign::sfLDOF,
  h0 = function(t) {     log(2)/12 },
  S0 = function(t) {     exp(-log(2)/12 * t) },
  h1 = function(t) {     log(2)/12 * 0.7 },
  S1 = function(t) {     exp(-log(2)/12 * 0.7 * t) },
  Fentry = function(t) {     (t/24)^2 * as.numeric(t < 24) + as.numeric(t >= 24) },
  drop.rate = 0.03/12,
  bd.p = NULL,
  alpha = 0.025,
  variance = "H1"
)
}
\arguments{
\item{T}{A vector of analysis times for IAs and FA, calculated from first subject in}

\item{n}{Total sample size for two arms}

\item{r}{Proportion of experimental arm subjects among all subjects. For 2:1 randomization, r = 2/3}

\item{sf}{Alpha spending function. Default sf = gsDesign::sfLDOF, i.e.,
LanDeMets implementation of O'Brien Fleming spending function.}

\item{h0}{Hazard function for control arm}

\item{S0}{Survival function for control arm}

\item{h1}{Hazard function for experimental arm}

\item{S1}{Survival function for experimental arm}

\item{Fentry}{Entry cumulative distribution function. For enrolment period
of E months, Fentry is an increasing function from 0 to E, and equals 1 after 1.}

\item{drop.rate}{Drop-off rate in month. For example, if drop-off rate is
3\% for every year of follow-up, then drop.rate = 0.03/12.}

\item{bd.p}{A vector of rejection bound in one-sided p value. When bd.p is
specified, sf and alpha are ignored.}

\item{alpha}{Overall alpha for the test in the group sequential design.
Default, alpha = 0.025. Alpha must be one-sided. When bd.p is specified, alpha
will be re-calculated. Required when bd.p is not specified.}

\item{variance}{Option for variance estimate. "H1" or "H0". Default H1,
which is usually more conservative than H0.}
}
\value{
\itemize{
\item events     Target events at each analysis
\item bd         Rejection ounds at each analysis
\item marg.power Marginal power at each analysis
\item cum.power  Cumulative power by each analysis
\item overall.power Overall power
\item CV         Critical value (minimum detectable difference)
\item corr       Correlation matrix
}
}
\description{
This function calculates the power for any customized survival function in
each arm and rejection boundary at each analysis.
}
\examples{

#(1) O'Brien Fleming alpha spending, variance under H0
flexPower(T = c(36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
                   h0=function(t){log(2)/12},
                   S0=function(t){exp(-log(2)/12*t)},
                   h1=function(t){log(2)/12*0.70}, 
                   S1=function(t){exp(-log(2)/12 * 0.7 * t)},
                   Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
                   drop.rate = 0.03/12,
                   bd.p=NULL,
                   alpha=0.025,
                   variance="H0")
flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
                   h0=function(t){log(2)/12},
                   S0=function(t){exp(-log(2)/12*t)},
                   h1=function(t){log(2)/12*0.70}, 
                   S1=function(t){exp(-log(2)/12 * 0.7 * t)},
                   Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
                   drop.rate = 0.03/12,
                   bd.p=NULL,
                   alpha=0.025,
                   variance="H0")
simplePower(events=c(132.4076, 267.0104), events0=NULL, events1=NULL, 
                        hr = 0.7, r = 0.5, 
                        bd.p=NULL, sf=gsDesign::sfLDOF, 
                        alpha=0.025, variance="H0")
#(2) O'Brien Fleming alpha spending, variance under H1
flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
                   h0=function(t){log(2)/12},
                   S0=function(t){exp(-log(2)/12*t)},
                   h1=function(t){log(2)/12*0.70}, 
                   S1=function(t){exp(-log(2)/12 * 0.7 * t)},
                   Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
                   drop.rate = 0.03/12,
                   bd.p=NULL,
                   alpha=0.025,
                   variance="H1")
simplePower(events=c(132.4076, 267.0104), events0=c(75.08895,147.45471), events1=c(57.31865,119.55565), 
                        hr = 0.7, r = 0.5, 
                        bd.p=NULL, sf=gsDesign::sfLDOF, 
                        alpha=0.025, variance="H1")

#(3) Customized rejection bounds, variance under H1
flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
                   h0=function(t){log(2)/12},
                   S0=function(t){exp(-log(2)/12*t)},
                   h1=function(t){log(2)/12*0.70}, 
                   S1=function(t){exp(-log(2)/12 * 0.7 * t)},
                   Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
                   drop.rate = 0.03/12,
                   bd.p=c(0.002, 0.023),
                   alpha=0.025,
                   variance="H1")
simplePower(events=c(132.4076, 267.0104), events0=c(75.08895,147.45471), events1=c(57.31865,119.55565), 
                        hr = 0.7, r = 0.5, 
                        bd.p=c(0.002, 0.023), sf=gsDesign::sfLDOF, 
                        alpha=0.025, variance="H1")
                        
}
