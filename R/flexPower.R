#' Power Calculation Based on Flexible Distributions in Group Sequential Design
#' 
#' This function calculates the power for any customized survival function in 
#' each arm and rejection boundary at each analysis. 
#' 
#' @param T  A vector of analysis times for IAs and FA, calculated from first subject in
#' @param n  Total sample size for two arms 
#' @param r  Proportion of experimental arm subjects among all subjects. For 2:1 randomization, r = 2/3
#' @param sf Alpha spending function. Default sf = gsDesign::sfLDOF, i.e., 
#' LanDeMets implementation of O'Brien Fleming spending function.
#' @param h0 Hazard function for control arm
#' @param S0 Survival function for control arm
#' @param h1 Hazard function for experimental arm
#' @param S1 Survival function for experimental arm
#' @param Fentry Entry cumulative distribution function. For enrolment period 
#' of E months, Fentry is an increasing function from 0 to E, and equals 1 after 1.
#' @param drop.rate Drop-off rate in month. For example, if drop-off rate is 
#' 3% for every year of follow-up, then drop.rate = 0.03/12.
#' @param bd.p  A vector of rejection bound in one-sided p value. When bd.p is 
#' specified, sf and alpha are ignored.
#' @param alpha Overall alpha for the test in the group sequential design. 
#' Default, alpha = 0.025. Alpha must be one-sided. When bd.p is specified, alpha
#' will be re-calculated. Required when bd.p is not specified.
#' @param variance Option for variance estimate. "H1" or "H0". Default H1, 
#' which is usually more conservative than H0.
#'  
#' @return 
#'  \itemize{
#'  \item events     Target events at each analysis
#'  \item bd         Rejection ounds at each analysis
#'  \item marg.power Marginal power at each analysis
#'  \item cum.power  Cumulative power by each analysis
#'  \item overall.power Overall power
#'  \item CV         Critical value (minimum detectable difference)    
#'  \item corr       Correlation matrix
#'  }
#' @examples
#' 
#' #(1) O'Brien Fleming alpha spending, variance under H0
#' flexPower(T = c(36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
#'                    h0=function(t){log(2)/12},
#'                    S0=function(t){exp(-log(2)/12*t)},
#'                    h1=function(t){log(2)/12*0.70}, 
#'                    S1=function(t){exp(-log(2)/12 * 0.7 * t)},
#'                    Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
#'                    drop.rate = 0.03/12,
#'                    bd.p=NULL,
#'                    alpha=0.025,
#'                    variance="H0")
#' flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
#'                    h0=function(t){log(2)/12},
#'                    S0=function(t){exp(-log(2)/12*t)},
#'                    h1=function(t){log(2)/12*0.70}, 
#'                    S1=function(t){exp(-log(2)/12 * 0.7 * t)},
#'                    Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
#'                    drop.rate = 0.03/12,
#'                    bd.p=NULL,
#'                    alpha=0.025,
#'                    variance="H0")
#' simplePower(events=c(132.4076, 267.0104), events0=NULL, events1=NULL, 
#'                         hr = 0.7, r = 0.5, 
#'                         bd.p=NULL, sf=gsDesign::sfLDOF, 
#'                         alpha=0.025, variance="H0")
#' #(2) O'Brien Fleming alpha spending, variance under H1
#' flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
#'                    h0=function(t){log(2)/12},
#'                    S0=function(t){exp(-log(2)/12*t)},
#'                    h1=function(t){log(2)/12*0.70}, 
#'                    S1=function(t){exp(-log(2)/12 * 0.7 * t)},
#'                    Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
#'                    drop.rate = 0.03/12,
#'                    bd.p=NULL,
#'                    alpha=0.025,
#'                    variance="H1")
#' simplePower(events=c(132.4076, 267.0104), events0=c(75.08895,147.45471), events1=c(57.31865,119.55565), 
#'                         hr = 0.7, r = 0.5, 
#'                         bd.p=NULL, sf=gsDesign::sfLDOF, 
#'                         alpha=0.025, variance="H1")
#' 
#' #(3) Customized rejection bounds, variance under H1
#' flexPower(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
#'                    h0=function(t){log(2)/12},
#'                    S0=function(t){exp(-log(2)/12*t)},
#'                    h1=function(t){log(2)/12*0.70}, 
#'                    S1=function(t){exp(-log(2)/12 * 0.7 * t)},
#'                    Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
#'                    drop.rate = 0.03/12,
#'                    bd.p=c(0.002, 0.023),
#'                    alpha=0.025,
#'                    variance="H1")
#' simplePower(events=c(132.4076, 267.0104), events0=c(75.08895,147.45471), events1=c(57.31865,119.55565), 
#'                         hr = 0.7, r = 0.5, 
#'                         bd.p=c(0.002, 0.023), sf=gsDesign::sfLDOF, 
#'                         alpha=0.025, variance="H1")
#'                         
#' @export
#' 
flexPower = function(T = c(24, 36), n = 450, r= 1/2, sf=gsDesign::sfLDOF, 
                      h0=function(t){log(2)/12},
                      S0=function(t){exp(-log(2)/12*t)},
                      h1=function(t){log(2)/12*0.70}, 
                      S1=function(t){exp(-log(2)/12 * 0.7 * t)},
                      Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)},
                      drop.rate = 0.03/12,
                      bd.p=NULL,
                      alpha=0.025,
                      variance="H1"){
  
  #Number of analyses
  K = length(T)
  

  #Density functions  
  f0 = function(t){h0(t)*S0(t)}
  f1 = function(t){h1(t)*S1(t)}
  hr = function(t){h1(t)/h0(t)}
  f = function(t){(1-r)*f0(t) + r*f1(t)}
  
  #Number of events
  e0 = e1 = mu = rep(NA, K)
  for (k in 1:K){
    e0[k] = fevents(DCO = T[k], f=f0, Fentry=Fentry, drop.rate = drop.rate, n=n*(1-r))
    e1[k] = fevents(DCO = T[k], f=f1, Fentry=Fentry, drop.rate = drop.rate, n=n*r)    
  }
  e = e0 + e1
  
  #(1) Mean of z statistics
  G = function(t){1-exp(-drop.rate*t)}

  for(k in 1:K){
    f.mu.top = function(t){
      -log(hr(t))*Fentry(T[k] - t) * (1-G(t)) * f(t)
    }
    f.mu.bot = function(t){
      Fentry(T[k] - t) * (1-G(t)) * f(t)
    }
    mu.top = integrate(f=f.mu.top,lower=0,upper=T[k])$value
    mu.bot = integrate(f=f.mu.bot,lower=0,upper=T[k])$value
    #######By default, use variance under H1 for more conservative estimate of power
    if (variance == "H0"){
      mu[k] = sqrt(n*r*(1-r))*mu.top/sqrt(mu.bot)
    } else {
      #Proportion of events
      re = e1[k]/e[k]
      mu[k] = sqrt(n*re*(1-re))*mu.top/sqrt(mu.bot)
    }
  }  
  
  #(2) Correlation matrix
  corr = matrix(1, nrow=K, ncol=K)
  if(K > 1){
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        corr[i,j] = corr[j,i] = sqrt(e[i]/e[j])
      }
    }  
  }
  
  #(3) Bound
  timing = e/e[K]
  
  if (K > 1){
    if(is.null(bd.p)){
      d = gsDesign::gsDesign(k=K,alpha=alpha,timing=timing,sfu=sf, test.type=1)$upper
      z <- d$bound   #z boundary
      bd.p = 1-pnorm(z) #p-value boundary
    } else{
      z = qnorm(1-bd.p)
      #the overall alpha needs update
      alpha = 1 - mvtnorm::pmvnorm(lower = rep(-Inf, K), 
                                 upper = z, mean=rep(0,K), 
                                 corr = corr, abseps = 1e-8, maxpts=100000)[1]
    }
  } else {
    if(is.null(bd.p)){
      bd.p = alpha
      z = qnorm(1-bd.p)
    } else{
      z = qnorm(1-bd.p)
      #the overall alpha needs update
      alpha = 1 - pnorm(z)
    }
  }
  
  #(4)Marginal Power
  marg.power = rep(NA, K)
  
  for (k in 1:K){
    marg.power[k] = 1-pnorm(qnorm(1-bd.p[k]), mean=mu[k])
  }
  
  #(4c)Incremental Power
  incr.power = rep(NA, K)
  cum.power =  rep(NA, K)
  cum.power[1] = incr.power[1] = marg.power[1]
  
  
  if(K > 1){
    for (k in 2:K){
      incr.power[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), z[k]), 
                                       upper = c(z[1:(k-1)], Inf), mean=mu[1:k], 
                                       corr = corr[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      cum.power[k] = cum.power[k-1] + incr.power[k]
    }
    overall.power = cum.power[K]
  } else{
    cum.power = incr.power = marg.power = overall.power = 1-pnorm(z, mean=mu)
  }
  
  #(5) Critical Values in HR
  #######By default, use variance under H1 for more conservative estimate of power
  if (variance == "H0"){
    cv = exp(-z/sqrt(r*(1-r)*e))
  } else {   
  re = e1/e
  cv = exp(-z/sqrt(re*(1-re)*e))
  }
  side = 1
  
  power=data.frame(cbind(timing,marg.power,incr.power,cum.power,overall.power))
  bd = data.frame(cbind(bd.p, z, side))
  o = list()
  o$overall.alpha = data.frame(alpha)
  o$events = data.frame(cbind(e0, e1, e))
  o$power = power
  o$bd = bd
  o$mu = mu
  o$CV = data.frame(cv)
  o$corr = corr
  
  return(o)
}  
