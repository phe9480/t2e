#' Power Calculation Based on Provided Number of Events And Rejection Boundary in Group Sequential Design
#' 
#' This function calculates the power for the given number of events and rejection boundary
#' at each analysis for group sequential design
#' 
#' @param events  A vector of target events for IAs and FA
#' @param events0 A vector of target events for IAs and FA for control arm. 
#' Required when variance option is H1.
#' @param events1 A vector of target events for IAs and FA for experimental arm
#' Required when variance option is H1.
#' @param hr      Hazard ratio (smaller better)
#' @param r       Proportion of experimental arm subjects among all subjects
#' @param bd.p  A vector of rejection bound in one-sided p value. When bd.p is 
#' specified, sf and alpha are ignored.
#' @param sf Alpha spending function. Default sf = gsDesign::sfLDOF, i.e., 
#' LanDeMets implementation of O'Brien Fleming spending function.
#' @param sfpar Parameter in alpha spending function when applicable. For example,
#'              HSD(-3) is specified as sf = gsDesign::sfHSD, sfupar = -3. 
#' LanDeMets implementation of O'Brien Fleming spending function.
#' @param alpha Overall alpha for the test in the group sequential design. 
#' Default, alpha = 0.025. Alpha must be one-sided. When bd.p is specified, alpha
#' will be re-calculated. Required when bd.p is not specified.
#' @param variance Option for variance estimate. "H1" or "H0". Default H1, which is usually more conservative than H0.
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
#' simplePower(events=c(126, 210), events0=NULL, events1=NULL, 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=NULL, sf=gsDesign::sfLDOF, 
#'                         alpha=0.025, variance="H0")
#' #(2) O'Brien Fleming alpha spending, variance under H1
#' simplePower(events=c(126, 210), events0=c(66,115), events1=c(60,95), 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=NULL, sf=gsDesign::sfLDOF, 
#'                         alpha=0.025, variance="H1")
#' 
#' #(3) Customized rejection bounds, variance under H1
#' simplePower(events=c(126, 210), events0=c(66,115), events1=c(60,95), 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=c(0.004,0.024), sf=NULL, 
#'                         alpha=NULL, variance="H1")
#' simplePower(events=c(208, 287), events0=c(66,115), events1=c(60,95), 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=c(0.016933/2,0.044898/2), sf=NULL, 
#'                         alpha=NULL, variance="H0")
#' simplePower(events=c(231, 287), events0=c(66,115), events1=c(60,95), 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=c(0.016933/2,0.044898/2), sf=NULL, 
#'                         alpha=NULL, variance="H0")
#'                         
#' #(4) HSD(-3) alpha spending, variance under H0
#' simplePower(events=c(126, 210), events0=NULL, events1=NULL, 
#'                         hr = 0.6, r = 0.5, 
#'                         bd.p=NULL, sf=gsDesign::sfHSD, sfpar=-3,
#'                         alpha=0.025, variance="H0")
#' @export
#'
simplePower = function(events=c(126, 210), events0=NULL, events1=NULL, 
                        hr = 0.6, r = 0.5, 
                        bd.p=NULL, sf=gsDesign::sfLDOF, sfpar=-3,
                        alpha=0.025, variance="H0"){
  
  #Number of analyses
  if (variance=="H1"){ 
    events = events1 + events0
    re = events1/events
  }
  K = length(events)
  
  #(1) Mean of z statistics
  #######By default, use variance under H1 for more conservative estimate of power
  if (variance == "H0"){
    mu = -log(hr) * sqrt(r*(1-r)*events)
  } else {
    #Proportion of events
    mu = -log(hr) * sqrt(re*(1-re)*events)
  }

  #(2) Correlation matrix
  corr = matrix(1, nrow=K, ncol=K)
  if(K > 1){
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        corr[i,j] = corr[j,i] = sqrt(events[i]/events[j])
      }
    }  
  }
  
  #(3) Bound
  timing = events/events[K]
  
  if (K > 1){
    if(is.null(bd.p)){
      d = gsDesign::gsDesign(k=K,alpha=alpha,timing=timing,sfu=sf, sfupar=sfpar, test.type=1)$upper
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

  #Population A and B
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
    cv = exp(-z/sqrt(r*(1-r)*events))
  } else {   
    cv = exp(-z/sqrt(re*(1-re)*events))
  }
  side = 1
  
  power=data.frame(cbind(timing,marg.power,incr.power,cum.power,overall.power))
  bd = data.frame(cbind(bd.p, z, side))
  o = list()
  o$overall.alpha = data.frame(alpha)
  o$events = data.frame(events)
  o$power = power
  o$bd = bd
  o$CV = data.frame(cv)
  o$corr = corr

  return(o)
}  
