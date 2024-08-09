#' Hazard ratio of overall population for heterogeneous population
#' 
#' This function calculates the hazard ratio function of overall population when
#' there is heterogeneous effect in subgroups. In general, given the hazard ratio
#' for each mutually exclusive subgroup, the overall hazard ratio in the overall
#' population combining all subgroups is a function of survival time. This function
#' calculates the average hazard ratio, which is defined as the expected value of
#' the overall hazard ratio function, i.e, average Hazard Ratio = E(HR(t)), which 
#' can be evaluated based on S0(t), however in common situations, it is generally
#' robust to the specification of S0(t).
#' 
#' @param S0 Control arm survival function. Optional. The average HR is generally
#' robust with respect to S0. If S0 is not provided, an exponential distribution 
#' with median 12 months will be used by default.
#' @param h0 Control arm hazard rate functio. Optional. For exponential distribution,
#' h0 is a constant. If h0 is not provided, constant hazard of log(2)/12 is used.
#' If S0 and h0 are provided, they must be consistent to each other, ie, S0*h0=f0.
#' @param HR A vector of hazard ratios for subgroups. 
#' @param p  A vector of proportions for subgroups. sum(p) must be 1. HR and p
#' must be in the same order of subgroups. 
#'
#' @return Average hazard ratio in the overall population combining all mutually
#'  exclusive subgroups. 
#'  
#' @examples
#' #Example 1. PD-L1+/- subgroups have HRs 0.65 and 0.75 with prevalence of 35% and 65% 
#' #respectively. The average hazard ratio of the overall population (PD-L1+/-) is estimated as 
#' HR.m12 = overallAHR(S0 = function(t){exp(-log(2)/12*t)}, h0=function(t){log(2)/12}, 
#' HR = c(0.65, 0.75), p = c(0.35, 0.65))
#' 
#' HR.m20 = overallAHR(S0 = function(t){exp(-log(2)/20*t)}, h0=function(t){log(2)/20}, 
#' HR = c(0.65, 1), p = c(0.70, 0.30))
#' 
#' HR.m20 = overallAHR(S0 = function(t){exp(-log(2)/20*t)}, h0=function(t){log(2)/20}, 
#' HR = c(0.55, 0.65, 1), p = c(0.45, 0.25, 0.30))
#' 
#' overallAHR(S0 = function(t){exp(-log(2)/20*t)}, h0=function(t){log(2)/20}, 
#' HR = c(0.65, 0.8), p = c(0.4, 0.6))
#' 
#' @export
#' 
overallAHR = function(S0 = list(function(t){exp(-log(2)/12*t)}, function(t){exp(-log(2)/12*t)}), 
                      h0=list(function(t){log(2)/12}, function(t){log(2)/12}),
                      HR = c(0.65, 0.75),
                      p = c(0.45, 0.55)){
  #Number of subgroups
  K = length(HR)
  #f0
  f0 = function(t){
    ss = 0
    for (j in 1:K) {ss = ss + S0[[j]](t)*h0[[j]](t)*p[j]}
    return(ss)
  }
  
  #Ratio of S(t) / S0(t), where S(t) = p1*S1+p2*S2+...
  rS = function(t){
    ss = 0
    for (j in 1:length(HR)){ss = ss + p[j]*S0(t)^(HR[j]-1)}
    return(ss)
  }
  hp1 = function(t){
    ss = 0
    for (j in 1:length(HR)){ss = ss + p[j]*S0(t)^(HR[j]-1)/rS(t)*HR[j]}
    return(ss)
  }
  hp0 = function(t){
    ss = 0
    for (j in 1:length(HR)){ss = ss + p[j]*S0(t)^(HR[j]-1)/rS(t)*HR[j]}
    return(ss)
  }
  
  #Average HR
  I0 = function(t){HRp(t)*S0(t)*h0(t)}
  aveHR = integrate(I0, lower=0, upper=1000, abs.tol=1e-8)$value
  
  return(aveHR)
}                    


