#' Expected Number of Events Over Time With Staggered Entry For Overlapping Populations
#' 
#' This function calculates the expected number of events under alternative hypothesis
#' at an analysis time, which is calculated from first subject in. The function
#' returns the expected number of events for each arm, based on the provided
#' enrollment distribution function and drop-off distribution if applicable.
#' If the total sample size is not provided, then only the corresponding probability of event
#' for each arm is provided. The expected number of events for each set of subjects is 
#' calculated including in both A and B, in A not B, in B not A, and in A or B.
#' 
#' @param T  Analysis time calculated from first subject randomization date.
#' @param r  Proportion of experimental subjects in both population A and B, 
#' A not B, B not A. Default, 1/2 for 1:1 randomization stratified by A and B
#' @param h0 Hazard function of control arm for subjects in both population A and B, 
#' A not B, B not A. h0(t) = log(2)/m0 means T~exponential distribution with 
#' median m0. For study design without considering heterogeneous effect in 
#' strata for control arm, then specify the same h0(t) function across strata.
#' @param S0 Survival function of control arm for subjects in both population A and B, 
#' A not B, B not A. h0(t) = log(2)/m0 means T~exponential distribution with 
#' median m0. For study design without considering heterogeneous effect in 
#' strata for control arm, then specify the same S0(t) function across strata.
#' The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm for subjects in both population A and B, 
#' A not B, B not A. For study design without considering heterogeneous effect in 
#' strata for the experimental arm, then specify the same h1(t) function across strata.
#' @param S1 Survival function of experimental arm for subjects in both population A and B, 
#' A not B, B not A. For study design without considering heterogeneous effect in 
#' strata for the experimental arm, then specify the same h1(t) function across strata.
#' @param F.entry Distribution function of enrollment. For uniform enrollment, 
#' F.entry(t) = (t/A) where A is the enrollment period, i.e., F.entry(t) = t/A for 0<=t<=A, and 
#' F.entry(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#' survival time is min(survival time, lost-to-follow-up time). For 3\% drop-off 
#' every year as a constant rate, then G.ltfu = 1-exp(-0.03/12*t), i.e., drop-off~exp distribution.
#' Default G.ltfu = 0 (no lost-to-followup)
#' @param n Total sample size for two arms for subjects in both population A and B, 
#' A not B, B not A. Default is NULL. 
#'
#' @return An object with dataframes below.
#'  \describe{
#'  \item{p.event}{Probability of event for each arm (p.event0: control group; p.event1: experimental group)}
#'  \item{n.events}{Expected number of events}
#'       \itemize{
#'       \item subgroup: Label of subgroups       
#'       \item n.events0: number of events for control group
#'       \item n.events1: number of events for experimental group
#'       \item n.events.total: total number of events for two groups
#'       \item n0: number of subjects in control arm
#'       \item n1: number of subjects in experimental arm
#'       \item n.total: total number of subjects
#'       \item maturity0: maturity in control arm, percent of events
#'       \item maturity1: maturity in experimental arm, percent of events
#'       \item maturity: maturity in both arms, percent of events
#'       }
#'  \item{param}{Parameters specified: T, r, and n}
#'  \item{param.fun}{f0, f1, F.entry, G.ltfu}
#'  }
#'  
#' @examples
#' #PD-L1+ subgroup and overall population tests. Control arm has exponential
#' #distribution with median 12 months in all strata. Experimental arm has
#' #proportional hazards with a hazard ratio of 0.7 in all strata. The randomization
#' #is stratified by PD-L1+ status. 1:1 randomization. 300 subjects in PD-L1+ and 
#' #450 total subjects in overall population. 3\% drop-off every year. Enrollment
#' #period is 18 months and weight 1.5. The expected number of events at 24 months is
#' #calculated as
#' corrEvents(T = 24, n = list(AandB = 300, AnotB=0, BnotA=450), 
#'    r = list(AandB=1/2, AnotB =0, BnotA = 1/2), 
#'    h0=list(AandB=function(t){log(2)/12}, AnotB=function(t){log(2)/12}, 
#'            BnotA=function(t){log(2)/12}), 
#'    S0=list(AandB=function(t){exp(-log(2)/12*t)}, AnotB=function(t){exp(-log(2)/12*t)},
#'            BnotA=function(t){exp(-log(2)/12*t)}),
#'    h1=list(AandB=function(t){log(2)/12*0.70},    AnotB=function(t){log(2)/12*0.70},
#'            BnotA=function(t){log(2)/12*0.70}), 
#'    S1=list(AandB=function(t){exp(-log(2)/12 * 0.7 * t)},AnotB=function(t){exp(-log(2)/12 * 0.7 * t)},
#'            BnotA=function(t){exp(-log(2)/12 * 0.7 * t)}),
#'    F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
#'    G.ltfu = function(t){1-exp(-0.03/12*t)})
#'
#' 
#' @export
#' 
corrEvents = function(T = 24, n = list(AandB = 300, AnotB=0, BnotA=450), 
    r = list(AandB=1/2, AnotB =0, BnotA = 1/2), 
    h0=list(AandB=function(t){log(2)/12}, AnotB=function(t){log(2)/12}, 
            BnotA=function(t){log(2)/12}), 
    S0=list(AandB=function(t){exp(-log(2)/12*t)}, AnotB=function(t){exp(-log(2)/12*t)},
            BnotA=function(t){exp(-log(2)/12*t)}),
    h1=list(AandB=function(t){log(2)/12*0.70},    AnotB=function(t){log(2)/12*0.70},
            BnotA=function(t){log(2)/12*0.70}), 
    S1=list(AandB=function(t){exp(-log(2)/12 * 0.7 * t)},AnotB=function(t){exp(-log(2)/12 * 0.7 * t)},
            BnotA=function(t){exp(-log(2)/12 * 0.7 * t)}),
    F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
    G.ltfu = function(t){0}){

    f.e = function(h00=h0$AandB, h11=h1$AandB, S00=S0$AandB, S11=S1$AandB, nn=n$AandB, rr=r$AandB){
      #Integrand of control arm and experimental arm for calculation of prob. of event
      I0 = function(t){F.entry(T-t) * (1 - G.ltfu(t)) * h00(t) * S00(t)}
      I1 = function(t){F.entry(T-t) * (1 - G.ltfu(t)) * h11(t) * S11(t)}
  
      #prob. of event for control and experimental arm
      p.event0 = integrate(I0, lower=0, upper=T)$value
      p.event1 = integrate(I1, lower=0, upper=T)$value
   
      p.event = data.frame(cbind(p.event0, p.event1))
  
      #expected number of events
  
      if(!is.null(nn)){
        n.events0 = nn * F.entry(T) * (1-rr) * p.event0
        n.events1 = nn * F.entry(T) * rr * p.event1
        n.events.total = n.events0 + n.events1
      }
  
      n.events = data.frame(cbind(n.events0, n.events1, n.events.total))
      ans = list()
      ans$p.event = p.event; ans$n.events = n.events
      return(ans)
    }

    ##Trick the calculation: if no subjects in A not B; set h and S function as 0
    z00 = function(t){0}
    if(n$AnotB == 0){h0$AnotB=h1$AnotB=S0$AnotB=S1$AnotB=z00}
    if(n$BnotA == 0){h0$BnotA=h1$BnotA=S0$BnotA=S1$BnotA=z00}
    
    AandB = f.e(h00=h0$AandB, h11=h1$AandB, S00=S0$AandB, S11=S1$AandB, nn=n$AandB)
    AnotB = f.e(h00=h0$AnotB, h11=h1$AnotB, S00=S0$AnotB, S11=S1$AnotB, nn=n$AnotB)
    BnotA = f.e(h00=h0$BnotA, h11=h1$BnotA, S00=S0$BnotA, S11=S1$BnotA, nn=n$BnotA)
    events = data.frame(rbind(AandB$n.events, AnotB$n.events, BnotA$n.events))
    events[4,] = events[1,]+events[2,]
    events[5,] = events[1,]+events[3,]
    events[6,]= events[1,]+events[2,]+events[3,]

    nA = n$AandB+n$AnotB
    nB = n$AandB+n$BnotA
    nAorB = n$AandB+n$AnotB+n$BnotA
    rA = (n$AandB * r$AandB + n$AnotB*r$AnotB)/nA
    rB = (n$AandB * r$AandB + n$AnotB*r$AnotB + n$BnotA*r$BnotA)/nB
    rAorB = (n$AandB * r$AandB + n$BnotA*r$BnotA)/nAorB
    n.total = c(n$AandB,n$AnotB,n$BnotA, nA,nB, nAorB)
    n1 = n.total*c(r$AandB,r$AnotB,r$BnotA, rA,rB,rAorB)
    n0 = n.total-n1
    
    maturity0 = events$n.events0 / n0
    maturity1 = events$n.events1 / n1
    maturity = events$n.events.total / n.total
    DCO = T
    subgroup = c("AandB", "AnotB", "BnotA", "A", "B", "AorB")
    out = data.frame(cbind(subgroup, DCO, events,n0, n1, n.total, maturity0, maturity1, maturity))
    
  return(out)
}
