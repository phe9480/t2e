#'  Power Calculation for Log-rank Tests in Overlapping Populations
#' 
#'  This function calculates the powers at specified analysis times based on the asymptotic
#'  distribution of the log-rank test statistics in overalapping populations under H1.
#'  For group sequential design, the power will be calculated for each analysis
#'  and overall study.
#'  
#' @param T  A vector of analysis times for interim and final analysis, calculated from first subject randomized, .
#' @param incr.alpha A vector of incremental alpha allocated to all analyses. sum(incr.alpha) = overall.alpha. 
#'           If sf is provided, then incr.alpha will be ignored. If sf is not provided,
#'           then incr.alpha is required. In detail, if the alpha spending function a(t) is used,
#'           then incr.alpha=c(a(t1), alpha2 = a(t2)-a(t1), ..., alphaK = a(tK)-a(t_{K-1})
#'           for timing = c(t1, ..., tK=1).
#' @param overall.alpha  Overall familywise one-sided alpha, default 0.025, for both tests.
#' @param sf Spending functions for tests A and B. Default sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF),
#' both tests are based on O'Brien Fleming spending boundary. Refer to gsDesign() function
#' for other choices of spending functions.
#' @param r A vector of proportions of experimental subjects in each of subgroup: Both in A and B, A not B, B not A.
#' For randomization stratified by A and B, then r$AandB = r$AnotB = r$BnotA. 
#' @param strat.ana stratified analysis flag, "Y" or "N". Default, "Y". The stratified analysis 
#' means that testing HA is stratified by B, and vice versa.
#' @param n Total sample size for two arms for subjects in both population A and B, 
#' A not B, B not A. Default is NULL. 
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
#' @param w A vector of proportions for type I error allocation to all primary hypotheses. 
#'          The sum of w must be 1.
#' @param epsA A vector of efficiency factors for testing Ha at all analyses. 
#' @param epsB A vector of efficiency factors for testing Hb at all analyses.
#' epsA and epsB are required when method = "Customized Allocation". At analysis k,
#' either epsA[k] or epsB[k] is required, but not both. The unspecified will be determined. 
#' For example, epsA = c(1, NA) and epsB = c(NA, 1): At the 1st analysis, Ha rejection boundary
#' is the same as the boundary based on alpha-splitting method, and the benefit from the correlation
#' is fully captured to Hb to improve its rejection boundary. At the 2nd analysis, Hb's rejection
#' boundary is the same as the alpha-splitting approach, and the benefit from the correlation is
#' fully captured to Ha. In order to insure the improved rejection boundary is no worse than the 
#' alpha-splitting method, epsA and epsB is must be at least 1, and also capped by the acceptable
#' value when the other one is 1.
#' @param method The method for alpha adjustment: "Balanced Allocation" or "Customized Allocation".
#'        "Balanced Allocation" = the adjustment is equally allocated to all 
#'        primary hypotheses. 
#'        "Customized Allocation" = the adjustment is made according to pre-specified levels
#'        for some hypotheses. Default "Balanced Allocation".
#' @param F.entry Distribution function of enrollment. For uniform enrollment, 
#' F.entry(t) = (t/A) where A is the enrollment period, i.e., F.entry(t) = t/A for 0<=t<=A, and 
#' F.entry(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#' survival time is min(survival time, lost-to-follow-up time). Default G.ltfu = 0 (no lost-to-followup)
#' @param variance Option for variance estimate. "H1" or "H0". Default H1, which is usually more conservative than H0.
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item overall.alpha: Family-wise type I error and each test's type I error
#'  \item  events: Number of events by analysis and by treatment arm
#'       \itemize{
#'       \item DCO: Data cutoff time for analysis, calculated from first subject randomized
#'       \item n.events0: Expected events in control arm
#'       \item n.events1: Expected events in experimental arm
#'       \item n.events.total: Expected events in both control and experimental arms
#'       \item n0: number of subjects in control arm
#'       \item n1: number of subjects in experimental arm
#'       \item n.total: total number of subjects
#'       \item maturity0: maturity in control arm, percent of events
#'       \item maturity1: maturity in experimental arm, percent of events
#'       \item maturity: maturity in both arms, percent of events
#'       }
#'  \item  power Power calculation including variables: 
#'       \itemize {
#'       \item  timingA: Information fraction for test A
#'       \item  marg.powerA:  Marginal power for test A regardless of the previous test results
#'       \item  incr.powerA: Incremental power for test A
#'              The sum of all incremental powers is the overall power.
#'       \item  cum.powerA:  Cumulative power for test A      
#'       \item  overall.powerA: Overall power for test A
#'       \item  marg.powerA0:  Marginal power for test A regardless of the previous test results using alpha-splitting method
#'       \item  incr.powerA0: Incremental power for test A using alpha-splitting method
#'              The sum of all incremental powers is the overall power.
#'       \item  cum.powerA0:  Cumulative power for test A using alpha-splitting method      
#'       \item  overall.powerA0: Overall power for test A using alpha-splitting method
#'       \item  timingB: Information fraction for test B
#'       \item  marg.powerB:  Marginal power for test B regardless of the previous test results
#'       \item  incr.powerB: Incremental power for test B
#'              The sum of all incremental powers is the overall power.
#'       \item  cum.powerB:  Cumulative power for test B     
#'       \item  overall.powerB: Overall power for test B
#'       \item  marg.powerB0:  Marginal power for test B regardless of the previous test results using alpha-splitting method
#'       \item  incr.powerB0: Incremental power for test B using alpha-splitting method
#'              The sum of all incremental powers is the overall power.
#'       \item  cum.powerB0:  Cumulative power for test B using alpha-splitting method     
#'       \item  overall.powerB0: Overall power for test B using alpha-splitting method
#'       }
#'  \item  bd: Rejection boundary in z value and p value including variables
#'       \itemize{
#'       \item timingA: Information fraction for test A
#'       \item incr.alphaA: Incremental alpha for test A
#'       \item cum.alphaA:  Cumulative alpha for test A
#'       \item bd.pA0: p value boundary for test A based on alpha-splitting method
#'       \item bd.zA0: z value boundary for test A based on alpha-splitting method
#'       \item bd.pA: p value boundary for test A
#'       \item bd.zA: z value boundary for test A
#'       \item epsA:   Efficiency factor for test A with correlation considered. Larger epsA indicates more improvement from the alpha-splitting method.
#'       \item timingB: Information fraction for test B
#'       \item incr.alphaB: Incremental alpha for test B
#'       \item cum.alphaB:  Cumulative alpha for test B
#'       \item bd.pB0: p value boundary for test B based on alpha-splitting method
#'       \item bd.zB0: z value boundary for test B based on alpha-splitting method
#'       \item bd.pB: p value boundary for test B
#'       \item bd.zB: z value boundary for test B
#'       \item epsB:   Efficiency factor for test B with correlation considered. Larger epsA indicates more improvement from the alpha-splitting method.
#'       }
#'  \item  CV: Critical Value in HR and median
#'  \item median: Medians by treatment arm (0 or 1) and test (A or B)
#'  \item max.eps: Upper bound of acceptable epsA and epsB values when method = "Customized Allocation".
#'  \item corr:   Correlation matrix of log-rank test statistics vector for K
#'  analyses: (zA1, zB1, zA2, zB2, ..., zAK, zBK)
#'  \item cov:   Covariance matrix of log-rank test score statistics vector for K
#'  analyses: (uA1, uB1, uA2, uB2, ..., uAK, uBK)
#'  \item method: method of improvement allocation
#'  \item strat.ana: stratified analysis flag (Y/N)
#'  }
#'  
#' @examples
#'  
#' #Example: 1:1 randomization, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming 3\% drop-off per 12 months of followup.
#' #250 PD-L1+ subjects, a total of 600 subjects.
#' #3 Analyses are planned: 24 mo, 36 mo, and 42 mo.
#' #Assumed HR: 0.60 for PD-L1+, and 0.80 for PD-L1-, so the HR for overall population is 0.71.
#' 
#' pow = corrPower(T = c(24, 36), n = list(AandB = 350, AnotB=0, BnotA=240), 
#'            r = list(AandB=1/2, AnotB =0, BnotA = 1/2), 
#'            sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'            h0=list(AandB=function(t){log(2)/12}, AnotB=function(t){log(2)/12}, 
#'                    BnotA=function(t){log(2)/12}), 
#'            S0=list(AandB=function(t){exp(-log(2)/12*t)}, AnotB=function(t){exp(-log(2)/12*t)},
#'                    BnotA=function(t){exp(-log(2)/12*t)}),
#'            h1=list(AandB=function(t){log(2)/12*0.6},    AnotB=function(t){log(2)/12*0.6},
#'                    BnotA=function(t){log(2)/12*0.80}), 
#'            S1=list(AandB=function(t){exp(-log(2)/12 * 0.6 * t)},AnotB=function(t){exp(-log(2)/12 * 0.6 * t)},
#'                    BnotA=function(t){exp(-log(2)/12 * 0.80 * t)}),
#'            strat.ana=c("Y", "N"),
#'            alpha=0.025, w=c(1/3, 2/3), epsilon = list(epsA = c(NA,NA), epsB=c(1,1)),
#'            method=c("Balanced Allocation", "Customized Allocation"), 
#'            F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
#'            G.ltfu = function(t){1-exp(-0.03/12*t)}, variance="H1")
#' 
#' @export
#' 
#' 
#' 
corrPower = function(T = c(24, 36), n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), 
      sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
      h0=list(AandB=function(t){log(2)/12}, AnotB=function(t){log(2)/12}, 
              BnotA=function(t){log(2)/12}), 
      S0=list(AandB=function(t){exp(-log(2)/12*t)}, AnotB=function(t){exp(-log(2)/12*t)},
      BnotA=function(t){exp(-log(2)/12*t)}),
      h1=list(AandB=function(t){log(2)/12*0.70},    AnotB=function(t){log(2)/12*0.70},
              BnotA=function(t){log(2)/12*0.70}), 
      S1=list(AandB=function(t){exp(-log(2)/12 * 0.7 * t)},AnotB=function(t){exp(-log(2)/12 * 0.7 * t)},
              BnotA=function(t){exp(-log(2)/12 * 0.7 * t)}),
      strat.ana=c("Y", "N"),
      alpha=0.025, w=c(1/3, 2/3), epsilon = list(epsA = c(NA,NA), epsB=c(1,1)),
      method=c("Balanced Allocation", "Customized Allocation"),
      F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}, 
      G.ltfu = function(t){1-exp(-0.03/12*t)}, variance="H1"){
  
  ########################
  #(1) General settings
  ########################
  
  #(1.1) Number of analyses
  K = length(T)
  #(1.2) rA and rB: proportions of exp. pts in A and B.
  nA = n$AandB+n$AnotB; 
  nB = n$AandB+n$BnotA; 
  N = n$AandB + n$AnotB + n$BnotA
  
  rA = (r$AandB*n$AandB+r$AnotB*n$AnotB)/nA
  rB = (r$AandB*n$AandB+r$BnotA*n$BnotA)/nB
  gamma =  n$AandB / N
  
  #(1.3) Mixture distributions
  qA = n$AandB / nA
  qB = n$AandB / nB
  if (n$AnotB > 0){
    S0A = function(t){qA*S0$AandB(t)+(1-qA)*S0$AnotB(t)}
    S1A = function(t){qA*S1$AandB(t)+(1-qA)*S1$AnotB(t)}
    f0A = function(t){qA*S0$AandB(t)*h0$AandB(t)+(1-qA)*S0$AnotB(t)*h0$AnotB(t)}
    f1A = function(t){qA*S1$AandB(t)*h1$AandB(t)+(1-qA)*S1$AnotB(t)*h1$AnotB(t)}
  } else {
    S0A = S0$AandB; S1A = S1$AandB
    f0A = function(t){S0$AandB(t)*h0$AandB(t)}
    f1A = function(t){S1$AandB(t)*h1$AandB(t)}
  }
  h0A = function(t){f0A(t) / S0A(t)}
  h1A = function(t){f1A(t) / S1A(t)}
  
  if (n$BnotA > 0){
    S0B = function(t){qB*S0$AandB(t)+(1-qB)*S0$BnotA(t)}
    S1B = function(t){qB*S1$AandB(t)+(1-qB)*S1$BnotA(t)}
    f0B = function(t){qB*S0$AandB(t)*h0$AandB(t)+(1-qB)*S0$BnotA(t)*h0$BnotA(t)}
    f1B = function(t){qB*S1$AandB(t)*h1$AandB(t)+(1-qB)*S1$BnotA(t)*h1$BnotA(t)}
  } else {
    S0B = S0$AandB; S1B = S1$AandB
    f0B = function(t){S0$AandB(t)*h0$AandB(t)}
    f1B = function(t){S1$AandB(t)*h1$AandB(t)}
  }  
  h0B = function(t){f0B(t) / S0B(t)}
  h1B = function(t){f1B(t) / S1B(t)}
  
  #(2) Determine number of events based on T
  E = list()
  eAandB = eAnotB = eBnotA = eA = eB = NULL
  for (k in 1:K){
    E[[k]] = tmp = corrEvents(T = T[k], n = n,r = r,h0 = h0, S0 = S0,
                            h1=h1, S1=S1, F.entry = F.entry, G.ltfu = G.ltfu)
    eAandB = c(eAandB, tmp$n.events.total[tmp$subgroup == "AandB"])
    eAnotB = c(eAnotB, tmp$n.events.total[tmp$subgroup == "AnotB"])
    eBnotA = c(eBnotA, tmp$n.events.total[tmp$subgroup == "BnotA"])
    eA = c(eA, tmp$n.events.total[tmp$subgroup == "A"])
    eB = c(eB, tmp$n.events.total[tmp$subgroup == "B"])
  }

  #(3) Determine the rejection boundary incorporating the correlation
  bd = corrBounds(sf=sf, eAandB = eAandB, eAnotB = eAnotB, eBnotA = eBnotA,
             r=r, rA=rA, rB=rB, gamma = gamma,
             strat.ana=strat.ana, alpha=alpha, w=w, epsA = epsilon$epsA, 
             epsB=epsilon$epsB, method=method)
  
  #(4) Non-centrality mean parameter using Schoenfeld method
  
  #(4.1)Pooled density function per stratum
  f.AandB = function(t){
    (1-r$AandB)*h0$AandB(t)*S0$AandB(t) + r$AandB * h1$AandB(t)*S1$AandB(t)
  }
  f.AnotB = function(t){
    (1-r$AnotB)*h0$AnotB(t)*S0$AnotB(t) + r$AnotB * h1$AnotB(t)*S1$AnotB(t)
  }
  f.BnotA = function(t){
    (1-r$BnotA)*h0$BnotA(t)*S0$BnotA(t) + r$BnotA * h1$BnotA(t)*S1$BnotA(t)
  }
  
  #eta and V in K analyses
  ##################
  #(4.2) Stratified 
  ##################
  if (strat.ana[1] == "Y") {
    eta.AandB = eta.AnotB = eta.BnotA = V.AnotB=V.AandB=V.BnotA=rep(0, K)
    for(k in 1:K){
      #(4.2.1) eta functions
      f.eta.AandB = function(t){
        return(log(h1$AandB(t)/h0$AandB(t))*F.entry(T[k]-t) * (1 - G.ltfu(t)) * f.AandB(t))
      }
      f.eta.AnotB = function(t){
        return(log(h1$AnotB(t)/h0$AnotB(t))*F.entry(T[k]-t) * (1 - G.ltfu(t)) * f.AnotB(t))
      }
      f.eta.BnotA = function(t){
        return(log(h1$BnotA(t)/h0$BnotA(t))*F.entry(T[k]-t) * (1 - G.ltfu(t)) * f.BnotA(t))
      }
      #(4.2.2) Calculate eta and V
      if(n$AandB>0){
        I.AandB = integrate(f.eta.AandB, lower=0, upper=T[k], abs.tol=1e-8)$value
        if(variance == "H0"){
          eta.AandB[k] = -n$AandB * F.entry(T[k])*r$AandB*(1-r$AandB)*I.AandB
          V.AandB[k] = r$AandB*(1-r$AandB)*eAandB[k]
        } else{
          rAandB.k = E[[k]]$n.events1[E[[k]]$subgroup == "AandB"]/eAandB[k]
          eta.AandB[k] = -n$AandB * F.entry(T[k])*rAandB.k*(1-rAandB.k)*I.AandB
          V.AandB[k] = rAandB.k*(1-rAandB.k)*eAandB[k]
        }
      }
      if(n$AnotB>0){
        I.AnotB = integrate(f.eta.AnotB, lower=0, upper=T[k], abs.tol=1e-8)$value
        if(variance == "H0"){
          eta.AnotB[k] = -n$AnotB * F.entry(T[k])*r$AnotB*(1-r$AnotB)*I.AnotB
          V.AnotB[k] = r$AnotB*(1-r$AnotB)*eAnotB[k]
        } else{
          rAnotB.k = E[[k]]$n.events1[E[[k]]$subgroup == "AnotB"]/eAnotB[k]
          eta.AnotB[k] = -n$AnotB * F.entry(T[k])*rAnotB.k*(1-rAnotB.k)*I.AnotB
          V.AnotB[k] = rAnotB.k*(1-rAnotB.k)*eAnotB[k]
        }
      }
      if(n$BnotA>0){
        I.BnotA = integrate(f.eta.BnotA, lower=0, upper=T[k], abs.tol=1e-8)$value
        if(variance == "H0"){
          eta.BnotA[k] = -n$BnotA * F.entry(T[k])*r$BnotA*(1-r$BnotA)*I.BnotA
          V.BnotA[k] = r$BnotA*(1-r$BnotA)*eBnotA[k]
        } else {
          rBnotA.k = E[[k]]$n.events1[E[[k]]$subgroup == "BnotA"]/eBnotA[k]
          eta.BnotA[k] = -n$BnotA * F.entry(T[k])*rBnotA.k*(1-rBnotA.k)*I.BnotA
          V.BnotA[k] = rBnotA.k*(1-rBnotA.k)*eBnotA[k]
        }
      }
    }
    eta.A = eta.AandB + eta.AnotB
    eta.B = eta.AandB + eta.BnotA
    
    V.A = V.AandB + V.AnotB
    V.B = V.AandB + V.BnotA
  
    #(4.2.3) mu
    muA = eta.A / sqrt(V.A)
    muB = eta.B / sqrt(V.B)
  } else {
    ##########################
    #(4.3) Unstratified test
    ##########################
    #(4.3.1) Pooled density function across all strata and both arms
    f.A = function(t){
      if (n$AnotB > 0){ans = qA*f.AandB(t) + (1-qA) * f.AnotB(t)} else {
        ans = f.AandB(t)
      }
      return(ans)
    }
    f.B = function(t){
      if (n$BnotA > 0){ans = qB*f.AandB(t) + (1-qB) * f.BnotA(t)} else {
        ans = f.AandB(t)
      }
      return(ans)
    }
    #(4.3.2) eta and V
    eta.A = eta.B = V.A = V.B = rep(0, K)
    for(k in 1:K){
      #eta functions
      f.eta.A = function(t){
        return(log(h1A(t)/h0A(t))*F.entry(T[k]-t) * (1 - G.ltfu(t)) * f.A(t))
      }
      f.eta.B = function(t){
        return(log(h1B(t)/h0B(t))*F.entry(T[k]-t) * (1 - G.ltfu(t)) * f.B(t))
      }
      I.A = integrate(f.eta.A, lower=0, upper=T[k], abs.tol=1e-8)$value
      I.B = integrate(f.eta.B, lower=0, upper=T[k], abs.tol=1e-8)$value
      
      if(variance == "H0"){
        eta.A[k] = -nA * F.entry(T[k])*rA*(1-rA)*I.A
        eta.B[k] = -nB * F.entry(T[k])*rB*(1-rB)*I.B
        
        V.A[k] = rA*(1-rA)*eA[k]
        V.B[k] = rB*(1-rB)*eB[k]
      } else{
        rA.k = E[[k]]$n.events1[E[[k]]$subgroup == "A"]/eA[k]
        rB.k = E[[k]]$n.events1[E[[k]]$subgroup == "B"]/eB[k]
        eta.A[k] = -nA * F.entry(T[k])*rA.k*(1-rA.k)*I.A
        eta.B[k] = -nB * F.entry(T[k])*rB.k*(1-rB.k)*I.B
        V.A[k] = rA.k*(1-rA.k)*eA[k]
        V.B[k] = rB.k*(1-rB.k)*eB[k]
      }
    }
    #(4.3.3) mu
    muA = eta.A / sqrt(V.A)
    muB = eta.B / sqrt(V.B)
  }

  #########################
  #(5) Power Calculation
  #########################
  
  #(5.1) Marginal Power
  marg.powerA = marg.powerB = marg.powerA0 = marg.powerB0 = rep(NA, K)
  
  #Population A and B
  for (k in 1:K){
    marg.powerA[k] = 1-pnorm(bd$bd$bd.zA[k], mean=muA[k])
    marg.powerB[k] = 1-pnorm(bd$bd$bd.zB[k], mean=muB[k])
    marg.powerA0[k] = 1-pnorm(bd$bd$bd.zA0[k], mean=muA[k])
    marg.powerB0[k] = 1-pnorm(bd$bd$bd.zB0[k], mean=muB[k])
  }
  
  #(5.2) Incremental Power
  ####Correlation matrix within A and within B
  corrA = corrB = matrix(1, nrow=K, ncol=K)
  if(K > 1){
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        corrA[i,j] = corrA[j,i] = eA[i]/eA[j]
        corrB[i,j] = corrB[j,i] = eB[i]/eB[j]
      }
    }  
  }
  incr.powerA = incr.powerB = incr.powerA0 = incr.powerB0 = rep(NA, K)
  cum.powerA = cum.powerB = cum.powerA0 = cum.powerB0 = rep(NA, K)
  
  cum.powerA[1] = incr.powerA[1] = marg.powerA[1]
  cum.powerB[1] = incr.powerB[1] = marg.powerB[1]
  cum.powerA0[1] = incr.powerA0[1] = marg.powerA0[1]
  cum.powerB0[1] = incr.powerB0[1] = marg.powerB0[1]
  
  #Population A and B
  if(K > 1){
    for (k in 2:K){
      incr.powerA[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), bd$bd$bd.zA[k]), 
                                        upper = c(bd$bd$bd.zA[1:(k-1)], Inf), mean=muA[1:k], 
                                        corr = corrA[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      incr.powerB[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), bd$bd$bd.zB[k]), 
                                        upper = c(bd$bd$bd.zB[1:(k-1)], Inf), mean=muB[1:k], 
                                        corr = corrB[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      incr.powerA0[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), bd$bd$bd.zA0[k]), 
                                        upper = c(bd$bd$bd.zA0[1:(k-1)], Inf), mean=muA[1:k], 
                                        corr = corrA[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      incr.powerB0[k] = mvtnorm::pmvnorm(lower = c(rep(-Inf, k-1), bd$bd$bd.zB0[k]), 
                                        upper = c(bd$bd$bd.zB0[1:(k-1)], Inf), mean=muB[1:k], 
                                        corr = corrB[1:k, 1:k], abseps = 1e-8, maxpts=100000)[1]
      cum.powerA[k] = cum.powerA[k-1] + incr.powerA[k]
      cum.powerB[k] = cum.powerB[k-1] + incr.powerB[k]
      cum.powerA0[k] = cum.powerA0[k-1] + incr.powerA0[k]
      cum.powerB0[k] = cum.powerB0[k-1] + incr.powerB0[k]
    }
    overall.powerA = cum.powerA[K]
    overall.powerB = cum.powerB[K]
    overall.powerA0 = cum.powerA0[K]
    overall.powerB0 = cum.powerB0[K]
  }
  
  ###########################################
  #(6) Critical Values (min detectable diff)
  ###########################################
  
  #(6.1) In terms of HR
  #######By default, use variance under H1 for more conservative estimate of power
  if (variance == "H0"){
    cvA = exp(-bd$bd$bd.zA/sqrt(rA*(1-rA)*eA))
    cvB = exp(-bd$bd$bd.zB/sqrt(rB*(1-rB)*eB))
    cvA0 = exp(-bd$bd$bd.zA0/sqrt(rA*(1-rA)*eA))
    cvB0 = exp(-bd$bd$bd.zB0/sqrt(rB*(1-rB)*eB))
  } else {
    #Use unstratified version for approximation of overall HR CV 
    rAe = rBe = rep(NA, K)
    for (k in 1:K){
      rAe[k] = E[[k]]$n.events1[E[[k]]$subgroup == "A"]/eA[k]
      rBe[k] = E[[k]]$n.events1[E[[k]]$subgroup == "B"]/eB[k]
    }
    cvA = exp(-bd$bd$bd.zA/sqrt(rAe*(1-rAe)*eA))
    cvB = exp(-bd$bd$bd.zB/sqrt(rBe*(1-rBe)*eB))
    cvA0 = exp(-bd$bd$bd.zA0/sqrt(rAe*(1-rAe)*eA))
    cvB0 = exp(-bd$bd$bd.zB0/sqrt(rBe*(1-rBe)*eB))
  }
  
  #(6.2) In terms of median
  f.m0A = function(t){S0A(t) - 0.5}
  f.m1A = function(t){S1A(t) - 0.5}
  f.m0B = function(t){S0B(t) - 0.5}
  f.m1B = function(t){S1B(t) - 0.5}
  
  m0A = uniroot(f.m0A, interval= c(1, 100), tol = 1e-8)$root
  m1A = uniroot(f.m1A, interval= c(1, 100), tol = 1e-8)$root
  m0B = uniroot(f.m0B, interval= c(1, 100), tol = 1e-8)$root
  m1B = uniroot(f.m1B, interval= c(1, 100), tol = 1e-8)$root
  
  cv.mA = m0A / cvA
  cv.mB = m0B / cvB
  cv.mA0 = m0A / cvA0
  cv.mB0 = m0B / cvB0

  ########################
  #(7) Output
  ########################
  timingA = bd$bd$timingA; timingB = bd$bd$timingB; 
  power=data.frame(cbind(timingA,marg.powerA,incr.powerA,cum.powerA,overall.powerA,
                         marg.powerA0,incr.powerA0,cum.powerA0,overall.powerA0,
                         timingB,marg.powerB,incr.powerB,cum.powerB,overall.powerB,
                         marg.powerB0,incr.powerB0,cum.powerB0,overall.powerB0))  
  o = list()
  o$overall.alpha = bd$overall.alpha  
  o$events = E
  o$bd = bd$bd
  o$mu = data.frame(cbind(muA, muB))
  o$power = power
  o$median = data.frame(cbind(m0A, m1A, m0B, m1B))
  o$CV = data.frame(cbind(cvA, cvB, cvA0, cvB0, cv.mA, cv.mB, cv.mA0, cv.mB0))
  o$max.eps = bd$max.eps
  o$corr = bd$corr
  o$cov = bd$cov
  method = method[1]; strat = strat.ana[1]; var = variance[1]
  o$options = data.frame(cbind(method, strat, var))
  
  return(o)
}

