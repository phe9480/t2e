#' Actual Rejection Boundary Based on Actual Number of Events in Group Sequential Design
#' 
#' This function calculates the actual rejection boundary at FA based on the 
#' alpha spending already occurred at interim analyses. The alpha spending for 
#' interim analyses are according to the information fractions relative to
#' the planned target events at final analysis. The rejection boundary at
#' the final analysis is calculated according to the alpha already spent at IAs
#' and the actual number of events observed at final analysis.
#' 
#' @param planned.events  A vector of planned target events for IAs and FA
#' @param actual.events  A vector of actual target events for IAs and FA
#' @param sf Alpha spending function. Default sf = gsDesign::sfLDOF, i.e., 
#' LanDeMets implementation of O'Brien Fleming spending function.
#' @param alpha Overall alpha for the test in the group sequential design. 
#' Default, alpha = 0.025. Alpha must be one-sided.
#'  
#' @return An object with values
#'  \itemize{
#'  \item planned.p Planned rejection bounds in p value (1-sided) based on planned number of events at each analysis
#'  \item planned.z Planned rejection bounds in z value  based on planned number of events at each analysis
#'  \item actual.p Actual rejection bounds in p value (1-sided) based on actual number of events at each analysis and the FA bound is adjusted by the alpha already spent at IAs and the actual number of events at FA
#'  \item actual.z Actual rejection bounds in z value  based on actual number of events at each analysis  and the FA bound is adjusted by the alpha already spent at IAs and the actual number of events at FA
#'  }
#' @examples
#' 
#' actualBounds(planned.events=c(126, 210), act.events=c(140, 250), sf=gsDesign::sfLDOF, alpha=0.025)
#' 
#' @export
#' 
actualBounds = function(planned.events=c(126, 210), act.events=c(140, 250), sf=gsDesign::sfLDOF, alpha=0.025){
  
  #Number of analyses
  M = length(planned.events)
  
  #Planned information fractions
  planned.frac = planned.events / planned.events[M]
  
  if (M == 1){
    p = p0 = alpha
    z = z0 = qnorm(1-alpha)
  } else {
    #Actual information fractions used for IAs alpha spending
    act.frac.IA = act.events[1:(M-1)] / planned.events[M]
    
    #Planned z and p boundaries
    d0 = gsDesign::gsDesign(k=M,alpha=alpha,timing=planned.frac,sfu=sf, test.type=1)$upper
    z0 <- d0$bound   #z boundary
    p0 = 1-pnorm(z0) #p-value boundary
    a0 = d0$spend    #incremental alpha spend
    
    #Actual z and p voundaries for IAs
    act.d = gsDesign::gsDesign(k=M,alpha=alpha,timing=c(act.frac.IA, 1),sfu=sf, test.type=1)$upper
    act.z <- act.d$bound   #z boundary
    act.p = 1-pnorm(act.z) #p-value boundary
    act.a = act.d$spend    #incremental alpha spend
    
    #incremental alpha for FA
    act.a.FA = act.a[M]

    #Actual information fractions used for correlation matrix calculation
    act.frac = act.events / act.events[M]
    
    #correlation matrix
    corr = matrix(1, nrow=M, ncol=M)
    for (i in 1:(M-1)){for (j in (i+1):M){
      corr[i, j] = corr[j, i] = sqrt(act.frac[i]/act.frac[j])
    }}
    
    #Calculate the actual boundary for FA
    f.c = function(c){
      #c is the FA boundary in z value
      I1 = mvtnorm::pmvnorm(lower = c(rep(-Inf, M-1),c), upper = c(act.z[1:(M-1)], Inf), 
                            corr = corr, abseps = 1e-8, maxpts=100000)[1]
      return(I1 - act.a.FA)
    }
    
    #FA bound
    z.FA = uniroot(f=f.c, interval=c(1, 10), tol=1e-8)$root
    p = c(act.p[1:(M-1)], 1-pnorm(z.FA))
    z = c(act.z[1:(M-1)], z.FA)
  }
  out = list()
  out$actual.z = z; out$actual.p = p; out$planned.p=p0; out$planned.z=z0
  return(out)
}  
