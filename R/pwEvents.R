#' Expected number of events for a piecewise enrolment distribution
#' 
#' This function calculates the expected number of events for piecewise enrolment
#' distribution. 
#' 
#' @param DCO  Data cut-off time, calculated from first subject in.
#' @param f  Density function of survival time
#' @param time A vector of time points to define the piecewise enrolment intervals
#' for example, time = c(1,2,3) defines enrolment intervals 0-1, 1-2, 2-3.
#' @param cum.p A vector of cumulative enrolment at each time point. It must have 
#' the same length of time. For example, cum.p = c(0.2, 0.5, 1.0).
#' @param xi A vector of weights for each piecewise enrolment distribution function.
#' It must have the same length as time and cum.p. For example, xi = c(1,1,1)
#' means uniform enrollment in all intervals.
#' @param drop.rate Drop-Off rate in month. If the drop-off rate is 3% per year, 
#' then drop.rate = 0.03/12. 
#' @param n Number of subjects
#'
#'  
#' @examples
#' 
#' pwEvents(DCO = 40, f=fa, time=1:5, cum.p = c(0.05, 0.15, 0.35, 0.65, 1.0), 
#' xi=rep(1,5), drop.rate = 0.03/12, n = 300)
#' 
#' pwEvents(DCO = 40, f=fa, time=5, cum.p = 1, xi=2, drop.rate = 0.03/12, n = 300)
#' 
#' 
#' @export
#' 
pwEvents = function(DCO = 40, f=fa, time=1:5, cum.p = c(0.05, 0.15, 0.35, 0.65, 1.0),
                     xi=rep(1,5), drop.rate = 0.03/12, n = 300){
 #Entry distribution
 Fentry = function(t, time=1:5, cum.p = c(0.05, 0.15, 0.35, 0.65, 1.0),
                    xi=rep(1,5)){
    #To faciliate integrate() function, needs to assume t is a vector for programming
    Lt =length(t)
    E = max(time)
    K = length(time) #number of pieces
    p = cum.p
    p[2:length(p)] = cum.p[2:length(p)]-cum.p[1:(length(p)-1)]
    
    #kth piece entry function
    Fe = function(t, L, R, xi){
      if (t < R){((t-L)/(R-L))^xi} else{1}
    }
    
    #Find the location of t
    loc = ans = rep(NA, Lt)
    for (i in 1:Lt){
      if(t[i] <= time[1]) {loc[i] = 1} else {
        for (k in 2:K){if(t[i] <= time[k] && t[i] > time[k-1]) loc[i] = k}
        if (t[i] > E){loc[i] = K+1}
      }
      
      if(loc[i] == 1){
        ans[i] = p[1]*Fe(t=t[i], L = 0, R = time[1], xi=xi[1])
      }else if (t[i] <= E){
        ans[i] = sum(p[1:(loc[i]-1)]) + p[loc[i]] * Fe(t=t[i], L = time[loc[i]-1], R = time[loc[i]], xi=xi[loc[i]])
      } else{ans[i] = 1}
    }
    return(ans)
 }
 
 #Expected number of events
 fevents = function(DCO = 40, f=fa, Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)}, drop.rate = 0.03/12, n = 300){
   G = function(t){1-exp(-drop.rate*t)}
   fp = function(t){
     Fentry(DCO - t) * (1-G(t)) * f(t)
   }
   p = integrate(f=fp,lower=0,upper=DCO)$value
   
   return(p*n)
 }
 
  if (length(time)==1){
    e=fevents(DCO = DCO, f=fa, 
              Fentry=function(t){(t/time[1])^xi[1]*as.numeric(t<time[1])+as.numeric(t>=time[1])}, 
              drop.rate = drop.rate, n = n)
  } else {
    Fe = function(t){
      Fentry(t, time=time, cum.p = cum.p, xi=xi)
    }
    e = fevents(DCO = DCO, f=f, Fentry=Fe, drop.rate = drop.rate, n = n)   
  }  
  return(e)
}
