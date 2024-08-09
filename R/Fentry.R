#' Creation of Piecewise Enrolment Distribution
#' 
#' This function creates a piecewise enrolment distribution to faciliate the
#' calculation of expected number of events.  
#' 
#' @param t  Enrolment time
#' @param time A vector of time points to define the piecewise enrolment intervals
#' for example, time = c(1,2,3) defines enrolment intervals 0-1, 1-2, 2-3.
#' @param cum.p A vector of cumulative enrolment at each time point. It must have 
#' the same length of time. For example, cum.p = c(0.2, 0.5, 1.0).
#' @param xi A vector of weights for each piecewise enrolment distribution function.
#' It must have the same length as time and cum.p. For example, xi = c(1,1,1)
#' means uniform enrollment in all intervals.
#'  
#' @examples
#' 
#' t = seq(0, 8, by=0.05)
#' fe1 = Fentry(t=t, time=1:5, cum.p = c(0.05, 0.15, 0.35, 0.65, 1.0),xi = rep(1,5))
#' fe2 = Fentry(t=t, time=1:5, cum.p = c(0.05, 0.15, 0.35, 0.65, 1.0),xi = rep(2.5,5))
#' plot(t, fe1, type="n", xlab="Enrollment time (month)", ylab="Enrolment Distribution Function")
#' lines(t, fe1, col=1, lty=1, lwd=2)  
#' lines(t, fe2, col=2, lty=2, lwd=2)  
#' legend(0, 1, c("Piecewise xi=1", "Piecewise xi=2.5"), lwd=rep(2,2), col=1:2, lty=1:2, cex=0.8, bty="n")
#' 
#' @export
#' 
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
