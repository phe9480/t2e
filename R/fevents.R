#' Expected number of events by DCO for each arm with specified distribution
#' 
#' This function calculates the expected number of events for a given enrolment
#' distribution. 
#' 
#' @param DCO  Data cut-off time, calculated from first subject in.
#' @param f  Density function of survival time
#' @param Fentry Enrolment distribution function. Must be non-decreasing valid 
#' distribution function on a fixed interval. For t > enrolment period, Fentry(t)=1
#' @param drop.rate Drop-Off rate in month. If the drop-off rate is 3% per year, 
#' then drop.rate = 0.03/12. 
#' @param n Number of subjects
#'
#'  
#' @examples
#' 
#' #Density function: exponential
#' fa = function(t, lam=log(2)/12){
#' lam*exp(-lam*t)
#' }
#' #Density function: piecewise exponential
#' fb = function(t, lam.a=log(2)/12, lam.b=log(2)/12*0.65, tau = 6){
#'   lam.a*exp(-lam.a*t)*as.numeric(t<6) + as.numeric(t>6) * lam.b*exp(-lam.a*tau)*exp(-lam.b*(t-tau))
#' }
#' #Density function: Mixture cure rate model
#' fc = function(t, k0=0.15, m=12){
#'   lam = log((1-k0)/(0.5-k0)) / m
#'   (1-k0)*lam*exp(-lam*t)
#' }
#' 
#' ####Plot event curves ####
#' DCO = seq(0, 50, by=1); ea = eb = ec = rep(NA, length(DCO))
#' for (i in 1:length(DCO)){
#'   ea[i] = fevents(DCO = DCO[i], f=fa, Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)}, drop.rate = 0.03/12, n = 300)
#'   eb[i] = fevents(DCO = DCO[i], f=fb, Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)}, drop.rate = 0.03/12, n = 300)
#'   ec[i] = fevents(DCO = DCO[i], f=fc, Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)}, drop.rate = 0.03/12, n = 300)
#' }
#' 
#' plot(DCO, ea, type="n", xlab = "Analysis Time (mo)", ylab="Expected number of events", ylim=c(0, max(ea,eb,ec)))
#' lines(DCO, ea, col=1, lwd=2)
#' lines(DCO, eb, col=2, lwd=2)
#' lines(DCO, ec, col=3, lwd=2)
#' legend(0,max(ea,eb,ec),c("Exp: median 12 mo","Piecewise exp","MCR with median 12 mo"), 
#' col=1:3, lwd=2, bty="n", cex=0.8)
#' 
#' @export
#' 
fevents = function(DCO = 40, f=fa, Fentry=function(t){(t/24)^2*as.numeric(t<24)+as.numeric(t>=24)}, drop.rate = 0.03/12, n = 300){
  G = function(t){1-exp(-drop.rate*t)}
  fp = function(t){
    Fentry(DCO - t) * (1-G(t)) * f(t)
  }
  p = integrate(f=fp,lower=0,upper=DCO)$value
  
  return(p*n)
}
