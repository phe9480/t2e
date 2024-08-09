#' Rejection Boundaries Of Correlated Tests in Group Sequential Design Using Time-To-Event Endpoints
#' 
#' This function calculates the rejection boundaries in p value 
#' (significance level) and z value in group sequential design based on the alpha spending function
#' for each test using the log-rank test (He et al 2021).
#' 
#' @param sf Spending functions for tests A and B. Default sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF),
#' both tests are based on O'Brien Fleming spending boundary. Refer to gsDesign() function
#' for other choices of spending functions.
#' @param eAandB Number of events in subjects included in both population A and population B.
#' For group sequential design, eAandB is a vector by each analysis.
#' @param eAnotB Number of events in subjects included in population A but not in population B.
#' @param eBnotA Number of events in subjects included in population B but not in population A.
#' @param r A vector of proportions of experimental subjects in each of subgroup: Both in A and B, A not B, B not A.
#' For randomization stratified by A and B, then r$AandB = r$AnotB = BcA. 
#' @param rA   Proportion of experimental subjects among those included in population A.
#' @param rB   Proportion of experimental subjects among those included in population B.
#' @param gamma Proportion of subjects included in both A and B among all subjects in A or B. 
#' rA, rB and gamma are required parameter for unstratified analysis; They are not required for stratified analysis.
#' @param strat.ana stratified analysis flag, "Y" or "N". Default, "Y". The stratified analysis 
#' means that testing HA is stratified by B, and vice versa.
#' @param alpha  Overall one-sided type I error for all primary hypothesis tests. Default 0.025.
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
#'
#' @return An object with values
#'  \itemize{
#'  \item overall.alpha Overall type I error allocated to the family-wise type I error, Ha, and Hb.
#'  \itemize{
#'     \item FW.alpha Family-wise type I error, default 0.025 one-sided.
#'     \item alphaA   Overall alpha for testing A
#'     \item alphaB   Overall alpha for testing B
#'     \item side     Side of test. Always one-sided, 1.
#'  }
#'  \item bd:   Data frame includes the rejection boundaries before and after improvements.
#'  \itemize{
#'     \item timingA Timing of analysis in terms of information fraction for testing A
#'     \item bd.pA0 P value bound for the alpha-splitting approach 
#'     \item bd.zA0 z value bound for the alpha-splitting approach
#'     \item bd.pA  Improved p value bound incorporating the correlation
#'     \item bd.zA  Improved z value bound incorporating the correlation
#'     \item epsA   Efficiency factor for testing A, defined as the ratio of improved p bound and the p bound from alpha-splitting method
#'     \item timingB Timing of analysis in terms of information fraction for testing B
#'     \item bd.pB0 P value bound for the alpha-splitting approach 
#'     \item bd.zB0 z value bound for the alpha-splitting approach
#'     \item bd.pB  Improved p value bound incorporating the correlation
#'     \item bd.zB  Improved z value bound incorporating the correlation
#'     \item epsB   Efficiency factor for testing B, defined as the ratio of improved p bound and the p bound from alpha-splitting method
#'  }
#'  \item max.eps: Data frame of maximum epsA and epsB at each analysis. In order to ensure the 
#'  improved bound is not worse than the alpha-splitting approach, epsA and epsB
#'  should be within the range of \[1, max.epsA\] and \[1, max.epsB\]. max.epsAk is
#'  obtained by passing all improvement of the bound to test A at analysis k while 
#'  keeping the bound for test B same as the alpha-splitting approach. max.epsBk is 
#'  also determined similarly. This range is useful for customized setting of epsA and epsB.
#'  \item corr: Correlation matrix of the logrank test statistics
#'  \item cov: Covariance matrix of the logrank scores for (Z_Ak, Z_Bk, Z_Ak', Z_Bk'), where
#'  k' > k
#'  \item method: Method for allocation of the improvements
#'  \item strat: Flag for stratified analysis "Y" or "N", as user input.
#'  }
#'  @references 
#'  He P, Ni P, Zhang F, and Yu C. Group Sequential Monitoring and Study Design for Time-to-Event Endpoints in Overlapping Populations. Manuscript submitted, 2021
#'  
#' @examples 
#' #Example 1. A subgroup (S) and overall population: at the same DCO, the numbers 
#' #of target events are 100(eS) and 150 (eT) respectively. 1:1 randomization.
#' #1/3 alpha is allocated to S and 2/3 allocated to overall population. 
#' 
#' #(1a) Default method: Balanced allocation; stratified analysis
#' #stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(100), eAnotB = c(0), eBnotA = c(50),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana=c("Y", "N"),alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method=c("Balanced Allocation", "Customized Allocation"))
#' 
#' #(1b) Default method: Balanced allocation; unstratified analysis
#' #For unstratified analysis, gamma is required.
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(100), eAnotB = c(0), eBnotA = c(50),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method=c("Balanced Allocation", "Customized Allocation"))
#'    
#' #(1c) Customized Allocation with the subgroup alpha fixed at the initial alpha
#' #stratified Analysis
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(100), eAnotB = c(0), eBnotA = c(50),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana=c("Y"),alpha=0.025, w=c(1/3, 2/3),epsA = c(1), epsB=c(NA),
#'    method=c("Customized Allocation"))
#' 
#' #(1d) Customized Allocation with the subgroup alpha fixed at the initial alpha
#' #unstratified Analysis
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(100), eAnotB = c(0), eBnotA = c(50),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana=c("N"),alpha=0.025, w=c(1/3, 2/3),epsA = c(1), epsB=c(NA),
#'    method=c("Customized Allocation"))
#'    
#' #####
#' #Example 2. Group sequential design. O'Brien Fleming spending function
#' #is used for both tests Ha and Hb. One IA and FA are performed. The number
#' #of events at IA and FA in each set of patients are: (126, 210), (0,0), (54,90) 
#' #for in A and B, in A not B, in B not A respectively. So the events ratio
#' #at IA for testing Ha is 126/180 = 0.7; and for testing Hb is 210 / 300 = 0.70.
#' #The overall type I error is split as 1/3 alpha and 2/3 alpha.
#' 
#' #(2a) stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana="Y",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method=c("Balanced Allocation", "Customized Allocation"))
#'    
#' ##(2b)unstratified Analysis    
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method="Balanced Allocation")
#'    
#' #(2c) Improve Ha only. stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana="Y",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method="Customized Allocation")
#'    
#' #(2d)Improve Ha only. unstratified Analysis    
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,NA), epsB=c(1,1),
#'    method="Customized Allocation")
#'    
#' #(2e) Improve Hb only. stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana="Y",alpha=0.025, w=c(1/3, 2/3),epsA = c(1,1), epsB=c(NA,NA),
#'    method="Customized Allocation")
#'    
#' #(2f)Improve Hb only. unstratified Analysis    
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(1,1), epsB=c(NA,NA),
#'    method="Customized Allocation")
#'    
#' #(2g) Improve Ha at IA and improve Hb at FA. stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana="Y",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,1), epsB=c(1,NA),
#'    method="Customized Allocation")
#'    
#' #(2h)Improve Ha at IA and improve Hb at FA. unstratified Analysis    
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(NA,1), epsB=c(1,NA),
#'    method="Customized Allocation")
#'    
#' #(2i) Improve Hb at IA and improve Ha at FA. stratified Analysis by default
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = NA,
#'    strat.ana="Y",alpha=0.025, w=c(1/3, 2/3),epsA = c(1,NA), epsB=c(NA,1),
#'    method="Customized Allocation")
#'    
#' #(2j)Improve Hb at IA and improve Ha at FA. unstratified Analysis    
#' corrBounds(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
#'    eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
#'    r=list(AandB = 1/2, AnotB=0, BnotA=1/2), rA=1/2, rB=1/2, gamma = 0.8,
#'    strat.ana="N",alpha=0.025, w=c(1/3, 2/3),epsA = c(1,NA), epsB=c(NA,1),
#'    method="Customized Allocation")
#'    
#' @export
corrBounds = function(sf=list(sfuA=gsDesign::sfLDOF, sfuB=gsDesign::sfLDOF), 
               eAandB = c(126, 210), eAnotB = c(0,0), eBnotA = c(54, 90),
               r=list(AandB = 1/2, AnotB=0, BnotA=1/2), 
               rA=1/2, rB=1/2, gamma = 0.8, strat.ana=c("Y", "N"),
      alpha=0.025, w=c(1/3, 2/3), epsA = c(NA,NA), epsB=c(1,1),
      method=c("Balanced Allocation", "Customized Allocation")){

  #Number of analyses
  M = length(eAandB)

  #Events for A and B GSD analyses
  eA = eAandB + eAnotB; eB = eAandB + eBnotA; e = eAandB + eAnotB + eBnotA

  if (M > 1) {  
    #z and p boundaries based on alpha-splitting approach
    dA0 = gsDesign::gsDesign(k=M,alpha=w[1]*alpha,timing=eA/eA[M],sfu=sf$sfuA, test.type=1)$upper
    zA0 <- dA0$bound   #z boundary
    pA0 = 1-pnorm(zA0) #p-value boundary
    aA0 = dA0$spend    #incremental alpha spend

    dB0 = gsDesign::gsDesign(k=M,alpha=w[2]*alpha,timing=eB/eB[M],sfu=sf$sfuB, test.type=1)$upper
    zB0 <- dB0$bound
    pB0 = 1-pnorm(zB0)
    aB0 = dB0$spend
  } else {
    aA0 = pA0 = w[1]*alpha   #p-value boundary
    zA0 = qnorm(1 - pA0)

    aB0 = pB0 = w[2]*alpha   #p-value boundary
    zB0 = qnorm(1 - pB0)
  }
  
  #sigma2 by stratum: for stratified analysis
  s2AB = r$AandB*(1-r$AandB)*eAandB; s2AcB = r$AnotB*(1-r$AnotB)*eAnotB; s2BcA = r$BnotA*(1-r$BnotA)*eBnotA
  strat.s2A = s2AB + s2AcB; strat.s2B = s2AB + s2BcA; 

  #sigma2 by population: for unstratified analysis
  unstrat.s2A = rA*(1-rA)*eA; unstrat.s2B = rB*(1-rB)*eB
  sAB = (rA*rB+(1-rA-rB)*r$AandB)*gamma*e
  
  #Find the correlation matrix OMEGA (2M X 2M) dimension. 
  #Arranged by block of analysis; cov[1:2, 1:2] is for 1st analysis; 
  #cov[3:4,3:4] is for 2nd analysis; cov[5:6, 5:6] is for 3rd analysis, etc.
  #cov[1:2, 3:4] is for between 1st and 2nd in group sequential testing within and between Za and Zb.
  
  #cov: covariance matrix; corr: correlation matrix
  
  cov = corr = matrix(NA, nrow=M*2, ncol=M*2)
  
  #(1)stratified analysis; 
  if (strat.ana[1] == "Y"){
    #For each analysis
    for (j in 1:M) {
      #diag: jth analysis covariance
      cov[2*j-1,2*j-1] = strat.s2A[j]; cov[2*j, 2*j] = strat.s2B[j]
      cov[2*j-1, 2*j] = cov[2*j, 2*j-1] = s2AB[j]
      #off-diag: jth analysis and kth analysis covariance
      if(M > 1 && j+1 <= M){
        for(k in (j+1):M){
          cov[2*j-1,2*k-1] = cov[2*k-1,2*j-1] = strat.s2A[j]; 
          cov[2*j-1,2*k] = cov[2*k, 2*j-1] = s2AB[j]
         
          cov[2*j, 2*k-1] = cov[2*k-1,2*j] = s2AB[j]; 
          cov[2*j, 2*k] = cov[2*k,2*j] = strat.s2B[j] 
        }
      }
    }
  } else {

    #(2)unstratified analysis
    #For each analysis
    for (j in 1:M) {
      #diag: jth analysis covariance
      cov[2*j-1,2*j-1] = unstrat.s2A[j]; cov[2*j, 2*j] = unstrat.s2B[j]
      cov[2*j-1, 2*j] = cov[2*j, 2*j-1] = sAB[j]

      #off-diag: jth analysis and kth analysis covariance
      if(M > 1 && j+1 <= M){
        for(k in (j+1):M){
          cov[2*j-1,2*k-1] = cov[2*k-1,2*j-1] = unstrat.s2A[j]; 
          cov[2*j-1,2*k] = cov[2*k, 2*j-1] = sAB[j]
          
          cov[2*j, 2*k-1] = cov[2*k-1,2*j] = sAB[j]; 
          cov[2*j, 2*k] = cov[2*k,2*j] = unstrat.s2B[j]; 
        }
      }
    }
  }
  
  for (i in 1:(2*M)){for (j in i:(2*M)){
    corr[i, j] = corr[j, i] = cov[i,j]/sqrt(cov[i,i]*cov[j,j])
  }}

  #Total alpha spending in GSD for both tests
  a0 = aA0 + aB0

  #Define eta function (Refer to the manuscript)
  f.eta = function(epsA, epsB){
    L = length(epsA)
    if (L==1){
      u = qnorm(c(1-epsA[1]*pA0[1], 1-epsB[1]*pB0[1]))
      I1 = mvtnorm::pmvnorm(lower = rep(-Inf, 2), upper = u, 
                            corr = corr[1:2, 1:2], abseps = 1e-8, maxpts=100000)[1]
      return(1 - I1 - a0[1])
    }
    if(L > 1){
      pu = NULL
      for (i in 1:L){pu = c(pu, c(1-epsA[i]*pA0[i], 1-epsB[i]*pB0[i]))}
      u = qnorm(pu)
      I1 = mvtnorm::pmvnorm(lower = rep(-Inf, 2*(L-1)), upper = u[1:2*(L-1)], 
                            corr = corr[1:(2*(L-1)), 1:(2*(L-1))], abseps = 1e-8, maxpts=100000)[1]
      I2 = mvtnorm::pmvnorm(lower = rep(-Inf, 2*L), upper = u[1:(2*L)], 
                            corr = corr[1:(2*L), 1:(2*L)], abseps = 1e-8, maxpts=100000)[1]
      return(I1 - I2 - a0[L])
    }      
  }

  #Rename epsA and espB
  epsA0 = epsA; epsB0 = epsB
  epsA = epsB = max.epsA = max.epsB = rep(NA, M)

  #Max epsA and epsB
  for (l in 1:M){
    if (l == 1){epsA.pre = NULL} else {epsA.pre = max.epsA[1:(l-1)]}
    if (l == 1){epsB.pre = NULL} else {epsB.pre = max.epsB[1:(l-1)]}
    
    f.epsA.max = function(eps){ f.eta(epsA=c(epsA.pre,eps), epsB=rep(1, l)) }
    max.epsA[l] = uniroot(f=f.epsA.max, interval=c(1, 10), tol=1e-8)$root
    
    f.epsB.max = function(eps){ f.eta(epsA=rep(1, l), epsB=c(epsB.pre,eps)) }
    max.epsB[l] = uniroot(f=f.epsB.max, interval=c(1, 10), tol=1e-8)$root
  }
  
  #Recursively solve for the rejection boundary for each analysis
  eps = rep(NA, M) #vector to solve
  for (l in 1:M){
    #Solve for eps according to the specified method
    if (l == 1){
      if (method[1] == "Balanced Allocation"){
        #Define the efficiency factor function
        f.eps = function(eps){f.eta(epsA=eps, epsB=eps)}
        eps[1] = uniroot(f=f.eps, interval=c(1, 10), tol=1e-8)$root
        #update epsA and epsB
        epsA[1] = eps[1]; epsB[1] = eps[1]
      }else if (method[1] == "Customized Allocation"){
        f.eps = function(eps){
          epsAi = ifelse(!is.na(epsA0[1]), epsA0[1], eps)
          epsBi = ifelse(!is.na(epsA0[1]) || (is.na(epsA0[1]) && is.na(epsB0[1])), eps, epsB0[1])
          f.eta(epsA=epsAi, epsB=epsBi)
        }
        eps[1] = uniroot(f=f.eps, interval=c(1, 10), tol=1e-8)$root
        #update epsA and epsB
        epsA[1] = ifelse(!is.na(epsA0[1]), epsA0[1], eps[1]);
        epsB[1] = ifelse(!is.na(epsA0[1]) || (is.na(epsA0[1]) && is.na(epsB0[1])), eps[1], epsB0[1])
      }
    } else {
      if (method[1] == "Balanced Allocation"){
        #Define the efficiency factor function
        f.eps = function(eps){f.eta(epsA=c(epsA[1:(l-1)],eps), epsB=c(epsB[1:(l-1)],eps))}
        eps[l] = uniroot(f=f.eps, interval=c(1, 10), tol=1e-8)$root
        #update epsA and epsB
        epsA[l] = eps[l]; epsB[l] = eps[l]
      }else if (method[1] == "Customized Allocation"){
        f.eps = function(eps){
          epsAi = ifelse(!is.na(epsA0[l]), epsA0[l], eps)
          epsBi = ifelse(!is.na(epsA0[l]) || (is.na(epsA0[l]) && is.na(epsB0[l])), eps, epsB0[l])
          f.eta(epsA=c(epsA[1:(l-1)],epsAi), epsB=c(epsB[1:(l-1)],epsBi))
        }
        eps[l] = uniroot(f=f.eps, interval=c(1, 10), tol=1e-8)$root
        #update epsA and epsB
        epsA[l] = ifelse(!is.na(epsA0[l]), epsA0[l], eps[l]);
        epsB[l] = ifelse(!is.na(epsA0[l]) || (is.na(epsA0[l]) && is.na(epsB0[l])), eps[l], epsB0[l])
      }
    }
  }  
    
    bd.pA = epsA * pA0
    bd.pB = epsB * pB0
    bd.zA = qnorm(1-bd.pA)
    bd.zB = qnorm(1-bd.pB)
    timingA = eA/eA[M]
    timingB = eB/eB[M]
    
    bd.pA0 = pA0
    bd.pB0 = pB0
    bd.zA0 = qnorm(1-bd.pA0)
    bd.zB0 = qnorm(1-bd.pB0)
    incr.alphaA = aA0
    incr.alphaB = aB0
    
    cum.alphaA = cum.alphaB = rep(NA, M)
    for (j in 1:M){
      cum.alphaA[j] = sum(aA0[1:j])
      cum.alphaB[j] = sum(aB0[1:j])
    }
    
    o=list()
    FW.alpha = alpha
    alphaA = alpha * w[1]
    alphaB = alpha * w[2]
    side = 1
    
    overall.alpha = data.frame(cbind(FW.alpha, alphaA, alphaB, side))
    bd = data.frame(cbind(timingA, incr.alphaA, cum.alphaA, bd.pA0, bd.zA0, epsA, bd.pA, bd.zA,
                          timingB, incr.alphaB, cum.alphaB, bd.pB0, bd.zB0, epsB, bd.pB, bd.zB))
    max.eps = data.frame(cbind(max.epsA, max.epsB))
    
    o$overall.alpha = overall.alpha
    o$bd = bd
    o$max.eps = max.eps
    o$corr = corr
    o$cov = cov
    o$method = method[1]
    o$strat = strat.ana[1]
    return(o)
}

