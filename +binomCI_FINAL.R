
binomCI <- function(k_p,
                    n_p,
                    k_se,
                    n_se,
                    k_sp,
                    n_sp,
                    conf = 0.95,
                    alt = "t",
                    toler = 0.0001) {
  # Function to calculate an adjusted Profile Likelihood Confidence Interval for a Prevalence when Sensitivity and / or Specificity is below 100%.
  #  
  # Args:
  #   n_p: The number of individuals that are tested.
  #   k_p: Number of individuals with positive test results in n_p.
  #   n_se: Sample size for a sensitivity validation study.
  #   k_se: Number of individuals with a positive test result in n_se.
  #   n_sp: Sample size for a specificity validation study.
  #   k_sp: Number of individuals with a positive test result in n_sp.
  #   conf.level: Confidence level.
  #   alt: alternative - "two.sided", "less", "greater".
  #   toler: Tolerable imprecision of the confidence limits.
  #   upto_p: Heuristic adjustment used to improve coverage until k_p = upto_p
  #   upto_se: Heuristic adjustment used to improve coverage until k_se = upto_se
  #   upto_sp: Heuristic adjustment used to improve coverage until k_sp = upto_sp
  #
  # Returns:
  #   The observed Sensitivity, Speicificity, Prevalence, the Rogan-Gladen true prevalence and the adjusted PLCI.
  
  PLCI2022 <- function(k_p, 
                       n_p, 
                       k_se, 
                       n_se, 
                       k_sp, 
                       n_sp, 
                       conf.level = 0.95, 
                       alt = "two.sided", 
                       toler = 0.0001) {
    
    # Function to calculate the Profile Likelihood Confidence Interval for a Prevalence when Sensitivity and / or Specificity is below 100%.
    #  
    # Args:
    #   n_p: The number of individuals that are tested.
    #   k_p: Number of individuals with positive test results in n_p.
    #   n_se: Sample size for a sensitivity validation study.
    #   k_se: Number of individuals with a positive test result in n_se.
    #   n_sp: Sample size for a specificity validation study.
    #   k_sp: Number of individuals with a positive test result in n_sp.
    #   conf.level: Confidence level.
    #   alt: alternative - "two.sided", "less", "greater".
    #   toler: Tolerable imprecision of the confidence limits.
    #   upto_p: Heuristic adjustment used to improve coverage until k_p = upto_p
    #   upto_se: Heuristic adjustment used to improve coverage until k_se = upto_se
    #   upto_sp: Heuristic adjustment used to improve coverage until k_sp = upto_sp
    #
    # Returns:
    #   A vector containing the LCL and UCL of the PLCI.
    
    # Observed relative frequencies
    obs.prev = k_p/n_p
    obs.sens = k_se/n_se
    obs.spec = k_sp/n_sp
    
    # Rogan-Gladen point estimate of true prevalence
    est.prev = (obs.prev+obs.spec-1)/(obs.sens+obs.spec-1)
    est.prev = min(1,max(0,est.prev))
    
    alpha <- 1 - conf.level
    alt <- match.arg(alt)
    
    # Warnings
    if (alt != "two.sided" & alt != "less" & alt != "greater") {
      stop("Please select alt from two.sided, less or greater.")
    }
    if (k_p > n_p | k_se > n_se | k_sp > n_sp) {
      stop("k must always be less than or equal to n")
    }
    
    # Likelihood function for parameters p, Se and Sp using their logit
    nlhood_logit_2 = function(theta) {
      p_odds = exp(theta[1]); se_odds = exp(theta[2]); sp_odds = exp(theta[3])
      p=p_odds/(1+p_odds); se=se_odds/(1+se_odds); sp=sp_odds/(1+sp_odds)
      p.app=se*p+(1-sp)*(1-p)
      loglhood_logit = dbinom(k_p, n_p, p.app, log = TRUE) + dbinom(k_se, n_se, se, log = TRUE) + dbinom(k_sp, n_sp, sp, log = TRUE)
      return(-loglhood_logit)
    } 
    # Profile likelihood function
    ProfLik.p <- function(p) {
      nlhood_logit_sesp = function(theta) {
        se_odds = exp(theta[1]); sp_odds = exp(theta[2])
        se=se_odds/(1+se_odds); sp=sp_odds/(1+sp_odds)
        p.app=se*p+(1-sp)*(1-p)
        loglhood_logit = dbinom(k_p, n_p, p.app, log = TRUE) + dbinom(k_se, n_se, se, log = TRUE) + dbinom(k_sp, n_sp, sp, log = TRUE)  ###RJ
        return(-loglhood_logit)
      }
      mle.sesp = optim(par = c(3, 3), fn = nlhood_logit_sesp)	
      return(-mle.sesp$value)
    } 
    
    mlbecsl=optim(par = c(0, 3, 3), fn = nlhood_logit_2, control=list(maxit=10000, reltol=1e-14))
    mlpar_oddsz=exp(mlbecsl$par[1]); mlpar=mlpar_oddsz/(1+mlpar_oddsz)
    
    if (alt == "two.sided") {
      kuszob = mlbecsl$value + qchisq(1-alpha, 1)/2
      
      # search for LCL
      balveg=0; jobbveg=mlpar
      repeat{
        #cat(balveg,jobbveg,"\n")
        ujpont=(balveg+jobbveg)/2
        if(ujpont < toler/2) {LCL=0; break}
        if(jobbveg-balveg < toler) {LCL=balveg; break}
        #print(-ProfLik.p(ujpont))
        if(-ProfLik.p(ujpont) < kuszob) jobbveg=ujpont else balveg=ujpont
      }
      
      # search for UCL
      balveg=mlpar; jobbveg=1
      repeat{
        #cat(balveg,jobbveg,"\n")
        ujpont=(balveg+jobbveg)/2
        if(ujpont > 1 - toler/2) {UCL=1; break}
        if(jobbveg-balveg < toler) {UCL=jobbveg; break}
        #print(-ProfLik.p(ujpont))
        if(-ProfLik.p(ujpont) > kuszob) jobbveg=ujpont else balveg=ujpont
      }
      
      #cat(LCL,UCL,"\n")
      CI=c(LCL,UCL)
    }
    if (alt == "less") {
      
      kuszob = mlbecsl$value + qchisq(1-2*alpha, 1)/2  
      LCL = 0
      
      # search for UCL
      balveg=mlpar; jobbveg=1
      repeat{
        #cat(balveg,jobbveg,"\n")
        ujpont=(balveg+jobbveg)/2
        if(ujpont > 1 - toler/2) {UCL=1; break}
        if(jobbveg-balveg < toler) {UCL=jobbveg; break}
        #print(-ProfLik.p(ujpont))
        if(-ProfLik.p(ujpont) > kuszob) jobbveg=ujpont else balveg=ujpont
      }
      
      CI=c(LCL,UCL)
    }
    if (alt == "greater") {
      
      kuszob = mlbecsl$value + qchisq(1-2*alpha, 1)/2  
      UCL = 1
      
      # search for LCL
      balveg=0; jobbveg=mlpar
      repeat{
        #cat(balveg,jobbveg,"\n")
        ujpont=(balveg+jobbveg)/2
        if(ujpont < toler/2) {LCL=0; break}
        if(jobbveg-balveg < toler) {LCL=balveg; break}
        #print(-ProfLik.p(ujpont))
        if(-ProfLik.p(ujpont) < kuszob) jobbveg=ujpont else balveg=ujpont
      }
      CI=c(LCL,UCL)
    }
    
    return(CI)
    # cat("\nSensitivity: ", round(obs.sens, 4), ", adjusted: ", round(obs.sens., 4),
    #     "\nSpecificity: ", round(obs.spec, 4), ", adjusted: ", round(obs.spec., 4),
    #     "\nObserved prevalence: ", round(obs.prev, 8), ", adjusted: ", round(obs.prev., 8),
    #     "\nRogan-Gladen true prevalence: ", round(est.prev,4), "\n",
    #     "Adjusted ", round(100*conflevel), "% CI: ", round(LCL,4), " - ", round(UCL,4), "\n",
    #     sep="")
  }
  
  # Observed relative frequencies
  obs.prev = k_p/n_p
  obs.sens = k_se/n_se
  obs.spec = k_sp/n_sp
  
  # Rogan-Gladen point estimate of true prevalence
  est.prev = (obs.prev+obs.spec-1)/(obs.sens+obs.spec-1)
  est.prev = min(1,max(0,est.prev))
  
  # Heuristic adjustment used to improve coverage
  upto_p=upto_se=upto_sp=0
  
  # Adjusting the original PLCI 
  CI=PLCI2022(k_p=k_p,n_p=n_p, k_se=k_se,n_se=n_se, k_sp=k_sp, n_sp=n_sp, conf=conf, toler=toler, alt=alt)
  if(k_p<=upto_p) {
    CInext=PLCI2022(k_p=k_p+1,n_p=n_p, k_se=k_se,n_se=n_se, k_sp=k_sp, n_sp=n_sp, conf=conf, toler=toler, alt=alt)
    if(CInext[1]<CI[1])CI[1]=(CI[1]+CInext[1])/2
    if(CInext[2]>CI[2])CI[2]=(CI[2]+CInext[2])/2
  }
  if(k_se>=n_se-upto_se) {
    CInext=PLCI2022(k_p=k_p,n_p=n_p, k_se=k_se-1,n_se=n_se, k_sp=k_sp, n_sp=n_sp, conf=conf, toler=toler, alt=alt)
    if(CInext[1]<CI[1])CI[1]=(CI[1]+CInext[1])/2
    if(CInext[2]>CI[2])CI[2]=(CI[2]+CInext[2])/2
  }
  if(k_sp>=n_sp-upto_sp) {
    CInext=PLCI2022(k_p=k_p+1,n_p=n_p, k_se=k_se,n_se=n_se, k_sp=k_sp-1, n_sp=n_sp, conf=conf, toler=toler, alt=alt)
    if(CInext[1]<CI[1])CI[1]=(CI[1]+CInext[1])/2
    if(CInext[2]>CI[2])CI[2]=(CI[2]+CInext[2])/2
  }
  
  #return(CI)
  cat("\nSensitivity: ", round(obs.sens, 4), 
      "\nSpecificity: ", round(obs.spec, 4), 
      "\nObserved prevalence: ", round(obs.prev, 8), 
      "\nRogan-Gladen true prevalence: ", round(est.prev,4), "\n",
      "Adjusted ", round(100*conf), "% CI: ", round(CI[1],4), " - ", round(CI[2],4), "\n",sep="")
}

binomCI(k_p=711,n_p=11862, k_se=8, n_se=10, k_sp=12, n_sp=12, alt="t", toler=.000001)
