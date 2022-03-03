# SRS function generates a simple random sample of size n from
# bivariate normal distribution and 
# type-1 Gumbel's bivarite exponential distribution functions



SRS<-function(n,ro,MU,VAR,Theta="",Dist="Bnorm"){
  library(MASS)    # mvtnorm 
  # n is sample size
  # ro is correlation coefficient
  # MU is a vector includes MUx adn MUy, respectively.
  # VAR is a vector includes VARx and VARy, respectively.
  # Theta is association parameter of type-1 Gumbel's bivariate exponential distribution
  # Dist can be select Bnorm (for bivariate normal distribution) and
  # Bexp (for type-1 Gumbel's bivariate exponential distribution)
  
  # Note: When Dist="Bnorm", n, ro, MU and VAR are used. On the other hand,
  # n and Theta are used when Dist="Bexp".
  
  if (Dist=="Bnorm"){
    SIGMA <- matrix(nrow=2,ncol=2) # Var-Cov matrix
    ########### Var-Cov matrix #######################################
    for ( i in 1:2) {                                                #
      for ( l in 1:2) {                                              #
        if (i==l){                                                   #
          SIGMA[i,i] <- VAR[i]                                       #
        }                                                            #
        else {                                                       #
          SIGMA[i,l] <- ro*sqrt(VAR[i])*sqrt(VAR[l])                 #
        }                                                            #
      }                                                              #
    }                                                                #
    ##################################################################
  }
  
  Sample<-switch (Dist,
                  Bnorm = mvrnorm(n = n,mu = MU,Sigma = SIGMA), 
                  Bexp = BiExp(n = n,theta = Theta)
  )
  
  
  SRSx<-Sample[,1] # This matrix includes selected concomitant variable X using PRSS
  SRSy<-Sample[,2] # This matrix includes measured interested variable Y using PRSS
  
  Results<-list(SRSx,SRSy)
  return(Results)
}

SRS(n = 20,Theta=1,Dist="Bexp")
SRS(n = 20,ro = 0.7,MU = c(1,2),VAR = c(2,1),Dist = "Bnorm")
# [[1]] is SRSx and [[2]] is SRSy.