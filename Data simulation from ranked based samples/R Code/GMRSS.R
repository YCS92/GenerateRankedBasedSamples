# GMRSS function generates a modified ranked set sample of size n=s*r from
# bivariate normal distribution and 
# type-1 Gumbel's bivarite exponential distribution functions

GMRSS<-function(s,c,ro,R,MU,VAR,Theta="",Dist="Bnorm"){
  library(MASS)    # mvtnorm 
  # s : set size 
  # c: cycle
  # ro is correlation coefficient of bivariate normal distribution
  # R: rank of measured unit R takes a value in the set {1,2,...,s}
  # MU is a vector includes MUx adn MUy
  # VAR is a vector includes VARx and VARy
  # Theta is association parameter of type-1 Gumbel's bivariate exponential distribution
  # Dist can be select Bnorm (for bivariate normal distribution) and
  # Bexp (for type-1 Gumbel's bivariate exponential distribution)
  
  RSSx<-matrix(nrow = c,ncol = s) # This matrix includes selected concomitant variable X using PRSS
  RSSy<-matrix(nrow = c,ncol = s) # This matrix includes measured interested variable Y using PRSS
  
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
  for (j in 1:c) {
    Sample<-switch (Dist,
                    Bnorm = mvrnorm(n = s^2,mu = MU,Sigma = SIGMA), 
                    Bexp = BiExp(n = s^2,theta = Theta)
    ) # Select s^2 units without replacement
    
    X<- matrix(nrow = s,ncol = s) # This matrix includes concomitant variable X
    Y<- matrix(nrow = s,ncol = s) # This matrix includes interested variable Y
    
    for (ii in 1:s) { # Sample is divided into s random sets of size s
      a<- sample(c(1:nrow(Sample)),size = s, replace = F)
      X[ii,]<- Sample[a,1] # In Sample, Second column includes concomitant variable X
      Y[ii,]<- Sample[a,2] # In Sample, First column includes interested variable Y
      RSSy[j,ii]<- Y[ii,order(X[ii,])][R] # Y is ranked from the smallest to the Rth largest using X and measure the Rth (the largest) ordered unit
      RSSx[j,ii]<- X[ii,order(X[ii,])][R] # select x values corresponding to measured Y.
      Sample<- Sample[-a,]
    }
    
  }
  Results<-list(RSSx,RSSy)
  return(Results)
}


GMRSS(s = 3,c = 5,R = 2,Theta = 0.5,Dist = "Bexp")
GMRSS(s = 4,c = 3,ro = 0.7,R = 1,MU = c(1,2),VAR = c(2,1),Dist = "Bnorm")
# [[1]] is RSSx and [[2]] is RSSy