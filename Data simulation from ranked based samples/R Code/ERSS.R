# ERSS function generates a extreme ranked set sample of size n=s*r from
# bivariate normal distribution and 
# type-1 Gumbel's bivarite exponential distribution functions

ERSS<-function(s,c,ro,R1,R2,MU,VAR,Theta="",Dist="Bnorm"){
  library(MASS)    # mvtnorm 
  # s : set size 
  # c: cycle
  # ro is correlation coefficient of bivariate normal distribution
  # R1=1 and R2=s for extreme ranked set sampling
  # Note: R1 and R2 can be selected other values in {1,2,...,s}. For example,
  # R1=3, R2=4 when the set size (s) is 4
  # MU is a vector includes MUx adn MUy
  # VAR is a vector includes VARx and VARy
  # Theta is association parameter of type-1 Gumbel's bivariate exponential distribution
  # Dist can be select Bnorm (for bivariate normal distribution) and
  # Bexp (for type-1 Gumbel's bivariate exponential distribution)
  
  ERSSx<-matrix(nrow = c,ncol = s) # This matrix includes selected concomitant variable X using PRSS
  OsetsX<- matrix(nrow = s,ncol = s) # This matrix includes ordered Sets X (Each row is a set)
  ERSSy<-matrix(nrow = c,ncol = s) # This matrix includes measured interested variable Y using PRSS
  OsetsY<- matrix(nrow = s,ncol = s) # This matrix includes ordered Sets Y (Each row is a set)
  
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
      Y[ii,]<- Sample[a,2] # In Sample, First column includes concomitant variable Y
      OsetsY[ii,]<- Y[ii,order(X[ii,])] # Y is ranked from the smallest to the largest using X 
      OsetsX[ii,]<- X[ii,order(X[ii,])] # select x values corresponding to measured Y.
      Sample<- Sample[-a,]
    }
    
    if (nrow(OsetsY) %% 2 == 0) {   # When the set size is even
      ERSSy[j,] <- c(OsetsY[1:((nrow(OsetsY))/2),R1], OsetsY[(((nrow(OsetsY))/2)+1):nrow(OsetsY),R2])
      ERSSx[j,] <- c(OsetsX[1:((nrow(OsetsX))/2),R1], OsetsX[(((nrow(OsetsX))/2)+1):nrow(OsetsX),R2])      
      Ranks<-c(rep(R1,((nrow(OsetsY))/2)),rep(R2,((nrow(OsetsY))/2)))     
    }
    
    else {   # When the set size is odd 
      Med<-(nrow(OsetsY)+1)/2
      ERSSy[j,] <- c(OsetsY[1:((nrow(OsetsY)-1)/2),R1],OsetsY[(nrow(OsetsY)+1)/2,Med], OsetsY[(((nrow(OsetsY)+1)/2)+1):nrow(OsetsY),R2])
      ERSSx[j,] <- c(OsetsX[1:((nrow(OsetsX)-1)/2),R1],OsetsY[(nrow(OsetsX)+1)/2,Med], OsetsY[(((nrow(OsetsX)+1)/2)+1):nrow(OsetsX),R2])  
      Ranks<-c(rep(R1,((nrow(OsetsY)-1)/2)),Med,rep(R2,((nrow(OsetsY)-1)/2)))
    }    
    
  }
  Results<-list(ERSSx,ERSSy,Ranks)
  return(Results)
}

ERSS(s = 4,c = 10,R1 = 1,R2 = 4,Theta = 0.5,Dist = "Bexp")
ERSS(s = 4,c = 2,ro = 0.6,R1 = 1,R2 = 4,MU = c(1,2),VAR = c(2,1),Dist = "Bnorm")
# [[1]] is RSSx, [[2]] is RSSy and [[3]] includes rank values of measured units.