# PRSS function generates a percentile ranked set sample of size n=s*r from
# bivariate normal distribution and 
# type-1 Gumbel's bivarite exponential distribution functions

PRSS<-function(s,c,ro,p,MU,VAR,Theta="",Dist="Bnorm"){
  library(MASS)    # mvtnorm 
  # s : set size 
  # c: cycle
  # ro is correlation coefficient of bivariate normal distribution
  # p is percentile value
  # MU is a vector includes MUx adn MUy
  # VAR is a vector includes VARx and VARy
  # Theta is association parameter of type-1 Gumbel's bivariate exponential distribution
  # Dist can be select Bnorm (for bivariate normal distribution) and
  # Bexp (for type-1 Gumbel's bivariate exponential distribution)
  q<- 1-p
  
  if (((s+1)*p)<=0.5){
    # When ((s+1)*p)<0.5, then r=0 and PRSS function gives error. 
    #This if function solves this problem.
    
    r<- 1  # measure rth ranked unit in the set
    rr<-s  # measure rrth ranked unit in the set
  }
  else {
    r<- round((s+1)*p)
    rr<- round((s+1)*q)
  }
  
  PRSSx<-matrix(nrow = c,ncol = s) # This matrix includes selected concomitant variable X using PRSS
  OsetsX<- matrix(nrow = s,ncol = s) # This matrix includes ordered Sets X (Each row is a set)
  PRSSy<-matrix(nrow = c,ncol = s) # This matrix includes measured interested variable Y using PRSS
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
      PRSSy[j,] <- c(OsetsY[1:((nrow(OsetsY))/2),r], OsetsY[(((nrow(OsetsY))/2)+1):nrow(OsetsY),rr])
      PRSSx[j,] <- c(OsetsX[1:((nrow(OsetsX))/2),r], OsetsX[(((nrow(OsetsX))/2)+1):nrow(OsetsX),rr])      
      Ranks<-c(rep(r,((nrow(OsetsY))/2)),rep(rr,((nrow(OsetsY))/2)))     
    }
    
    else {   # When the set size is odd 
      Med<-(nrow(OsetsY)+1)/2
      PRSSy[j,] <- c(OsetsY[1:((nrow(OsetsY)-1)/2),r],OsetsY[(nrow(OsetsY)+1)/2,Med], OsetsY[(((nrow(OsetsY)+1)/2)+1):nrow(OsetsY),rr])
      PRSSx[j,] <- c(OsetsX[1:((nrow(OsetsX)-1)/2),r],OsetsY[(nrow(OsetsX)+1)/2,Med], OsetsY[(((nrow(OsetsX)+1)/2)+1):nrow(OsetsX),rr])  
      Ranks<-c(rep(r,((nrow(OsetsY)-1)/2)),Med,rep(rr,((nrow(OsetsY)-1)/2)))
    }    
    
  }
  Results<-list(PRSSx,PRSSy,Ranks)
  return(Results)
}

PRSS(s = 4,c = 5,p = 0.4,Theta = 0.4,Dist = "Bexp")
PRSS(s = 4,c = 2,ro = 0.7,p = 0.6,MU = c(1,2),VAR = c(2,1),Dist = "Bnorm")
# [[1]] is RSSx, [[2]] is RSSy and [[3]] includes rank values of measured units.