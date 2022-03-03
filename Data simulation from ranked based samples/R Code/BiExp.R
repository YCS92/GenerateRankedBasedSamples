## BiExp function generates data from type-1 Gumbel's bivarite exponential distribution

BiExp<-function(n="",theta=""){
  X<-as.vector(matrix(data = c(rep(0,n)),ncol = n))
  Y<-as.vector(matrix(data = c(rep(0,n)),ncol = n))
  for (i in 1:n){
    while(X[i]==0 & Y[i]==0){
      X[i]<- rexp(n = 1,rate = 1)
      Y[i]<- rexp(n = 1,rate = 1)
      U<-runif(n = 1,min = 0,max = 1)
      
      if(U<=(theta/(1+(theta*X[i])))){
        E<-rexp(n = 1,rate = 1)
        Y[i]<- Y[i]+E
      }
    }
  }
  
  return(cbind(X,(Y/(1+(theta*X)))))
  
}

# Note: BiExp function is requiered for SRS, RSS, ERSS, PRSS and GMRSS 
#R functions when sample data is generated from type-1 Gumbel's 
#bivariate exponential distribution.