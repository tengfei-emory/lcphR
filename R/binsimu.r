binsimu <- function(mu,gamma){
  len <- length(mu)
  y <- rep(0,len)
  y[1] <- rbinom(1,1,mu[1])
  for (i in 2:len){
    lambda <- mu[i] + gamma*(y[i-1] - mu[i-1])*sqrt((mu[i]*(1-mu[i]))/(mu[i-1]*(1-mu[i-1])))
    y[i] <- rbinom(1,1,lambda)
  }
  y
}

poisimu <- function(mu,v){
  muz <- log(mu^2/sqrt(diag(v)-mu+mu^2))
  covz <- matrix(0,length(muz),length(muz))
  diag(covz) <- log((diag(v)-mu)/(mu^2)+1)
  
  for (i in 1:nrow(covz)){
    for (j in 1:ncol(covz)){
      if(i != j){
        covz[i,j] = log(v[i,j]/(mu[i]*mu[j])+1)
      }
    }
  }
  
  z <- rmvnorm(n=1,mean=muz,sigma=covz)
  z <- exp(z)
  y <- rep(0,length(mu))
  for (i in 1:length(mu)){
    y[i] <- rpois(1,z[i])
  }
  y
}