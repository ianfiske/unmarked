# Fits model Mh - closed population model with heterogeneity
# Assumes balanced data right now. 
# 
#this will be the engine for a poisson reg model for spec richness
Mh<-
function(nx){

##nx<-c(34, 16, 10, 4, 2, 2,0,0,0,0,0,0,0,0)
nind<-sum(nx)
J<-length(nx)

Mhlik<-function(parms){
mu<-parms[1]
sigma<-exp(parms[2])
n0<-exp(parms[3])

il<-rep(NA,J+1)
for(k in 0:J){
il[k+1]<-integrate( 
function(x){
 dbinom(k,J,exp(x)/(1+exp(x)))*dnorm(x,mu,sigma) 
},
lower=-Inf,upper=Inf)$value
}

-1*(    lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(n0,nx)*log(il)))
}


tmp<-nlm(Mhlik,c(0,0,0),hessian=TRUE)

return(tmp)


}


