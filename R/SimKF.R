# All must be matrix data type!
#X  : T   by s   matrix representing regressors
#Gam: p+q by q   matrix representing parameters of states (H' F')'
#Bet: p+q by s   matrix ... parameters of regressors (A' B')'
#G  : p+q by p+q matrix ... parameters of errors (R' Q')'
#h0 : q   by 1   matrix ... initial state
#X0 : s   by 1   matrix ... x at time step 0
#by Y. K. Yang
#output
#Z  : T   by p   observations
#h  : T   by q   underlying process
SimKF<-function(X,Gam,Bet,G,h0,X0)
{
  T=dim(X)[1]
  q=dim(Gam)[2]
  p=dim(Gam)[1]-q
  s=dim(Bet)[2]
  if(T<2) return(0)
  if(dim(X)[2]!=s) return(0)
  if(dim(Bet)[1]!=(p+q)) return(0)
  if(dim(G)[1]!=(p+q)) return(0)
  if(dim(G)[2]!=(p+q)) return(0)
  if(dim(h0)[1]!=q) return(0)
  if(dim(h0)[2]!=1) return(0)
  if(dim(X0)[1]!=s) return(0)
  if(dim(X0)[2]!=1) return(0)
  
  
  h=array(0,dim=c(T,q))# true states
  e=array(rnorm(T*(p+q)),dim=c(T,(p+q)))# innovations
  Z=array(0,dim=c(T,p))# data
  
  #Start
  
  temp=c(Gam%*%h0+Bet%*%X0+G%*%rnorm(p+q))
  h[1,]=temp[(p+1):(p+q)]
  for(i in 1:(T-1)){
    temp=c(Gam%*%h[i,]+Bet%*%X[i,]+G%*%e[i,])
    Z[i,]=temp[1:p]
    h[(i+1),]=temp[(p+1):(p+q)]
  }
  temp=c(Gam%*%h[T,]+Bet%*%X[T,]+G%*%e[T,])
  Z[T,]=temp[1:p]
  
  return(list(Z=Z,h=h))
}
