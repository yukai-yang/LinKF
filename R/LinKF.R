#' Linear Kalman Filter
#'
#'
# All must be matrix data type!
#Z  : T   by p   matrix representing observations
#X  : T   by s   matrix representing regressors
#Gam: p+q by q   matrix representing parameters of states (H' F')'
#Bet: p+q by s   matrix ... parameters of regressors (A' B')'
#G  : p+q by p+q matrix ... parameters of errors (R' Q')'
#h0 : q   by 1   matrix ... accurate initial state
#P0 : q   by q   matrix ... accurate initial state's cov
#X0 : s   by 1   matrix ... x at time step 0
#by Y. K. Yang
#output
#h1 : T   by q   expectation of h a priori
#h2 : T   by q   expectation of h a posteriori
#P1 : Tqq array  variance of h a priori
#P2 : Tqq array  variance of h a posteriori
#v  : T   by p   residuals in measurement eq.
#S  : Tpp array  covariance of z
#K  : Tqp array  Kalman gain
LinKF<-function(Z,X,Gam,Bet,G,h0,P0,X0)
{
  
  #Initial
  T=dim(Z)[1]
  p=dim(Z)[2]
  q=dim(Gam)[2]
  s=dim(Bet)[2]
  if(T<2) return(0)
  if(dim(Gam)[1]!=(p+q)) return(0)
  if(dim(Bet)[1]!=(p+q)) return(0)
  if(dim(G)[1]!=(p+q)) return(0)
  if(dim(G)[2]!=(p+q)) return(0)
  if(dim(X)[1]!=T) return(0)
  if(dim(X)[2]!=s) return(0)
  if(dim(h0)[1]!=q) return(0)
  if(dim(h0)[2]!=1) return(0)
  if(dim(P0)[1]!=q) return(0)
  if(dim(P0)[2]!=q) return(0)
  if(dim(X0)[1]!=s) return(0)
  if(dim(X0)[2]!=1) return(0)
  
  h1=array(0,dim=c(T,q))# h_t|t-1
  P1=array(0,dim=c(T,q,q))# P_t|t-1
  
  v=array(0,dim=c(T,p))# residual in z eq.
  S=array(0,dim=c(T,p,p))# cov of residual
  K=array(0,dim=c(T,q,p))# Kalman gain
  h2=array(0,dim=c(T,q))# h_t|t
  P2=array(0,dim=c(T,q,q))# P_t|t
  
  F=Gam[(p+1):(p+q),]
  B=Bet[(p+1):(p+q),]
  temp=G[(p+1):(p+q),]
  dim(temp)=c(q,(p+q))
  QQ=temp%*%t(temp)
  H=Gam[1:p,]
  A=Bet[1:p,]
  temp=G[1:p,]
  dim(temp)=c(p,(p+q))
  RR=temp%*%t(temp)
  
  #Start
  
  #Predict
  h1[1,]=F%*%h0+B%*%X0
  P1[1,,]=F%*%P0%*%t(F)+QQ
  #Update
  v[1,]=Z[1,]-H%*%h1[1,]-A%*%X[1,]
  S[1,,]=H%*%P1[1,,]%*%t(H)+RR
  K[1,,]=P1[1,,]%*%t(H)%*%solve(S[1,,])
  h2[1,]=h1[1,]+c(K[1,,]%*%v[1,])
  P2[1,,]=(diag(1,q)-K[1,,]%*%H)%*%P1[1,,]
  
  for(i in 2:T){
    #Predict
    h1[i,]=F%*%h2[(i-1),]+B%*%X[(i-1),]
    P1[i,,]=F%*%P2[(i-1),,]%*%t(F)+QQ
    #Update
    v[i,]=Z[i,]-H%*%h1[i,]-A%*%X[i,]
    S[i,,]=H%*%P1[i,,]%*%t(H)+RR
    K[i,,]=P1[i,,]%*%t(H)%*%solve(S[i,,])
    h2[i,]=h1[i,]+c(K[i,,]%*%v[i,])
    P2[i,,]=(diag(1,q)-K[i,,]%*%H)%*%P1[i,,]
  }
  
  return(list(h1=h1,h2=h2,P1=P1,P2=P2,v=v,S=S,K=K))
}
