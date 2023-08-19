#' Linear Kalman Filter
#'
# the model:
# (z_t', h_{t+1}') = Gam h_t + Bet x_t + G e_t
# where Gam = (H', F')', Bet = (A', B')', G = (R', Q')'
#'
# All must be matrix data type!
#mZ  : iT   by ip   matrix representing observations
#mX  : iT   by iu   matrix representing regressors
#Gam: ip+iq by iq   matrix representing parameters of states (H' F')'
#Bet: ip+iq by iu   matrix ... parameters of regressors (A' B')'
#mG  : ip+iq by ip+iq matrix ... parameters of errors (R' Q')'
#h0 : iq   by 1   matrix ... accurate initial state
#P0 : iq   by iq   matrix ... accurate initial state's cov
#mX0 : iu   by 1   matrix ... x at time step 0
#
#output
#h1 : iT   by iq   expectation of h a priori
#h2 : iT   by iq   expectation of h a posteriori
#P1 : Tqq array  variance of h a priori
#P2 : Tqq array  variance of h a posteriori
#v  : iT   by ip   residuals in measurement eq.
#S  : Tpp array  covariance of z
#K  : Tqp array  Kalman gain
LinKF<-function(mZ,mX,Gam,Bet,mG,h0,P0,mX0)
{
  
  #Initial
  iT=dim(mZ)[1]
  ip=dim(mZ)[2]
  iq=dim(Gam)[2]
  iu=dim(Bet)[2]
  if(iT<2) return(0)
  if(dim(Gam)[1]!=(ip+iq)) return(0)
  if(dim(Bet)[1]!=(ip+iq)) return(0)
  if(dim(mG)[1]!=(ip+iq)) return(0)
  if(dim(mG)[2]!=(ip+iq)) return(0)
  if(dim(mX)[1]!=iT) return(0)
  if(dim(mX)[2]!=iu) return(0)
  if(dim(h0)[1]!=iq) return(0)
  if(dim(h0)[2]!=1) return(0)
  if(dim(P0)[1]!=iq) return(0)
  if(dim(P0)[2]!=iq) return(0)
  if(dim(mX0)[1]!=iu) return(0)
  if(dim(mX0)[2]!=1) return(0)
  
  h1=array(0,dim=c(iT,iq))# h_t|t-1
  P1=array(0,dim=c(iT,iq,iq))# P_t|t-1
  
  v=array(0,dim=c(iT,ip))# residual in z eq.
  S=array(0,dim=c(iT,ip,ip))# cov of residual
  K=array(0,dim=c(iT,iq,ip))# Kalman gain
  h2=array(0,dim=c(iT,iq))# h_t|t
  P2=array(0,dim=c(iT,iq,iq))# P_t|t
  
  F=Gam[(ip+1):(ip+iq),]
  B=Bet[(ip+1):(ip+iq),]
  temp=mG[(ip+1):(ip+iq),]
  dim(temp)=c(iq,(ip+iq))
  QQ=temp%*%t(temp)
  H=Gam[1:ip,]
  A=Bet[1:ip,]
  temp=mG[1:ip,]
  dim(temp)=c(ip,(ip+iq))
  RR=temp%*%t(temp)
  
  #Start
  
  #Predict
  h1[1,]=F%*%h0+B%*%mX0
  P1[1,,]=F%*%P0%*%t(F)+QQ
  #Update
  v[1,]=mZ[1,]-H%*%h1[1,]-A%*%mX[1,]
  S[1,,]=H%*%P1[1,,]%*%t(H)+RR
  K[1,,]=P1[1,,]%*%t(H)%*%solve(S[1,,])
  h2[1,]=h1[1,]+c(K[1,,]%*%v[1,])
  P2[1,,]=(diag(1,iq)-K[1,,]%*%H)%*%P1[1,,]
  
  for(i in 2:iT){
    #Predict
    h1[i,]=F%*%h2[(i-1),]+B%*%mX[(i-1),]
    P1[i,,]=F%*%P2[(i-1),,]%*%t(F)+QQ
    #Update
    v[i,]=mZ[i,]-H%*%h1[i,]-A%*%mX[i,]
    S[i,,]=H%*%P1[i,,]%*%t(H)+RR
    K[i,,]=P1[i,,]%*%t(H)%*%solve(S[i,,])
    h2[i,]=h1[i,]+c(K[i,,]%*%v[i,])
    P2[i,,]=(diag(1,iq)-K[i,,]%*%H)%*%P1[i,,]
  }
  
  return(list(h1=h1,h2=h2,P1=P1,P2=P2,v=v,S=S,K=K))
}
