# All must be matrix data type!
#mX  : iT   by iu   matrix representing regressors
#Gam: ip+iq by iq   matrix representing parameters of states (H' F')'
#Bet: ip+iq by iu   matrix ... parameters of regressors (A' B')'
#mG  : ip+iq by ip+iq matrix ... parameters of errors (R' Q')'
#h0 : iq   by 1   matrix ... initial state
#mX0 : iu   by 1   matrix ... x at time step 0
#u  : iT   by ip+iq smoothed residuals obtained from Smoothing
#
#output
#Z  : iT   by ip   smoothed observations
#h  : iT   by iq   smoothed underlying process
SmoothProc<-function(mX,Gam,Bet,mG,h0,mX0,u)
{
  iT=dim(mX)[1]
  iq=dim(Gam)[2]
  ip=dim(Gam)[1]-iq
  iu=dim(Bet)[2]
  if(iT<2) return(0)
  if(dim(mX)[2]!=iu) return(0)
  if(dim(Bet)[1]!=(ip+iq)) return(0)
  if(dim(mG)[1]!=(ip+iq)) return(0)
  if(dim(mG)[2]!=(ip+iq)) return(0)
  if(dim(h0)[1]!=iq) return(0)
  if(dim(h0)[2]!=1) return(0)
  if(dim(mX0)[1]!=iu) return(0)
  if(dim(mX0)[2]!=1) return(0)
  
  
  h=array(0,dim=c(iT,iq))
  e=u# smoothed innovations
  Z=array(0,dim=c(iT,ip))# data
  
  #Start
  
  temp=c(Gam%*%h0+Bet%*%mX0+mG%*%rnorm(ip+iq))
  h[1,]=temp[(ip+1):(ip+iq)]
  for(i in 1:(iT-1)){
    temp=c(Gam%*%h[i,]+Bet%*%mX[i,]+mG%*%e[i,])
    Z[i,]=temp[1:ip]
    h[(i+1),]=temp[(ip+1):(ip+iq)]
  }
  temp=c(Gam%*%h[iT,]+Bet%*%mX[iT,]+mG%*%e[iT,])
  Z[iT,]=temp[1:ip]
  
  return(list(Z=Z,h=h))
}
