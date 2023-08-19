# All must be matrix data type!
#
# the model:
# (z_t', h_{t+1}') = Gam h_t + Bet x_t + G e_t
# where Gam = (H', F')', Bet = (A', B')', G = (R', Q')'
#
#mX  : iT   by iu   matrix representing regressors
#Gam : ip+iq by iq   matrix representing parameters of states (H' F')'
#Bet : ip+iq by iu   matrix ... parameters of regressors (A' B')'
#mG  : ip+iq by ip+iq matrix ... parameters of errors (R' Q')'
#h0  : iq   by 1   matrix ... initial state
#mX0 : iu   by 1   matrix ... x at time step 0
#
#output
#Z  : iT   by ip   observations
#h  : iT   by iq   underlying process
SimKF<-function(mX,Gam,Bet,mG,h0,mX0)
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
  
  
  h=matrix(0,iT,iq)# true states
  e=matrix(rnorm(iT*(ip+iq)),iT,(ip+iq))# innovations
  mZ=matrix(0,iT,ip)# data
  
  #Start
  
  temp=c(Gam%*%h0+Bet%*%mX0+mG%*%rnorm(ip+iq))
  h[1,]=temp[(ip+1):(ip+iq)]

  for(i in 1:(iT-1)){
    temp=c(Gam%*%h[i,]+Bet%*%mX[i,]+mG%*%e[i,])
    mZ[i,]=temp[1:ip]
    h[(i+1),]=temp[(ip+1):(ip+iq)]
  }

  temp=c(Gam%*%h[iT,]+Bet%*%mX[iT,]+mG%*%e[iT,])
  mZ[iT,]=temp[1:ip]
  
  return(list(Z=mZ,h=h))
}
