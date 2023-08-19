# Need LinKF function preloaded
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
#KF : list       output from LinKF
#u  : iT   by ip+iq E[eps |Z X]
#u0 : ip+iq by 1
#r  : iT   by iq
#r0 : iq   by 1
#e  : iT   by ip
#D  : Tpp array
#N  : Tqq array
#N0 : iq   by iq
#C  : T(ip+iq)(ip+iq) array
#C0 : ip+iq  by ip+iq
#CN : T(ip+iq)iq array
#L  : Tqq array
#L0 : iq    by iq
#M  : Tq(p+q) array
#M0 : iq    by ip+iq
Smoother<-function(mZ,mX,Gam,Bet,mG,h0,P0,mX0)
{
  
  KF=LinKF(mZ,mX,Gam,Bet,mG,h0,P0,mX0)
  if(class(KF)=="numeric") return(0)
  
  iT=dim(mZ)[1]
  ip=dim(mZ)[2]
  iq=dim(Gam)[2]
  iu=dim(Bet)[2]
  F=Gam[(ip+1):(ip+iq),]
  dim(F)=c(iq,iq)
  B=Bet[(ip+1):(ip+iq),]
  dim(B)=c(iq,iu)
  Q=mG[(ip+1):(ip+iq),]
  dim(Q)=c(iq,(ip+iq))
  QQ=Q%*%t(Q)
  dim(QQ)=c(iq,iq)
  H=Gam[1:ip,]
  dim(H)=c(ip,iq)
  A=Bet[1:ip,]
  dim(A)=c(ip,iu)
  R=mG[1:ip,]
  dim(R)=c(ip,(ip+iq))
  RR=R%*%t(R)
  dim(RR)=c(ip,ip)
  
  e=array(0,dim=c(iT,ip))
  r=array(0,dim=c(iT,iq))
  u=array(0,dim=c(iT,(ip+iq)))
  D=array(0,dim=c(iT,ip,ip))
  N=array(0,dim=c(iT,iq,iq))
  C=array(0,dim=c(iT,(ip+iq),(ip+iq)))
  CN=array(0,dim=c(iT,(ip+iq),iq))
  L=array(0,dim=c(iT,iq,iq))
  M=array(0,dim=c(iT,iq,(ip+iq)))
  
  tool=cbind(rbind(H,F),rbind(R,Q))
  for(i in iT:2){
    fk=F%*%KF$K[i,,]
    dim(fk)=c(iq,ip)
    L[i,,]=F-fk%*%H
    M[i,,]=Q-fk%*%R
    SS=solve(KF$S[i,,])
    e[i,]=SS%*%KF$v[i,]-t(fk)%*%r[i,]
    temp=t(tool)%*%c(e[i,],r[i,])
    r[(i-1),]=temp[1:iq,]
    u[i,]=temp[(iq+1):(ip+2*iq),]
    D[i,,]=SS+t(fk)%*%N[i,,]%*%fk
    temp=-N[i,,]%*%fk
    dim(temp)=c(iq,ip)
    temp=cbind(rbind(D[i,,],temp),rbind(t(temp),N[i,,]))
    temp=t(tool)%*%temp%*%tool
    N[(i-1),,]=temp[1:iq,1:iq]
    C[i,,]=temp[(iq+1):(ip+2*iq),(iq+1):(ip+2*iq)]
    CN[i,,]=temp[(iq+1):(ip+2*iq),1:iq]
  }
  
  fk=F%*%KF$K[1,,]
  dim(fk)=c(iq,ip)
  M[1,,]=Q-fk%*%R
  L[1,,]=F-fk%*%H
  SS=solve(KF$S[1,,])
  e[1,]=SS%*%KF$v[1,]-t(fk)%*%r[1,]
  temp=t(tool)%*%c(e[1,],r[1,])
  r0=temp[1:iq,]
  u[1,]=temp[(iq+1):(ip+2*iq),]
  D[1,,]=SS+t(fk)%*%N[1,,]%*%fk
  temp=-N[1,,]%*%fk
  dim(temp)=c(iq,ip)
  temp=cbind(rbind(D[1,,],temp),rbind(t(temp),N[1,,]))
  temp=t(tool)%*%temp%*%tool
  N0=temp[1:iq,1:iq]
  C[1,,]=temp[(iq+1):(ip+2*iq),(iq+1):(ip+2*iq)]
  CN[1,,]=temp[(iq+1):(ip+2*iq),1:iq]
  
  u0=t(Q)%*%r0
  C0=t(Q)%*%N0%*%Q
  M0=Q
  L0=F
  
  return(list(KF=KF,u=u,u0=u0,r=r,r0=r0,e=e,D=D,N=N,N0=N0,C=C,C0=C0,CN=CN,L=L,L0=L0,M=M,M0=M0))
}
