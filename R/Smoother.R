Smoother<-function(Z,X,Gam,Bet,G,h0,P0,X0)
# Need LinKF function preloaded
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
#KF : list       output from LinKF
#u  : T   by p+q E[eps |Z X]
#u0 : p+q by 1
#r  : T   by q
#r0 : q   by 1
#e  : T   by p
#D  : Tpp array
#N  : Tqq array
#N0 : q   by q
#C  : T(p+q)(p+q) array
#C0 : p+q  by p+q
#CN : T(p+q)q array
#L  : Tqq array
#L0 : q    by q
#M  : Tq(p+q) array
#M0 : q    by p+q
{

KF=LinKF(Z,X,Gam,Bet,G,h0,P0,X0)
if(class(KF)=="numeric") return(0)

T=dim(Z)[1]
p=dim(Z)[2]
q=dim(Gam)[2]
s=dim(Bet)[2]
F=Gam[(p+1):(p+q),]
dim(F)=c(q,q)
B=Bet[(p+1):(p+q),]
dim(B)=c(q,s)
Q=G[(p+1):(p+q),]
dim(Q)=c(q,(p+q))
QQ=Q%*%t(Q)
dim(QQ)=c(q,q)
H=Gam[1:p,]
dim(H)=c(p,q)
A=Bet[1:p,]
dim(A)=c(p,s)
R=G[1:p,]
dim(R)=c(p,(p+q))
RR=R%*%t(R)
dim(RR)=c(p,p)

e=array(0,dim=c(T,p))
r=array(0,dim=c(T,q))
u=array(0,dim=c(T,(p+q)))
D=array(0,dim=c(T,p,p))
N=array(0,dim=c(T,q,q))
C=array(0,dim=c(T,(p+q),(p+q)))
CN=array(0,dim=c(T,(p+q),q))
L=array(0,dim=c(T,q,q))
M=array(0,dim=c(T,q,(p+q)))

tool=cbind(rbind(H,F),rbind(R,Q))
for(i in T:2){
fk=F%*%KF$K[i,,]
dim(fk)=c(q,p)
L[i,,]=F-fk%*%H
M[i,,]=Q-fk%*%R
SS=solve(KF$S[i,,])
e[i,]=SS%*%KF$v[i,]-t(fk)%*%r[i,]
temp=t(tool)%*%c(e[i,],r[i,])
r[(i-1),]=temp[1:q,]
u[i,]=temp[(q+1):(p+2*q),]
D[i,,]=SS+t(fk)%*%N[i,,]%*%fk
temp=-N[i,,]%*%fk
dim(temp)=c(q,p)
temp=cbind(rbind(D[i,,],temp),rbind(t(temp),N[i,,]))
temp=t(tool)%*%temp%*%tool
N[(i-1),,]=temp[1:q,1:q]
C[i,,]=temp[(q+1):(p+2*q),(q+1):(p+2*q)]
CN[i,,]=temp[(q+1):(p+2*q),1:q]
}

fk=F%*%KF$K[1,,]
dim(fk)=c(q,p)
M[1,,]=Q-fk%*%R
L[1,,]=F-fk%*%H
SS=solve(KF$S[1,,])
e[1,]=SS%*%KF$v[1,]-t(fk)%*%r[1,]
temp=t(tool)%*%c(e[1,],r[1,])
r0=temp[1:q,]
u[1,]=temp[(q+1):(p+2*q),]
D[1,,]=SS+t(fk)%*%N[1,,]%*%fk
temp=-N[1,,]%*%fk
dim(temp)=c(q,p)
temp=cbind(rbind(D[1,,],temp),rbind(t(temp),N[1,,]))
temp=t(tool)%*%temp%*%tool
N0=temp[1:q,1:q]
C[1,,]=temp[(q+1):(p+2*q),(q+1):(p+2*q)]
CN[1,,]=temp[(q+1):(p+2*q),1:q]

u0=t(Q)%*%r0
C0=t(Q)%*%N0%*%Q
M0=Q
L0=F

return(list(KF=KF,u=u,u0=u0,r=r,r0=r0,e=e,D=D,N=N,N0=N0,C=C,C0=C0,CN=CN,L=L,L0=L0,M=M,M0=M0))
}
