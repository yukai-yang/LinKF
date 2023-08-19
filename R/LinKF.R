#' Linear Kalman filtering
#'
#' \code{LinKF} implements the the linear Kalman filtering algorithm.
#'
#' The function takes several arguments that specify the linear state space model as well as the data
#' and returns a list consisting of the estimates from the Kalman filtering algorithm.
#'
#' The corresponding linear state space model is as follows:
#'
#' \deqn{(z_t', h_{t+1}') = Gam h_t + Bet x_t + G e_t,}
#' where \eqn{Gam = (H', F')'}, \eqn{Bet = (A', B')'}, and \eqn{G = (R', Q')'}
#'
#' @param mZ a T by p matrix representing observations
#' @param mX a T by u matrix representing explanatory variables
#' @param Gam a p+q by q matrix representing the coefficients $(H' F')'$ of the latent states
#' @param Bet a p+q by u matrix representing the coefficients $(A' B')'$ of $X$
#' @param mG a p+q by p+q matrix representing the coefficients of errors $(R' Q')'$
#' @param h0 a q by 1 matrix representing the initial state
#' @param P0 a q by q matrix representing the initial state's covariance matrix
#' @param mX0 a u by 1  matrix representing $x_0$ at time 0
#'
#' @returns A list consisting of the following components:
#'
#' \item{h1}{T by q the estimated conditional expectations of the latent states given t-1 (a priori)}
#' \item{h2}{T by q the estimated conditional expectations of the latent states given t (a posteriori)}
#' \item{P1}{T by q by q array of covariances of the a priori states}
#' \item{P2}{T by q by q array of covariances of the a posteriori states}
#' \item{v}{T by p residuals in the measurement eq.}
#' \item{S}{T by p by p array of covariances of z}
#' \item{K}{T by q by p array of Kalman gains}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' @export
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
