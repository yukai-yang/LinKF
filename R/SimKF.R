#' Simulate from a linear state space model.
#'
#' \code{SimKF} implements the simulation from a linear state space model.
#'
#' The function takes several arguments that specify the linear state space model
#' and returns a list consisting of the simulated latent process and the observations.
#'
#' The corresponding linear state space model is as follows:
#'
#' \deqn{(z_t', h_{t+1}') = Gam h_t + Bet x_t + G e_t,}
#' where \eqn{Gam = (H', F')'}, \eqn{Bet = (A', B')'}, and \eqn{G = (R', Q')'}
#'
#' The seed can be set before running the function.
#'
#' @param mX  a T by u matrix representing the explanatory variable $X$
#' @param Gam a p+q by q matrix representing the coefficients $(H' F')'$ of the latent states
#' @param Bet a p+q by u matrix representing the coefficients $(A' B')'$ of $X$
#' @param mG a p+q by p+q matrix representing the coefficients of errors $(R' Q')'$
#' @param h0 a q by 1 matrix representing the initial state
#' @param mX0 a u by 1  matrix representing $x_0$ at time 0
#'
#' @returns A list consisting of the following components:
#'
#' \item{Z}{T by p  the observations}
#' \item{h}{T by q the latent process of the states}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' @export
SimKF <- function(mX,Gam,Bet,mG,h0,mX0)
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
  #e=matrix(rnorm(iT*(ip+iq)),iT,(ip+iq))# innovations
  e=matrix(rnorm((ip+iq)*(iT+1)),ip+iq,iT+1)# innovations
  mZ=matrix(0,iT,ip)# data

  #Start

  temp=c(Gam%*%h0+Bet%*%mX0+mG%*%e[,1])
  h[1,]=temp[(ip+1):(ip+iq)]

  for(i in 1:(iT-1)){
    temp=c(Gam%*%h[i,]+Bet%*%mX[i,]+mG%*%e[,i+1])
    mZ[i,]=temp[1:ip]
    h[(i+1),]=temp[(ip+1):(ip+iq)]
  }

  temp=c(Gam%*%h[iT,]+Bet%*%mX[iT,]+mG%*%e[,iT+1])
  mZ[iT,]=temp[1:ip]

  return(list(Z=mZ,h=h))
}
