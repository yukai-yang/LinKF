# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
List mcalc(const int iT, const int ip, const int iq, const arma::mat& mX, const arma::mat& me, const arma::mat& Gam, const arma::mat& Bet, const arma::mat& mG, const arma::vec& h0, const arma::vec& mX0) {
  // me must be (ip+iq) * (iT+1)
  arma::mat mh, mZ;
  mh.set_size(iq, iT);
  mZ.set_size(ip, iT);

  arma::vec tmp = Gam*h0 + Bet*mX0 + mG*me.col(0);
  mh.col(0) = tmp.tail(iq); // tmp[(ip+1):(ip+iq)]

  for(int iter=0; iter < iT-1; iter++){
    tmp = Gam*mh.col(iter) + Bet*mX.col(iter) + mG*me.col(iter+1);
    //temp=c(Gam%*%h[i,]+Bet%*%mX[i,]+mG%*%e[i,])
    mZ.col(iter) = tmp.head(ip);
    //mZ[i,]=temp[1:ip]
    mh.col(iter+1) = tmp.tail(iq);
    //h[(i+1),]=temp[(ip+1):(ip+iq)]
  }

  tmp = Gam*mh.col(iT-1) + Bet*mX.col(iT-1) + mG*me.col(iT);
  //temp=c(Gam%*%h[iT,]+Bet%*%mX[iT,]+mG%*%e[iT,])
  mZ.col(iT-1) = tmp.head(ip);
  //mZ[iT,]=temp[1:ip]

  return List::create(
    Named("Z") = mZ.t(),
    Named("h") = mh.t() );
}
