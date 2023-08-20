# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
List mcalc_sim(const int iT, const int ip, const int iq, const arma::mat& mX, const arma::mat& me, const arma::mat& Gam, const arma::mat& Bet, const arma::mat& mG, const arma::vec& h0, const arma::vec& mX0) {
  // me must be (ip+iq) * (iT+1)
  arma::mat mh(iq, iT), mZ(ip, iT);

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


// [[Rcpp::export]]
List mcalc_kf(const int iT, const int ip, const int iq, const arma::mat& mZ, const arma::mat& mX, const arma::mat& mF, const arma::mat& mB, const arma::mat& mH, const arma::mat& mA, const arma::mat& RR, const arma::mat& QQ, const arma::vec& h0, const arma::vec& mX0, const arma::mat& P0) {
  // mZ must be ip * iT
  // mX must be iq * iT
  arma::mat h1(iq, iT), h2(iq, iT), mv(ip,iT);
  arma::cube P1(iq,iq,iT), P2(iq,iq,iT), aS(ip,ip,iT), aK(iq,ip,iT);

  // Predict
  h1.col(0) = mF * h0 + mB * mX0;
  P1.slice(0) = mF * P0 * mF.t() + QQ;

  //Update
  mv.col(0) = mZ.col(0) - mH * h1.col(0) - mA * mX.col(0);
  aS.slice(0) = mH * P1.slice(0) * mH.t() + RR;
  aK.slice(0) = P1.slice(0) * mH.t() * inv_sympd(aS.slice(0));
  h2.col(0) = h1.col(0) + aK.slice(0) * mv.col(0);
  P2.slice(0) = (eye(iq,iq) - aK.slice(0) * mH ) * P1.slice(0);

  for(int iter=1; iter < iT; iter++){
    // Predict
    h1.col(iter) = mF * h2.col(iter-1) + mB * mX.col(iter-1);
    //h1[i,]=F%*%h2[(i-1),]+B%*%mX[(i-1),]
    P1.slice(iter) = mF * P2.slice(iter-1) * mF.t() + QQ;
    //P1[i,,]=F%*%P2[(i-1),,]%*%t(F)+QQ

    //Update
    mv.col(iter) = mZ.col(iter) - mH * h1.col(iter) - mA * mX.col(iter);
    //v[i,]=mZ[i,]-H%*%h1[i,]-A%*%mX[i,]
    aS.slice(iter) = mH * P1.slice(iter) * mH.t() + RR;
    //S[i,,]=H%*%P1[i,,]%*%t(H)+RR
    aK.slice(iter) = P1.slice(iter) * mH.t() * inv_sympd(aS.slice(iter));
    //K[i,,]=P1[i,,]%*%t(H)%*%solve(S[i,,])
    h2.col(iter) = h1.col(iter) + aK.slice(iter) * mv.col(iter);
    //h2[i,]=h1[i,]+c(K[i,,]%*%v[i,])
    P2.slice(iter) = (eye(iq,iq) - aK.slice(iter) * mH ) * P1.slice(iter);
    //P2[i,,]=(diag(1,iq)-K[i,,]%*%H)%*%P1[i,,]
  }

  return List::create(
    Named("h1") = h1.t(),
    Named("h2") = h2.t(),
    Named("P1") = P1,
    Named("P2") = P2,
    Named("v") = mv.t(),
    Named("S") = aS,
    Named("K") = aK );
}