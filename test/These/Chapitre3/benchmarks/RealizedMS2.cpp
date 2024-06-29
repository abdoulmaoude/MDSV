#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,const bool &LEVIER, const int &K){
  
  int l=0;
  if(LEVIER) l = 2;
  Eigen::MatrixXd P(K,K);
  Eigen::VectorXd Vol(K);
  Eigen::VectorXd probapi(K);
  double temp=0;

  for(int k=0; k<K; k++){
    Vol(k)=exp(para_tilde(l+k));   //volatility
  }

  for(int i=0; i<K; i++){
    for(int j=0; j<(K-1); j++){
    temp=temp+exp(para_tilde(l+K+i*(K-1)+j));
    }
    for(int j=0; j<(K-1); j++){
      P(i,j)=exp(para_tilde(l+K+i*(K-1)+j))/(1+temp);      //pij
    }
    P(i,(K-1))=1/(1+temp);
    temp=0;
  }

  Eigen::MatrixXd A = (Eigen::MatrixXd::Identity(K,K)-P+Eigen::MatrixXd::Ones(K,K)).transpose();
  probapi=A.colPivHouseholderQr().solve(Eigen::VectorXd::Ones(K));
  
    Rcpp::List output =
      Rcpp::List::create(
        Rcpp::Named("P") = P,
        Rcpp::Named("Vol") = Vol,
        Rcpp::Named("probapi") = probapi
      );
    
    if(LEVIER){
      Eigen::VectorXd para(2);
      para(0)=exp(para_tilde(0)); //l
      para(1)=1/(1+exp(para_tilde(1))); //theta
      
      output =
        Rcpp::List::create(
          Rcpp::Named("para") = para,
          Rcpp::Named("P") = P,
          Rcpp::Named("Vol") = Vol,
          Rcpp::Named("probapi") = probapi
        );
    }
    return output;
}


// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para, const Eigen::Map<Eigen::VectorXd> &Vol,
                        const Eigen::Map<Eigen::MatrixXd> &P,const bool &LEVIER, const int &K){
  int l=0;
  if(LEVIER) l=2;
  Eigen::VectorXd para_tilde(l+K*K);

  if(LEVIER){
    para_tilde[0]=log(para[0]); //l
    para_tilde[1]=log((1/para[1])-1); //theta
  }

  for(int k=0; k<K; k++){
    para_tilde[l+k]=log(Vol(k));            //volatility
  }

  for(int i=0; i<K; i++){
    for(int j=0; j<(K-1); j++){
      para_tilde[l+K+i*(K-1)+j]=log((P(i,j))/(P(i,K-1)));      //pij
    }
  }
  return para_tilde;
}

// // [[Rcpp::export]]
// NumericVector proba_pi(const double &omega, const int &N_D){
//   NumericVector probaPi=dbinom(seq(0,N_D-1),N_D-1,omega);
//   return probaPi;
// }

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){

  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[0]*pow(para[1],i);
  for(int t=Nl;t<ech.size();t++){
    levier=1;
    for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
    Levier[t]=levier;
  }
  levier=1;
  for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("Levier") = Levier,
      Rcpp::Named("levier") = levier
    );
  return output;
}
// [[Rcpp::export]]
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
double logLik(const Eigen::MatrixXd &ech,
                    const Eigen::Map<Eigen::VectorXd> &para_tilde,
                    const bool &LEVIER,
                    const int &K,
                    const int &Nl=70){
  
  Rcpp::List para=workNat(para_tilde,LEVIER,K);
  int n=ech.rows();
  Eigen::VectorXd sigma = para["Vol"];
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=para["probapi"];
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  double x=ech(0,1);
  double y=ech(0,0);
  Eigen::MatrixXd matP = para["P"];
  Eigen::VectorXd Res;
  
  Eigen::VectorXd temp;
  aj=r_dens(x,sigma).array();
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      x=ech(i,1);
      y=ech(i,0);
      Res=r_dens(x,sigma).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }else{
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Eigen::VectorXd opara = para["para"];
    Rcpp::List L = levierVolatility(ech1,Nl,opara);
    NumericVector Levier = L["Levier"];
    Eigen::VectorXd sigma1;
    
    for(int i(1); i<n;i++){
      x = ech(i,1);
      y = ech(i,0);
      sigma1 = sigma*Levier[i];
      Res = r_dens(x,sigma1).array();
      aj = (((w.transpose())*matP).array())*(Res.transpose().array());
      a = aj.sum();
      lik += log(a);
      w = aj/a;
    }
  }
  return (-lik);
}


/*** R



*/
