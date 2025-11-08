#include <Rcpp.h>

using namespace Rcpp;

//' Compute Generalized Variation of Information (VI) distance between two partitions
//'
//' @param c1 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
//' @param c2 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
//' @param K1 An integer, specifying the number of unique clusters in \code{cl1}.
//' @param K2 An integer, specifying the number of unique clusters in \code{cl2}.
//' @param a A parameter used in generalized VI, takes value between 0 and 2 and 1 by default.
//' @return The Generalized VI distance between \code{cl1} and \code{cl2}.
// [[Rcpp::export]]
double VI_compute_Rcpp(NumericVector c1, NumericVector c2, int K1, int K2, double a = 1.0) {
  int n=c1.length();
  
  NumericMatrix ctab(K1,K2);
  NumericVector r(K1);
  NumericVector c(K2);
  for(int i = 0; i < n; i++){
    ctab(c1[i],c2[i])++;
  }
  for(int k1 = 0; k1 < K1; k1++){
    for(int k2 = 0; k2 < K2; k2++){
      r[k1] += ctab(k1,k2);
      c[k2] += ctab(k1,k2);
    }
  }
  
  double f = 0.0;
  for(int k1 = 0; k1 < K1; k1++){
    if(r[k1] > 0){
      f += a * r[k1]/n*(log2(r[k1])-log2(n));
    }
  }
  for(int k2 = 0; k2 < K2; k2++){
    if(c[k2]>0){
      f += ( 2 - a ) *c[k2]/n*(log2(c[k2])-log2(n));
    }
  }
  for(int k1 = 0; k1 < K1; k1++){
    for(int k2 = 0; k2 < K2; k2++){
      if(ctab(k1,k2)>0){
        f += -2*ctab(k1,k2)/n*(log2(ctab(k1,k2))-log2(n));
      }
    }
  }
  return(f);
}

//' Compute Generalized Variation of Information (VI) distance between two groups of partitions
//'
//' @param cls1 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
//' @param cls2 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
//' @param K1s An integer vector, specifying the number of clusters in each row of \code{cls1}. Must have length equal to the number of rows in \code{cls1}.
//' @param K2s An integer vector, specifying the number of clusters in each row of \code{cls2}. Must have length equal to the number of rows in \code{cls2}.
//' @param a A parameter used in generalized VI, takes value between 0 and 2 and 1 by default.
//' @return A matrix containing the Generalized VI distance between each row of \code{cls1} and \code{cls2}.
// [[Rcpp::export]]
NumericMatrix VI_Rcpp(NumericMatrix cls1,NumericMatrix cls2,NumericVector K1s,NumericVector K2s, double a = 1.0) {

  if(cls1.ncol() != cls2.ncol()){
    stop("number of data points must be the same");
  }
  
  int M1=cls1.nrow();
  int M2=cls2.nrow();
    
  NumericMatrix output(M1,M2);
  
  for (int i = 0; i < M1; i++) {
    for(int j = 0; j < M2; j++){
      output(i,j) = VI_compute_Rcpp(cls1(i,_),cls2(j,_),K1s[i],K2s[j], a);
    }
  }
  
  return output;
}



