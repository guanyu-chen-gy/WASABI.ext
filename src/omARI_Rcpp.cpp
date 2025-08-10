#include <Rcpp.h>
using namespace Rcpp;

//' This function computes the number of ways to choose 2 items from x items.
//'
//' @param x An integer, the number of items.
//' @return The number of ways to choose 2 items from x items, or 0 if x is less than or equal to 1.
// [[Rcpp::export]]
int comb2 (int x) 
{ return x > 1 ? x * (x - 1) / 2.0 : 0.0; }

//' Compute "one minus adjusted Rand index" between two partitions
//'
//' @param c1 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
//' @param c2 An integer vector of cluster labels for $n$ items. Should be 0-indexed.
//' @param K1 An integer, specifying the number of unique clusters in \code{cl1}.
//' @param K2 An integer, specifying the number of unique clusters in \code{cl2}.
//' @return The VI distance between \code{cl1} and \code{cl2}.
// [[Rcpp::export]]
double omARI_compute_Rcpp(NumericVector c1, NumericVector c2, int K1, int K2){
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
    for(int k2 = 0; k2 < K2; k2++){
       f += comb2(ctab(k1,k2));
    }
  }  
  
  double g = 0.0;
  for(int k1 = 0; k1 < K1; k1++){
    g += comb2(r[k1]);
  }

  double h = 0.0;
  for(int k2 = 0; k2 < K2; k2++){
    h += comb2(c[k2]);
  }

  double p = comb2(n); 
  double AR = (f - g * h / p) / ((g + h) * 0.5 - g * h / p); 
  return 1 - AR;
}

//' omARI distance between two groups of partitions
//'
//' @param cls1 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
//' @param cls2 A matrix of $n$ columns, where each row contains the 0-indexed cluster labels for n items. 
//' @param K1s An integer vector, specifying the number of clusters in each row of \code{cls1}. Must have length equal to the number of rows in \code{cls1}.
//' @param K2s An integer vector, specifying the number of clusters in each row of \code{cls2}. Must have length equal to the number of rows in \code{cls2}.
//' @return A matrix containing the VI distance between each row of \code{cls1} and \code{cls2}.
// [[Rcpp::export]]
NumericMatrix omARI_Rcpp(NumericMatrix cls1,NumericMatrix cls2,NumericVector K1s,NumericVector K2s) {

  if(cls1.ncol() != cls2.ncol()){
    stop("number of data points must be the same");
  }
  
  int M1=cls1.nrow();
  int M2=cls2.nrow();
    
  NumericMatrix output(M1,M2);
  
  for (int i = 0; i < M1; i++) {
    for(int j = 0; j < M2; j++){
      output(i,j) = omARI_compute_Rcpp(cls1(i,_),cls2(j,_),K1s[i],K2s[j]);
    }
  }
  
  return output;
}