#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector edab(int win, NumericVector x, NumericVector y){
  double r = 0;
  double g = 0;
  for(int a = 1; a <= win-1; a++){
    for(int b = a+1; b <= win; b++){
      r+= (std::max((y(b-1)-y(a-1))-(x(b-1)-x(a-1)), 0.0)) * (a * (win+1-b));
      g+= (y(b-1)-y(a-1))* (a * (win+1-b));
    }
  }

  return NumericVector::create(Named("r",r), Named("g",g));
}