#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector test(int n){
  NumericVector v1 , v2;
  v1 = seq(1, n-10);
  v2 = seq(n - 9, n);

  //NumericVector X = NumericVector::create(v1, v2);
  return v1;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
test(100)
*/
