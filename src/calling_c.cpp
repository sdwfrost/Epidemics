#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector calling_c(NumericVector x, NumericVector y) {
  // calling c()
  Function c("c");

  NumericVector c_xy = c(x,y);
  return c_xy;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
calling_c(1:10, 11:20)
*/
