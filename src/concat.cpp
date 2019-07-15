#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector concat(NumericVector x, NumericVector y) {
  int n_x = x.size();
  int n_y = y.size();
  NumericVector c_xy(n_x + n_y);

  for(int i = 0; i < n_x; i++){
    c_xy[i] = x[i];
  }

  for(int j = 0; j < n_y; j++){
    c_xy[n_x + j] = y[j];
  }

  return c_xy;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x = runif(5e4)
y = runif(5e4)
bench::mark(concat(x, y), c(x,y))
*/
