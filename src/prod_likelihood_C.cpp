#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
using namespace Rcpp;

// [[Rcpp::export]]
double prod_part_infC(NumericVector inf_times, NumericVector rem_times, NumericMatrix B,
                       bool log = true){
  int n = inf_times.size();
  double Inf = std::numeric_limits<double>::infinity();
  // Which individuals are infected
  LogicalVector is_infected(n);

  for(int i = 0; i < n; ++i){
    if(inf_times[i] < Inf){
      is_infected[i] = true;
    } else{
      is_infected[i] = false;
    }
  }




  return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
