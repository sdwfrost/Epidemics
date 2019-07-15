#include <Rcpp.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::plugins(unwindProtect)]]
// [[Rcpp::export]]
List event_epidemicsC(NumericVector individual_inf_rate, double gamma, int Y, NumericVector U){


  // Return Variables
  int event;
  int ID_index;
  //Check whether U has been provided
  //If not draw Unif(0,1)
  LogicalVector U_not_provided = is_nan(U);
  //std::cout << U_not_provided << "\n";
  if(U_not_provided[0]){
    U = runif(1, 0, 1);
  }
  //std::cout << U << "\n";
  int X = individual_inf_rate.size();

  NumericVector rates(X + Y);
  NumericVector removal_rates = rep(gamma, Y);

  //rates = c(individual_inf_rate, removal_rates);
  for(int i = 0; i < X; ++i){
    rates[i] = individual_inf_rate[i];
  }
  for(int j = 0; j < Y; ++j){
    rates[X + j] = removal_rates[j];
  }
  //std::cout << rates << "\n";
  NumericVector cumsum_rates = cumsum(rates);
  //std::cout << cumsum_rates[0] << "\n";
  double total_rate = sum(rates);

  // Finding Index of individual event occurred to
  for(int i = 0; i < X + Y; ++i){
    //std::cout << i << "   " << cumsum_rates[i]/total_rate << "\n";
    if(cumsum_rates[i]/total_rate > U[0]){
      ID_index = i + 1;
      break;
    }
  }

  // Adjusting position depending on the current of
  // Individual
  if(ID_index <= X){
    event = 0;
  } else{
    event = 1;
    ID_index = ID_index - X;
  }
  return List::create(Named("event") = event, Named("ID_index") = ID_index);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

# X = 95
# Y = 5
# gamma = 0.5
# beta = 0.001
# individual_inf_rate = rep(0.001*Y, X)
# U = runif(1)
# # example = event_epidemicsC(individual_inf_rate, gamma, Y, if(missing(U)){U = NaN} else(U))
# # # #rbenchmark::benchmark(event_epidemicsC(individual_inf_rate, gamma, Y, U),
# # # #                      event.epidemics(individual_inf_rate, gamma, Y, U), replications = 7e6)
#  bench::mark(event_epidemicsC(individual_inf_rate, gamma, Y, U),
#              event.epidemics(individual_inf_rate, gamma, Y, U))
*/
