#include <Rcpp.h>
#include "MethylPatternData.h"
using namespace Rcpp;

//' Run mpms
//'
//' @param df data frame
//' @export
// [[Rcpp::export]]
DataFrame mpms(const DataFrame& df){
    Rcout << "mpms" << endl;
    MethylPatternData d(df);
    Rcout << d.avg_m() << endl;
    Rcout << d.get_na() << endl;
    Rcout << d.get_pattern_length() << endl;
    Rcout << d.get_number_of_samples() << endl;
    // d.print_samples(cerr);
    d.print_samples(Rcout);
    // Rcout << d << endl;
    return(df);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
