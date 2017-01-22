#include <Rcpp.h>
#include "MethylPatternData.h"
using namespace Rcpp;


//' Run Methyl Pattern Mixture Scanner (MPMS)
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




