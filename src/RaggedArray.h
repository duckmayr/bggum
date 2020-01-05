#ifndef BGGUM_RAGGEDARRAY_H
#define BGGUM_RAGGEDARRAY_H

#include <Rcpp.h>

/* tau storage
 * We do not need to take advantage of the potential for data heterogeneity
 * that Rcpp::List was designed for; so, even though on the R side, tau will
 * be stored in a list, on the C++ side we will use ragged arrays.
 */

typedef std::vector<Rcpp::NumericVector> RaggedArray;
typedef std::vector<std::vector<Rcpp::NumericVector> > NestedRaggedArray;

inline RaggedArray clone(RaggedArray x) {
    int N = x.size();
    RaggedArray result(N);
    for ( int i = 0; i < N; ++i ) {
        result[i] = Rcpp::clone(x[i]);
    }
    return result;
}

#endif
