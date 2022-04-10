#include <Rcpp.h>
#include <datatableAPI.h>
using namespace Rcpp;

// from https://github.com/Rdatatable/data.table/issues/4643

// [[Rcpp::export]]
SEXP dtsubset(SEXP x, SEXP rows, SEXP cols) { return DT_subsetDT(x, rows, cols); }



