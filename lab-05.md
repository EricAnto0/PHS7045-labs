Lab 05 - Rcpp
================

# Learning goals

- Use the different data types in Rcpp.
- Learn some fundamentals about C++ optimization.
- Practice your GitHub skills.

# Lab description

For this lab, we will create a function for propensity score matching.
The goal is simple: write out a C++ function with Rcpp and measure how
faster it is compared to the following R implementation:

``` r
ps_matchR <- function(x) {
  
  match_expected <- as.matrix(dist(x))
  diag(match_expected) <- .Machine$integer.max
  indices <- apply(match_expected, 1, which.min)
  
  list(
    match_id = as.integer(unname(indices)),
    match_x  = x[indices]
  )
  
}
```

## Question 1: Create a simple function

Use the following pseudo-code template to get started:

``` cpp
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List ps_match1(const NumericVector & x) {
     
    int n = x.size();
    integerVector values(n);

    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < x.size(); ++i) {
            if (...the closests so far...) {
                ...update the optimum...
            }
        }
        
    }

    return [a list like the R function]

}
```

``` rcpp
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List ps_match1(const NumericVector & x) {
  
  int n = x.size();
  IntegerVector indices(n);
  NumericVector values(n);
  
  for (int i = 0; i < n; ++i) {
    int best_n = 0;
    double best_dist = std::numeric_limits< double >::max();
    for (int j = 0; j < n; ++j) {
      if(i==j)
        continue;
      double tmp_dist = abs(x[i]-x[j]);
      
      if (tmp_dist < best_dist) {
          best_dist = tmp_dist;
          best_n =j;
      }
    }
      indices[i] = best_n;
      values[i] = x[best_n];

  }

  return List::create(
    _["match_id"] = indices,
    _["match_x"] = values);
  
}
```

``` r
set.seed(1231)
ps_matchR(runif(5))
```

    $match_id
    [1] 4 5 1 5 2

    $match_x
    [1] 0.212552784 0.039240952 0.467915119 0.039240952 0.007559486

``` r
ps_match1(runif(5))
```

    $match_id
    [1] 4 2 1 2 0

    $match_x
    [1] 0.08654733 0.67072842 0.64116116 0.67072842 0.30754703

## Question 2: Things can be done faster

In the previous question, we have a double loop running twice over the
full set of observations. We need you to write the C++ so that the
computational complexity goes below `n^2`. (hint: Distance is symmetric)
