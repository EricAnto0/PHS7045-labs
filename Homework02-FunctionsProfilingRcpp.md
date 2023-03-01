Homework 02 - Functions, Profiling, and Rcpp
================

# Due date

Tuesday, February 28

# Background

For this assignment, you’ll be quested with speeding up some code using
what you have learned about vectorization and Rcpp.

## Part 1: Vectorizing code

The following functions can be written to be more efficient without
using parallel computing:

1.  This function generates a `n x k` dataset with all its entries
    distributed Poisson with mean `lambda`.

``` r
fun1 <- function(n = 100, k = 4, lambda = 4) {
  x <- NULL
  
  for (i in 1:n)
    x <- rbind(x, rpois(k, lambda))
  
  return(x)
}

fun1alt <- function(n = 100, k = 4, lambda = 4) {
  #sapply(1:k, function(g)rpois(n, lambda))
  matrix(rpois(k*n, lambda), nrow = n)
}

# Benchmarking
bench::mark(
  fun1(),
  fun1alt(), relative = TRUE, check = FALSE
)
```

    # A tibble: 2 × 6
      expression   min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr> <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 fun1()      16.7   31.7       1        62.4     1.48
    2 fun1alt()    1      1        39.9       1       1   

2.  Like before, speed up the following functions (it is OK to use
    StackOverflow)

``` r
# Total row sums
fun1 <- function(mat) {
  n <- nrow(mat)
  ans <- double(n) 
  for (i in 1:n) {
    ans[i] <- sum(mat[i, ])
  }
  ans
}

fun1alt <- function(mat){
  rowSums(mat)
}

# Cumulative sum by row
fun2 <- function(mat){
  n <- nrow(mat)
  k <- ncol(mat)
  ans <- mat
  for (i in 1:n) {
    for (j in 2:k) {
      ans[i,j] <- mat[i, j] + ans[i, j - 1]
    }
  }
  ans
}

fun2alt <- function(mat) {
  #t(cumsum(data.frame(t(dat))))
  t(apply(mat, 1, cumsum))
  #t(sapply(1:nrow(mat), function(k)cumsum(mat[k, ])))
}

# Use the data with this code
set.seed(2315)
dat <- matrix(rnorm(200 * 100), nrow = 200)
# Test for the first
bench::mark(
  fun1(dat),
  fun1alt(dat), relative = TRUE
)
```

    # A tibble: 2 × 6
      expression     min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr>   <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 fun1(dat)     5.18   8.04      1         196.     5.21
    2 fun1alt(dat)  1      1         8.21        1      1   

``` r
# Test for the second
bench::mark(
  fun2(dat),
  fun2alt(dat), relative = TRUE#, check = FALSE
)
```

    # A tibble: 2 × 6
      expression     min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr>   <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 fun2(dat)     4.49   3.01      1         1         NaN
    2 fun2alt(dat)  1      1         2.37      5.31      Inf

3.  Find the column max (hint: Check out the function `max.col()`).

``` r
# Find each column's max value
fun2 <- function(x){
  apply(x, 2, max)
}

fun2alt <- function(x){
  matrixStats::colMaxs(x)
  #diag(t(x)[, max.col(t(x))])
}
# Data Generating Process (10 x 10,000 matrix)
set.seed(1234)
x <- matrix(rnorm(1e4), nrow=10)

# Benchmarking
bench::mark(
  fun2(x),
  fun2alt(x), relative = TRUE, check = F
)
```

    # A tibble: 2 × 6
      expression   min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr> <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 fun2(x)     21.7   22.0       1        1        5.92
    2 fun2alt(x)   1      1        22.4      1.40     1   

## Part 2: Rcpp code

As we saw in the Rcpp week, vectorization may not be the best solution.
For this part, you must write a function using Rcpp that implements the
propensity score matching algorithm. You can use [Week 5’s
lab](https://github.com/UofUEpiBio/PHS7045-advanced-programming/issues/8#issuecomment-1424974938)
as a starting point for the problem. Your C++ file should look something
like the following:

In this C++ function, a control (i.e., untreated) is matched with more
than one treated individual such that each treated individual is matched
to only one control,

``` rcpp
#include<Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]
List psmatch(
  NumericVector pscores,
  LogicalVector is_treated
)
{
  /*... setup the problem creating the output...*/
 int n = static_cast<int>(is_treated.size());
  IntegerVector indices(n);
  NumericVector values(n);
    for (int i = 0; i < n; ++i) {
    //  Treated (is_treated == true) must be matched to control (is_treated == false)
    if(is_treated[i] == true){
    // Maximum value
    double cur_best = std::numeric_limits< double >::max();
    int cur_i = 0;
    
    for (int j = 0; j < n; ++j) {
      
        // grab treated == false and compare
        if (is_treated[j] == false){
          
        // If it is lower, then update
        if (std::abs(pscores[i] - pscores[j]) < cur_best) {
          
          cur_best = std::abs(pscores[i] - pscores[j]);
          cur_i    = j;
          
        }
      }
          
    }
    
    // register the result
    indices[i] = cur_i;
    values[i]  = pscores[cur_i];
    } 
  }
  
  // Returning
  return List::create(
   _["match_id"] = indices + 1, // We add one to match R's indices
    _["match_x"]  = values 
  );

}
```

## Example

``` r
set.seed(1001)
pscore = runif(20);
trt = rbinom(20, 1, 0.5)
res <- psmatch(pscore, trt)
cbind(pscore, trt, id = res$match_id, untreatedpscore = res$match_x)
```

               pscore trt id untreatedpscore
     [1,] 0.985688783   1 16      0.88129134
     [2,] 0.412628483   0  1      0.00000000
     [3,] 0.429539246   1  5      0.42650656
     [4,] 0.419172236   1  2      0.41262848
     [5,] 0.426506559   0  1      0.00000000
     [6,] 0.887797565   1 16      0.88129134
     [7,] 0.006096034   1  8      0.08121576
     [8,] 0.081215761   0  1      0.00000000
     [9,] 0.288657362   1 15      0.26666074
    [10,] 0.765342145   0  1      0.00000000
    [11,] 0.442924182   1 14      0.43987619
    [12,] 0.138363040   1 17      0.18208169
    [13,] 0.862473909   1 16      0.88129134
    [14,] 0.439876188   0  1      0.00000000
    [15,] 0.266660742   0  1      0.00000000
    [16,] 0.881291338   0  1      0.00000000
    [17,] 0.182081694   0  1      0.00000000
    [18,] 0.362768221   1  2      0.41262848
    [19,] 0.055500135   1  8      0.08121576
    [20,] 0.775185871   0  1      0.00000000
