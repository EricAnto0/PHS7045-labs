Lab 06 - Parallel computing
================

# Learning goals

- Learn from HW 1 and revisit vectorizing
- Parallel computing for simulations

# Lab task

**Part 1 (20 minutes)**:

Write an outline of [lab 2 / HW 1
solution](https://github.com/UofUEpiBio/PHS7045-advanced-programming-solutions).

Two people will present their outline. After the first presenter, the
second presenter will add any additional details/insights.

``` r
set.seed(NULL)
sample(c("Daniel", "Eric", "Haojia", "Hyejung", "Kline", "Linda", "Ransmond", "Ravi"),2)
```

**Part 2 (30 minutes)**:

Re-write design 1 and/or design 2 based upon what you’ve
learned/outlined from the solution. Compare efficiency between your
submission and your revision using the `bench` package.

**Part 3 (25 minutes)**:

Using either submitted design in HW 1 or revisions made in Part 2, use a
parallelization method to simulate 1000 replicates. Compare the timing
to simulating the same 1000 replicates without parallelization.

Extra time? Take any next steps on completing the simluation from the
original homework assignment:

- Assess the Type I error (increase the number of replicates)
- Refine $\delta$ as needed
- Simulate under the alternative hypothesis

------------------------------------------------------------------------

## Part 1

1.  Generate outcomes based on a set of treatment effects for each arm
    apriori

2.  For design one, the total number of samples are randomly assigned
    such that we have equal number of participants in each treatment arm
    (i.e., each arm will have equal allocation).

3.  Generate draws from posterior to compute and calculate probability
    each arm is greater than control

4.  Finally report the arm the maximum probablity greater than control

5.  For design two, 10 patients per arm for the first 40 patients are
    assigned.

6.  All have interims occurring at N = 40, 80, 120, 160, and 200
    patients.

7.  After each interim analysis, Generate draws from posterior and
    calculate pmax.

8.  Compute the randomization probabilities for the experimental arms as
    proportional to the probability that each experimental arm has the
    maximal responder rate.

``` r
# echo:false
N <- 228
arms <- 4
h0 <- 0.35
h1 <- 0.35
h2 <- 0.35
h3 <- 0.35



# Generate draws from posterior and calculate pmax
postDraws <- function(y,nMcmc,h0,h1,h2,h3,n0,n1,n2,n3){
  
  #------------------------------
  # Generate draws from posterior
  #------------------------------
  postDraws0 <- rbeta(n=nMcmc, shape1 = h0 + sum(y[1:n0,"0"]==1), shape2 = (1-h0) + n0 - sum(y[1:n0,"0"]==0))
  postDraws1 <- rbeta(n=nMcmc, shape1 = h1 + sum(y[1:n1,"1"]==1), shape2 = (1-h1) + n1 - sum(y[1:n1,"1"]==0))
  postDraws2 <- rbeta(n=nMcmc, shape1 = h2 + sum(y[1:n2,"2"]==1), shape2 = (1-h2) + n2 - sum(y[1:n2,"2"]==0))
  postDraws3 <- rbeta(n=nMcmc, shape1 = h3 + sum(y[1:n3,"3"]==1), shape2 = (1-h3) + n3 - sum(y[1:n3,"3"]==0))
  
  #-----------------------------------
  # Calculate allocation probabilities
  #-----------------------------------
  # v1-v3 probability each arm is best
  # v0 see RMatch Viele et al. 2020
  pBest <- pmax(postDraws0,postDraws1,postDraws2,postDraws3)
  
  v1 <- mean(pBest==postDraws1)
  v2 <- mean(pBest==postDraws2)
  v3 <- mean(pBest==postDraws3)
  
  v0 <- min(sum( c(v1,v2,v3) * (c( n1, n2, n3) + 1) / (n0 + 1), max(v1, v2, v3)) )
  
  # Standardize
  V0 <- v0 / (sum(v0,v1,v2,v3))
  V1 <- v1 / (sum(v0,v1,v2,v3))
  V2 <- v2 / (sum(v0,v1,v2,v3))
  V3 <- v3 / (sum(v0,v1,v2,v3))
  
  # Calculate probability each arm is greater than control
  p1 <- mean(postDraws1 > postDraws0)
  p2 <- mean(postDraws2 > postDraws0)
  p3 <- mean(postDraws3 > postDraws0)
  
  # Report maximum probablity an arm is greater than control
  pMax <- max(p1,p2,p3)
  # unname n0 objects for consistent object names in output
  n0 <- unname(n0)
  n1 <- unname(n1)
  n2 <- unname(n2)
  n3 <- unname(n3)
  out <- c(V0=V0,V1=V1,V2=V2,V3=V3,p1=p1,p2=p2,p3=p3,pMax=pMax,n0=n0,n1=n1,n2=n2,n3=n3)
  return(out)
  
}



design1 <- function(y, nMcmc=10000, h0, h1, h2, h3){
  
  # By the end of the study, each arm will have equal allocation
  postDraws(y=y,nMcmc=nMcmc,
            n0=nrow(y)/ncol(y),
            n1=nrow(y)/ncol(y),
            n2=nrow(y)/ncol(y),
            n3=nrow(y)/ncol(y),
            h0=h0, h1=h1, h2=h2, h3=h3)
  
}

design2 <- function(y, N, nInterim, h0, h1, h2, h3){
  
  #-----------------------------------------------
  # Set up parameters to inform interim monitoring
  #-----------------------------------------------

  # arms and looks as derived from y and look attributes
  arms         <- ncol(y)
  looks        <- floor(N / nInterim)
  
  # n at time t with n = 10 for each arm at t = 1
  nt           <- matrix(NA,nrow=looks,ncol=arms)
  colnames(nt) <- 0:(arms-1)
  nt[1,]       <- rep(nInterim / arms,arms)
  
  # Number of observations between each look
  # Account for possibility of no residual number of nInterim looks
  if(N %% nInterim != 0) {
    residual <- N %% nInterim
  } else {
    residual <- NULL
  }
  size         <- c(rep(nInterim,looks - 1),residual)
  
  
  
  #---------------------------------------------------------------
  # Update allocation probabilities and nt for each look iteration
  #---------------------------------------------------------------
  for(i in seq(looks-1)){
    
    alloProbs <- postDraws(y=y, nMcmc = 1000, h0=h0,h1=h1,h2=h2,h3=h3,
                           n0 = nt[i,"0"],
                           n1 = nt[i,"1"],
                           n2 = nt[i,"2"],
                           n3 = nt[i,"3"])
    
    nt[i+1, ] <- nt[i, ] + c(rmultinom(n = 1, size = size[i], 
                                    prob = alloProbs[c("V0","V1","V2","V3")]))
    
  }

  #--------------------------------------------------------
  # Report final probabilities of best arm and sample sizes
  #--------------------------------------------------------
  post <- postDraws(y=y, nMcmc = 1000, h0=h0,h1=h1,h2=h2,h3=h3,
                    n0 = nt[i,"0"],
                    n1 = nt[i,"1"],
                    n2 = nt[i,"2"],
                    n3 = nt[i,"3"])
  
  return(post)

}
```

``` r
# echo:false
getprobs <- function(norm_v, vdat, K = K, alpha_t, beta_t, realloc = TRUE){
  vdat <- as.data.frame(vdat)
  fprobs <- sapply(1:length(norm_v), function(m){
    rbeta(K, (alpha_t + sum(subset(vdat, trt_arm == m, select = y))), 
          (beta_t + nrow(subset(vdat, trt_arm == m))-sum(subset(vdat, 
                                      trt_arm == m, select = y))))})
  n_t <- as.vector(table(vdat[, "trt_arm"])[2:4])
  n_0 <- as.vector(table(vdat[, "trt_arm"])[1])
  V_t <- apply(fprobs, 1, function(x){which.max(x)})
  pt_max <- c(mean(V_t==1),
              mean(V_t==2),
              mean(V_t==3),
              mean(V_t==4)
  )
  if(realloc == TRUE){
    V_0 <- min(sum(pt_max[-1]*((n_t+1)/(n_0+1))), max(pt_max[-1]))
    norm_v <- sapply(c(V_0, pt_max[-1]), function(r){r/sum(c(V_0, pt_max[-1]))})
    return(norm_v)
  }else{
    return(which.max(pt_max)-1)
  }
}


sample_allocation <- function(norm_v, n0, N, alloc.type = c("F25", "RMATCH"), probs = probs, dat = dat){
  trt <- switch(alloc.type,
                "F25" = rep(c(1:length(norm_v)), N/length(norm_v)),
                "RMATCH" = sample(1:4, size = n0, replace = TRUE, prob = norm_v))
  ff <- lapply(unique(trt), function(x)cbind(x, sum(x==trt)))
  return(do.call(rbind, lapply(ff, function(m){cbind(trt_arm = m[1], 
                                                    y= sample(dat[, m[1]], m[2]))})))
}

mymodel <- function(threshold, N = N, K = 1000, design = 1, interim_n = c(40, 80, 120, 160, 200),
                    alpha_t = .35, beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65)){
  #set.seed(seed)
  mydatt <- NULL
  norm_v <- c(0.25, 0.25, .25, .25)
  names(norm_v) <- c("1", "2", "3", "4")
  dat <- do.call(cbind, lapply(1:length(norm_v), function(h){rbinom(1000, 1, probs[h])}))
  # dat <- data.frame(id = 1:10000, y = rbinom(10000, 1, .35))
  if(design == 1){
    #df <- lapply(split(sample(nrow(dat), N), c(0, 1, 2, 3)), function(i){dat[i, ]})
    mydatt <- sample_allocation(norm_v = norm_v, n0 = N, N = N, alloc.type = "F25", probs = probs, dat=dat)
  }else{
    #trt <- apply(t(rmultinom(40, size = 1, prob = as.vector(norm_v))), 1, function(x) which(x==1)-1)
    for(i in 1:(length(interim_n)+1)){
      norm_v <- c(0.25, 0.25, .25, .25)
      names(norm_v) <- c("1", "2", "3", "4")
      n0 <- ifelse(i<=(length(interim_n)), 40, N-(interim_n[i-1]))
      mydat <- sample_allocation(norm_v = norm_v, n0 = n0, N = N, alloc.type = "RMATCH", probs = probs, dat = dat)
      mydatt <- rbind(mydatt, mydat) 
      norm_v <- getprobs(norm_v = norm_v, vdat = mydatt, K = K, alpha_t = alpha_t, beta_t = beta_t, realloc = TRUE)
    }
  }
 mprobs <- getprobs(norm_v = norm_v, vdat = mydatt, K = K, alpha_t = alpha_t, beta_t = beta_t, realloc = FALSE)
 return(mprobs) 
}
```

Part 2 ::: {.cell}

``` r
N <- 228
arms <- 4
h0 <- 0.35
h1 <- 0.35
h2 <- 0.35
h3 <- 0.35
# Generate outcomes under null
y0 <- rbinom(N,size=1,prob=h0)
y1 <- rbinom(N,size=1,prob=h1)
y2 <- rbinom(N,size=1,prob=h2)
y3 <- rbinom(N,size=1,prob=h3)

y <- cbind("0" = y0, "1" = y1,"2" = y2,"3" = y3)


library(bench)
bench::mark("solution"=design1(y=y,h0=h0,h1=h1,h2=h2,h3=h3),
"mycode" =mymodel(threshold = .9912, N = 228, K = 1000, design = 1, alpha_t =.35,
        beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65)),
relative = TRUE, check = FALSE)
```

<div class="cell-output cell-output-stdout">

    # A tibble: 2 × 6
      expression   min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr> <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 solution    2.20   2.73      1         1.10      NaN
    2 mycode      1      1         2.59      1         Inf

</div>

``` r
bench::mark("solution" =design2(y=y,  N=N, nInterim = 40,
                                     h0=h0,h1=h1,h2=h2,h3=h3),
"mycode" =mymodel(threshold = .9912, N = 228, K = 1000, design = 2, alpha_t =.35,
        beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65)),
relative = TRUE, check = FALSE)
```

<div class="cell-output cell-output-stdout">

    # A tibble: 2 × 6
      expression   min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr> <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 solution    1      1         6.48      1        1   
    2 mycode      6.20   7.33      1         4.01     4.54

</div>

:::

# Part three

``` r
cores <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(cores)
# parallel::clusterExport(cl, list("threshold", "N", "K", "design","alpha_t", "beta_t", "probs", "mymodel", "sample_allocation", "getprobs"))
parallel::clusterExport(cl, list("mymodel", "sample_allocation", "y", "N", "getprobs", "postDraws", "design2"))

bench::mark("replicate" = replicate(100, {mymodel(threshold = .9912, N = 228, K = 1000, design = 2, alpha_t =.35, beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65))}),
            "parallelize" = parallel::parSapply(cl, 1:100, FUN = function(i)mymodel(threshold = .9912, N = 228, K = 1000, design = 2, alpha_t =.35, beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65))),
"solution" = parallel::parSapply(cl, 1:100, FUN = function(i)design2(y=y, N=N, nInterim = 40, h0=.35,h1=.35,h2=.35,h3=.35)), relative = TRUE, check = FALSE)
```

    Warning: Some expressions had a GC in every iteration; so filtering is
    disabled.

    # A tibble: 3 × 6
      expression    min median `itr/sec` mem_alloc `gc/sec`
      <bch:expr>  <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    1 replicate   68.0   61.0       1       891.        Inf
    2 parallelize  7.45   6.69      9.13      1.45      NaN
    3 solution     1      1        60.3       1         NaN

``` r
parallel::stopCluster(cl)
```
