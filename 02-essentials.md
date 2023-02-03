Lab 02/Lab 04 - R Essentials, Debugging and Profiling
================

# Learning goals

- Functions, loops, and vectorized programming
- Writing and simulation a function from a statistical paper
- Debugging, profiling, and improving efficiency

# Lab task 4

We did not have enough time for most to finish [lab
2](https://uofuepibio.github.io/PHS7045-advanced-programming/week-02-lab.html),
and today’s lab time is dedicating to finishing the lab.

If you complete the lab 2, profiling your functions and see if there are
rooms to improve efficiency.

Devote any extra time to working on [homework
1](https://uofuepibio.github.io/PHS7045-advanced-programming/hw-01-essentialsSimulations.html)
which is a continuation of the lab.

Today’s lab is the first steps toward designing a response-adaptive
randomization (RAR) trial. RAR designs are used in precision medicine
trials, such as the [BATTLE
trial](https://aacrjournals.org/cancerdiscovery/article/1/1/44/2198?casa_token=pK1gZcX-FgkAAAAA:KmsD6qnoaOMxqHJlg0VGlmqr2nqIl49Xupuh0FX7nnJXNjtdBwVWsdmVtUIXKdEWQ_e5i9pG),
to gather early evidence of treatment arms that work best for a given
biomarker. Throughout RAR, the treatment allocation adjusts depending on
which treatment arm looks most promising. We will focus on the initial
steps of coding this design.

The lab is motivated by the paper by Kurt Viele: [Comparison of methods
for control allocation in multiple arm studies using response adaptive
randomization](https://journals.sagepub.com/doi/pdf/10.1177/1740774519877836).
It practices:

- Pre-allocating vectors
- Using a loop
- Writing a function

# Notation

Notation and criteria for study to successfully declare a treatment as
efficacious:

- $i = 1, \dots, N$ participants
- $t = 0, 1, 2, 3$ study arms ($t = 0$ is control)
- $Y_i \mid t \sim Bern(p_t)$ and $y_t$ is a vector of $n_t$ observed
  outcomes on arm $t$
- The prior on $p_t \sim Beta(\alpha_t, \beta_t)$

Posterior Distribution
$Pr(p_t | y_t) \sim Beta(\alpha_t + \sum y_t, \beta_t + n_t - \sum y_t)$

Quoting from the paper: The trial is considered successful at the final
analysis if there is a high posterior probability that at least one arm
has a higher rate than control.

$\text{max}_t Pr( p_t > p_0 ) > \delta$

where $\delta$ is a threshold chosen to maintain familywise type I error
for the study at one-sided 2.5%.

# Consider 2 different designs:

1.  Equal allocation to four arms throughout design.
2.  RAR where the allocation probability is updated at an interim
    analysis as follows:

- $V_t = P_t$ (Max)
- $V_0 = min\{\sum V_t \frac{(n_t + 1)}{(n_0 + 1)}, max(V_1, V_2, V_3) \}$

$V_0, V_1, V_2,$ and $V_3$ are renormalized to sum to 1 and are
allocation probabilities.

Note: A way to estimate $P_t$(Max) is to `cbind` K = \[1000\] draws from
the posterior distribution of each arm and to see how frequently (across
the K draws from each arm) each arm is drawn to be the largest.

# Lab task:

Write a function for each study design to simulate one trial.

- N = 228 with interim analyses after every 40th participant starting at
  40.
- Use equal allocation for first 40 patients for both designs.
- Assume a setting where treatment effect is 0.35 for each study arm
  (the null scenario). (But allow flexibility in function for other
  treatment effects).
- $\alpha_t = 0.35$ for all $t$ and $\beta_t = 0.65$ for all arms.
- Use the following $\delta$ thresholds to determine a successful trial:
  - Design 1, $\delta = 0.9912$
  - Design 2, $\delta = 0.9892$

For simplicity, have your function return a list of at least the
following output:

1.  The probability that the best treatment arm is better than control.
2.  The number of patients assigned to each treatment arm.

``` r
mymodel <- function(threshold, N = N, K = 1000, design, 
                    alpha_t =.35, beta_t = .65, seed = 1500){
  set.seed(seed)
  interim_n <- c(40, 80, 120, 160, 200)
  
  mydatt <- NULL
  norm_v <- c(0.25, 0.25, .25, .25)
  names(norm_v) <- c("0", "1", "2", "3")
  
  dat <- data.frame(id = 1:10000, y = rbinom(10000, 1, .35))
  if(design ==1){
    df <- lapply(split(sample(nrow(dat), N), c(0, 1, 2, 3)), function(i){dat[i, ]})
    mydatt <- cbind(do.call(rbind, df), 
                                trt_arm = rep(names(df), 
                                              sapply(df, function(y){dim(y)[1]})))
  }else{
  norm_v <- c(0.25, 0.25, .25, .25)
  names(norm_v) <- c("0", "1", "2", "3")
trt <- apply(t(rmultinom(40, size = 1, prob = as.vector(norm_v))), 1, function(x) which(x==1)-1)
  for(i in 1:(length(interim_n)+1)){
    if(i<=(length(interim_n))){
       n0 <- 40
       df <- lapply(split(sample(nrow(dat), n0), trt), function(j){dat[j, ]})  
      }else{
        n0 <- N-(interim_n[i-1])
        df <- lapply(split(sample(nrow(dat), n0), trt), function(j){dat[j, ]})  
      }
    mydat <- cbind(do.call(rbind, df), 
                                trt_arm = rep(names(df), 
                                              sapply(df, function(y){dim(y)[1]})))
    post_probs <- sapply(df, 
                         function(m){rbeta(K, 
                                           (alpha_t + sum(m[, "y"])),
                                           (beta_t + nrow(m) - sum(m[, "y"])))})
    
    n_t <- as.vector(sample_per_trt <- sapply(df, function(m){nrow(m)})[2:4])
    #n0 <- as.vector(sample_per_trt <- sapply(df, function(m){nrow(m)})[-(2:4)])
    V_t <- apply(post_probs, 1, function(x){which.max(x)})
    pt_max <- c(mean(V_t==1),
              mean(V_t==2),
              mean(V_t==3),
              mean(V_t==4)
              )
    
    V_0 <- min(sum(pt_max[-1]*((n_t+1)/(n0+1))), max(pt_max[-1]))
    norm_v <- sapply(c(V_0, pt_max[-1]), function(r)r*1/sum(c(V_0, pt_max[-1])))
 mydatt <- rbind(mydatt, mydat)  
    }
  }
  
  fprobs <- sapply(0:3, function(m){rbeta(K, (alpha_t + sum(mydatt["trt_arm"==m, "y"])),
            (beta_t + nrow(mydatt["trt_arm"==m, ]) - sum(mydatt["trt_arm"==m, "y"])))})
  sample_per_trt <- table(mydatt$trt_arm)
  names(sample_per_trt) <- c("Control 1", "Arm 1", "Arm 3", "Arm 4")
  colnames(fprobs) <- Hmisc::Cs(p_0, p_1, p_2, p_3)
  trsholdprob <- apply(fprobs[, 2:4], 2, function(r){
    mean(r > fprobs[, 1])})#; print(trsholdprob) 
  return(list(paste("The probability that the best treatment arm",
                    c(1:3)[trsholdprob==max(trsholdprob)], 
                    "is better than control is",
                    max(trsholdprob)), 
              ifelse(max(trsholdprob) > threshold, 
                     "Trial was successsful", 
                     "Trial was not successsful"),
              sample_per_trt))
}
mymodel(threshold = .9912, N = 228, K = 1000, design = 2,
        alpha_t =.35, beta_t = .65)
```

    [[1]]
    [1] "The probability that the best treatment arm 3 is better than control is 0.507"

    [[2]]
    [1] "Trial was not successsful"

    [[3]]
    Control 1     Arm 1     Arm 3     Arm 4 
           63        34        64        67 

``` r
bench::mark(mymodel(threshold = .9912, N = 228, K = 1000, design = 2,
        alpha_t =.35, beta_t = .65))[c("expression", "min", "median", "itr/sec", "n_gc")]
```

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    Warning in split.default(sample(nrow(dat), n0), trt): data length is not a
    multiple of split variable

    # A tibble: 1 × 4
    # … with 4 more variables: expression <bch:expr>, min <bch:tm>,
    #   median <bch:tm>, `itr/sec` <dbl>

# If you have more time

- Replicate the design many (10K) times. Calculate the Type I error.
- Find $\delta$ for each design (supposing you didn’t already know it).
- Replicate the study design assuming treatment effects of
  - $p_0 = p_1 = p_2 = 0.35$
  - $p_3 = 0.65$
