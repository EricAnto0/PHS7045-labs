R Essentials: loops, functions, and \*apply functions
================

# Background

Response-adaptive randomization (RAR) has been used in precision
medicine trials, such as the trials
[I-SPY-2](https://www.nejm.org/doi/pdf/10.1056/NEJMoa1513750?articleTools=true),
[BATTLE-1](https://aacrjournals.org/cancerdiscovery/article/1/1/44/2198?casa_token=pK1gZcX-FgkAAAAA:KmsD6qnoaOMxqHJlg0VGlmqr2nqIl49Xupuh0FX7nnJXNjtdBwVWsdmVtUIXKdEWQ_e5i9pG),
and
[BATTLE-2](https://cdn.amegroups.cn/journals/amepc/files/journals/12/articles/6846/public/6846-PB2-R2.pdf)
to gather early evidence of treatment arms that work best for a given
biomarker. Throughout RAR, the treatment allocation adjusts depending on
which treatment arm looks most promising. RAR is criticized for the
following reasons
[(Korn)](https://watermark.silverchair.com/djx013.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAvswggL3BgkqhkiG9w0BBwagggLoMIIC5AIBADCCAt0GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMUKqYjuX4J38OrUAJAgEQgIICrkuWgqA03alHxi4xfUz3ybPTepY21z2PaOZFRSmRgKhEf4ROr5tIP8MoeLwYxIexRpjXLTjKQ_Yn6s2oYBYUTJEG_gVqeevFGwvxyi8zas8R1GyiQWDyjnhBsjKcAsKCl6M63_K-LXEPLj7nyw4FZEefXgMSO91jz_bKseNLmhDEmZ7oPvO_xlBzWptEW4gE8I7mrv8h3NbFrM0BrbDIeOkIr36F3X-F3ip5a8lEDirZwrKyLCb598ivER1JVlr7iCzNV9OYHv6P79E6YxKGxJPDNk1WJ-ClcUb6aA87-Gae81SVipvqhu2npbLCFtWVikZKliWQ5lPE9CieOB9hnd93HLsp7Da9U3D0jodTMbk9RE28I7kzuoWHO8BEoPK8p1_1nTWw6pV5CbS02RFxtTr4-aCMm364vJAHrgP3rNmvi-KEGB4rxufA_hebdi2yS20nll2HdAVT81m_cm87x3YvVP7Ut9H-0zTCAVYACBJJUjKS7w7QKkIN-FqdQ7N7u_YJ83bl69sMF0fddIVTA6WFAJWpVmiij5ou34auXuzDtEKLN5nmQ-467HZtpMAkoyFxcyCdtI_ieuyCp1jlP3IC07JZAWvA_8A2W8hEFaRBimq0STpclvxEgt2OCDtE9kTXvS8v6ArS4U6ZkuM7Uw3PSG1ISoYAWPPvWzGMKWU4W4nCS0rT-TPWPJkPf8TGkUK9Q1cDemuYlWZBrYN-RUJ8X204vbqa-p9jsRKIfsQ0QDrW97qHYxf-659kKav288WZ8d5wmESsLKyHEc8JatY2dvmsIMIdeK3mKA68qDUqS3C1yO7X0T4BJZpJPGjjPzFrrqaLAmLdhdqITE5AcT_Eiy3upRuYw3TdqbznbDog3k2K2yIgQUZwu2HC9q_fJhqmxh0WkHQDyYeiXyqb):

1.  Possible bias: Time trends in participant enrollment (example:
    healthier participants earlier) may lead to a bias in estimating
    which treatment is superior
2.  Possible inefficiency: Unequal allocation may lead to statistical
    inefficiencies compared to equal allocation
3.  Possibly unethical: A moderately large sample size could be enrolled
    onto an arm worse than control

It’s not always bad for a method to be criticized – it means there are
open questions to be addressed by new/improved methods. For example, the
manuscript [Comparison of methods for control allocation in multiple arm
studies using response adaptive
randomization](https://journals.sagepub.com/doi/pdf/10.1177/1740774519877836)
develops and compares methods to improve RAR efficiency compared to
equal allocation by maintaining a reasonable number of participants
allocated to the control arm.

# Assignment

Replicate the results in Table 2 of [Comparison of methods for control
allocation in multiple arm studies using response adaptive
randomization](https://journals.sagepub.com/doi/pdf/10.1177/1740774519877836)
for the designs RMatch and F25. Use only `base` R functions and use
functions to avoid duplicating multiple lines of code (i.e., this is
called modularizing your code).

In this assignment, you will practice writing functions, loops, and
using the \*apply functions.

# Guidance

1.  The manuscript uses a minimally informative, normal prior for the
    log-odds. Rather than this prior, use a Beta(0.35, 0.65) prior on
    the response rate for each arm. This is minimally informative
    favoring the null hypothesis and allows you to use the Beta-Binomial
    conjugacy to obtain a Beta posterior distribution. (Refer to [lab
    2](https://uofuepibio.github.io/PHS7045-advanced-programming/week-02-lab.html)
    for notation of the posterior distribution).

2.  A way to estimate $P_t$(Max) is to `cbind` K = \[1000+\] draws from
    the posterior distribution of each arm and to see how frequently
    (across the K draws from each arm) each arm is drawn to be the
    largest.

3.  The value of $\delta$ was found through simulation under the null
    scenario such that an efficacious trial was declared to be found
    only 2.5% of the time. The manuscript provides values for $\delta$
    though verify and modify as needed (there could be slight
    differences by using a different prior).

4.  The manuscript uses 100K replicates to determine $\delta$ and
    estimate the values in Table 2. Start with a small number of
    replicates (1K or 10K) to make sure the code is running correctly
    and your results are in the ball park of being similar to the
    manuscript. Ramp up the number of replicates as feasible (for this
    assignment feasible is a run time of no longer than 20 minutes).

5.  Please talk to us if you are having difficulties with this
    assignment. We realize it is in a space that may be new to you, and
    it is not the intent for this assignment to take more than 7-10
    hours over the course of 2 weeks.

``` r
#### message: false
getprobs <- function(norm_v, vdat, K = K, alpha_t, beta_t, realloc = TRUE){
  vdat <- as.data.frame(vdat)
  fprobs <- sapply(1:length(norm_v), function(m){
    rbeta(K, (alpha_t + sum(subset(vdat, trt_arm == m, select = y))), 
          (beta_t + nrow(subset(vdat, trt_arm == m))-sum(subset(vdat, trt_arm == m, select = y))))})
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
    return(list(fprobs, (which.max(pt_max)-1)))
  }
}


sample_allocation <- function(norm_v, n0, N, alloc.type = c("F25", "RMATCH"), probs = probs, dat){
  trt <- switch(alloc.type,
                "F25" = rep(c(1:length(norm_v)), N/length(norm_v)),
                "RMATCH" = sample(1:4, size = n0, replace = TRUE, prob = norm_v))
  ff <- (lapply(unique(trt), function(x)cbind(x, sum(x==trt))))
  df <- do.call(rbind, lapply(ff, function(m){cbind(trt_arm = m[1], 
                                                    y= sample(dat[, m[1]], m[2]))}))
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
    mydatt <- sample_allocation(norm_v = norm_v, n0 = N, N = N, alloc.type = "F25", probs = probs, dat)
  }else{
    #trt <- apply(t(rmultinom(40, size = 1, prob = as.vector(norm_v))), 1, function(x) which(x==1)-1)
    for(i in 1:(length(interim_n)+1)){
      norm_v <- c(0.25, 0.25, .25, .25)
      names(norm_v) <- c("1", "2", "3", "4")
      n0 <- ifelse(i<=(length(interim_n)), 40, N-(interim_n[i-1]))
      mydat <- sample_allocation(norm_v = norm_v, n0 = n0, N = N, alloc.type = "RMATCH", probs = probs, dat)
      mydatt <- rbind(mydatt, mydat) 
      norm_v <- getprobs(norm_v = norm_v, vdat = mydatt, K = K, alpha_t = alpha_t, beta_t = beta_t, realloc = TRUE)
    }
    
    
    mprobs <- getprobs(norm_v = norm_v, vdat = mydatt, K = K, alpha_t = alpha_t, beta_t = beta_t, realloc = FALSE)
    return(mprobs[[2]])
  }
  
}
getres1 <- replicate(1000, {mymodel(threshold = .9912, N = 228, K = 1000, design = 1,
alpha_t =.35, beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65))})

getres2 <- replicate(1000, {mymodel(threshold = .9912, N = 228, K = 1000, design = 2, alpha_t =.35, beta_t = .65, probs = c(0.35, 0.45, 0.55, 0.65))})
 

getr1 <- cumsum(rev(prop.table(table(getres1))))*100 
getr12 <- cumsum(rev(prop.table(table(getres2))))*100
tab <- rbind(getr1[1:3], getr12[1:3])
colnames(tab)  <- c("Pr(pick arm 3) = Pr(pick best arm) (%)",
"Pr(pick arm 2 or better) (%)", "Pr(pick arm 1 or better) (%)")
rownames(tab) <- c("F25", "RMatch")
```

## Performance of each design on arm selection in the mixed scenario.

           Pr(pick arm 3) = Pr(pick best arm) (%) Pr(pick arm 2 or better) (%)
    F25                                      12.5                         25.0
    RMatch                                   84.8                         99.1
           Pr(pick arm 1 or better) (%)
    F25                            37.5
    RMatch                        100.0
