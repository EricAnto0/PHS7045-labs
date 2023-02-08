Lab 03 - Functions and data.table
================

# Learning goals

- Used advanced features of functions in R.
- Use the `merge()` function to join two datasets.
- Deal with missings and data imputation data.
- Identify relevant observations using `quantile()`.
- Practice your GitHub skills.

# Lab description

For this lab, we will deal with the meteorological dataset downloaded
from the NOAA, the `met`. We will use `data.table` to answer some
questions regarding the `met` data set, and practice our Git+GitHub
skills.

This markdown document should be rendered using `gfm` document.

# Part 1: Setup the Git project and the GitHub repository

1.  Go to your documents (or wherever you are planning to store the
    data) in your computer, and create a folder for this project, for
    example, “PHS7045-labs”

2.  In that folder, save [this
    template](https://raw.githubusercontent.com/UofUEpiBio/PHS7045-advanced-programming/main/labs/03-functions-and-datatable/03-functions-and-datatable.qmd)
    as “README.qmd.” This will be the markdown file where all the magic
    will happen.

3.  Go to your GitHub account and create a new repository, hopefully of
    the same name this folder has, i.e., “PHS7045-labs”.

4.  Initialize the Git project, add the “README.qmd” file, and make your
    first commit.

5.  Add the repo you just created on GitHub.com to the list of remotes,
    and push your commit to `origin` while setting the upstream.

Most of the steps can be done using the command line:

``` sh
# Step 1
cd ~/Documents
mkdir PHS7045-labs
cd PHS7045-labs

# Step 2
wget https://raw.githubusercontent.com/UofUEpiBio/PHS7045-advanced-programming/main/labs/03-functions-and-datatable/03-functions-and-datatable.qmd
mv 03-functions-and-datatable.qmd README.qmd

# Step 3
# Happens on github

# Step 4
git init
git add README.Rmd
git commit -m "First commit"

# Step 5
git remote add origin git@github.com:[username]/PHS7045
git push -u origin master
```

You can also complete the steps in R (replace with your paths/username
when needed)

``` r
# Step 1
setwd("~/Documents")
dir.create("PHS7045-labs")
setwd("PHS7045-labs")

# Step 2
download.file(
  "https://raw.githubusercontent.com/UofUEpiBio/PHS7045-advanced-programming/main/labs/03-functions-and-datatable/03-functions-and-datatable.qmd",
  destfile = "README.qmd"
  )

# Step 3: Happens on Github

# Step 4
system("git init && git add README.Rmd")
system('git commit -m "First commit"')

# Step 5
system("git remote add origin git@github.com:[username]/PHS7045-labs")
system("git push -u origin master")
```

Once you are done setting up the project, you can now start working on
the lab.

# Part 2: Advanced functions

## Question 1: **Ellipsis**

Write a function using the ellipsis argument (`...`) with the goal of
(i) retrieving the list of arguments passed to it, (ii) printing
information about them using `str()`, and (iii) printing the environment
where they belong and the address of the object in memory using
`data.table::address()`.

``` r
library(data.table)
library(R.utils)
```

    Loading required package: R.oo

    Loading required package: R.methodsS3

    R.methodsS3 v1.8.2 (2022-06-13 22:00:14 UTC) successfully loaded. See ?R.methodsS3 for help.

    R.oo v1.25.0 (2022-06-12 02:20:02 UTC) successfully loaded. See ?R.oo for help.


    Attaching package: 'R.oo'

    The following object is masked from 'package:R.methodsS3':

        throw

    The following objects are masked from 'package:methods':

        getClasses, getMethods

    The following objects are masked from 'package:base':

        attach, detach, load, save

    R.utils v2.12.2 (2022-11-11 22:00:03 UTC) successfully loaded. See ?R.utils for help.


    Attaching package: 'R.utils'

    The following object is masked from 'package:utils':

        timestamp

    The following objects are masked from 'package:base':

        cat, commandArgs, getOption, isOpen, nullfile, parse, warnings

``` r
#plot_x <- function(x, y, ...){
#plot(x = x, y = y, main = "Experiment", col = "red")
#}
my_e <- function(...){
  val <- list(...)
  lapply(val ,str)
  lapply(val, \(x){
    print(data.table(environment(x)))})
  lapply(val, \(x){
    print(data.table::address(x))})
  invisible()
}
my_e(123, pi)
```

     num 123
     num 3.14
    Null data.table (0 rows and 0 cols)
    Null data.table (0 rows and 0 cols)
    [1] "00000223832ed900"
    [1] "0000022383f6a238"

Knit the document, commit your changes, and push them to GitHub.

## Question 2: **Lazy evaluation**

A concept we did not review was lazy evaluation. Write a function with
two arguments (`a` and `b`) that only uses one of them as an integer,
and then call the function passing the following arguments
`(1, this_stuff)`

``` r
lazy_eval <- function(x, y){
  print(x)
}
lazy_eval(1, this_stuff)
```

    [1] 1

Knit the document, commit your changes, and push them to GitHub.

## Question 3: **Putting all together**

Write a function that fits a linear regression model and saves the
result to the global environment using the `assign()` function. The name
of the output must be passed as a symbol using lazy evaluation.

``` r
lreg <- function(x, y, ...){
 #innerf <- list(...)
  assign("lfit", lm(y ~ x), envir = .GlobalEnv)
 #ls(environment(innerf))
}
x <- rnorm(100)
y <- rbinom(100, 1, 0.5)
lreg(x, y)  
```

Knit the document, commit your changes, and push them to GitHub.

# Part 3: Data.table

## Setup in R

1.  Load the `data.table` (and the `dtplyr` and `dplyr` packages if you
    plan to work with those).

2.  Load the met data from
    https://raw.githubusercontent.com/USCbiostats/data-science-data/master/02_met/met_all.gz,
    and the station data. For the latter, you can use the code we used
    during the lecture to pre-process the stations’ data:

``` r
# Download the data
stations <- fread("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv")
metdata <- fread("https://raw.githubusercontent.com/USCbiostats/data-science-data/master/02_met/met_all.gz")
stations[, USAF := as.integer(USAF)]

# Dealing with NAs and 999999
stations[, USAF   := fifelse(USAF == 999999, NA_integer_, USAF)]
stations[, CTRY   := fifelse(CTRY == "", NA_character_, CTRY)]
stations[, STATE  := fifelse(STATE == "", NA_character_, STATE)]

# Selecting the three relevant columns, and keeping unique records
stations <- unique(stations[, list(USAF, CTRY, STATE)])

# Dropping NAs
stations <- stations[!is.na(USAF)]

# Removing duplicates
stations[, n := 1:.N, by = .(USAF)]
stations <- stations[n == 1,][, n := NULL]

datta <- merge(
  # Data
  x     = metdata,      
  y     = stations, 
  # List of variables to match
  by.x  = "USAFID",
  by.y  = "USAF", 
  # Which obs to keep?
  all.x = TRUE,      
  all.y = FALSE
  )
dim(datta)
```

3.  Merge the data as we did during the lecture.

## Question 1: Representative station for the US

What is the median station in terms of temperature, wind speed, and
atmospheric pressure? Look for the three weather stations that best
represent the continental US using the `quantile()` function. Do these
three coincide?

``` r
# datta[, .(median_temp = quantile(temp, 0.5, na.rm = TRUE),
#         median_wind.sp = quantile(wind.sp, 0.5, na.rm = TRUE),
#         median_atm.press = quantile(atm.press, 0.5, na.rm = TRUE)),
#     keyby = .(STATE)][order(STATE)]
```

Knit the document, commit your changes, and Save it on GitHub. Don’t
forget to add `README.md` to the tree, the first time you render it.

## Question 2: Representative station per state

Identify what the most representative (the median) station per state is.
Instead of looking at one variable at a time, look at the euclidean
distance. If multiple stations show in the median, select the one at the
lowest latitude.

Knit the doc and save it on GitHub.

## (optional) Question 3: In the middle?

For each state, identify the closest station to the mid-point of the
state. Combining these with the stations you identified in the previous
question, use `leaflet()` to visualize all \~100 points in the same
figure, applying different colors for those identified in this question.

Knit the doc and save it on GitHub.

## (optional) Question 4: Means of means

Using the `quantile()` function, generate a summary table that shows the
number of states included, average temperature, wind speed, and
atmospheric pressure by the variable “average temperature level,” which
you’ll need to create.

Start by computing the states’ average temperature. Use that measurement
to classify them according to the following criteria:

- low: temp \< 20
- Mid: temp \>= 20 and temp \< 25
- High: temp \>= 25

Once you are done with that, you can compute the following:

- Number of entries (records),
- Number of NA entries,
- Number of stations,
- Number of states included, and
- Mean temperature, wind speed, and atmospheric pressure.

All by the levels described before.

Knit the document, commit your changes, and push them to GitHub. If
you’d like, you can take this time to include the link of [the issue of
the
week](https://github.com/UofUEpiBio/PHS7045-advanced-programming/issues/5)
so that you let us know when you are done, e.g.,

``` bash
git commit -a -m "Finalizing lab 3 https://github.com/UofUEpiBio/PHS7045-advanced-programming/issues/5"
```
