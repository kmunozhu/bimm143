class07
================
Kiana Munoz
1/28/2020

Functions are here

``` r
is_both_na <- function(x, y) {
  if (length(x) != length(z)) {
    stop("your inputs are not the same length")
  }
  is.na(x) & is.na(y)
}


rescale_new <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)

  message("Found ", na.number, " NA's at position(s):", 
          paste(na.which, collapse=", ") ) 
  
  return( list(number=na.number, which=na.which) )
}


#Functions
function(x, na.rm=TRUE, plot=FALSE, ...) {
  # Our rescale function from the end of lecture 9

  if(na.rm) {
    rng <-range(x, na.rm=TRUE)
  } else {
    rng <-range(x)
  }

  answer <- (x - rng[1]) / (rng[2] - rng[1])
  if(plot) { 
    plot(answer, ...) 
  }

  return(answer)
}
```

    ## function(x, na.rm=TRUE, plot=FALSE, ...) {
    ##   # Our rescale function from the end of lecture 9
    ## 
    ##   if(na.rm) {
    ##     rng <-range(x, na.rm=TRUE)
    ##   } else {
    ##     rng <-range(x)
    ##   }
    ## 
    ##   answer <- (x - rng[1]) / (rng[2] - rng[1])
    ##   if(plot) { 
    ##     plot(answer, ...) 
    ##   }
    ## 
    ##   return(answer)
    ## }

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
z <- c(1,1)
#is_both_na(x, z)
```

``` r
# Start with a simple version of the problem
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
 exp=c(2,1,1),
 stringsAsFactors=FALSE)
df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
 exp=c(-2, NA, 1, 2),
 stringsAsFactors=FALSE)
```

``` r
x <- df1$IDs
y <- df2$IDs

intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
# x in y 
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
#y in x
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

``` r
y[y %in% x]
```

    ## [1] "gene2" "gene3"

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],  y[y %in% x])
}

gene_intersect2 <- function(df1, df2) {
  cbind(df1[df1$IDs %in% df2$IDs, ],  
        df2[df2$IDs %in% df1$IDs, ], "exp")
}

gene_intersect3 <- function(df1, df2, gene.colname = "IDs") {
  cbind(df1[df1[, gene.colname] %in% df2[, gene.colname],],  
        df2[df2[, gene.colname] %in% df1[, gene.colname],], "exp")
}
  
gene_intersect4 <- function(df1, df2, gene.colname = "IDs") {
    
    df1.name <- df1[, gene.colname]
    df2.name <- df2[, gene.colname]
    
    df1.inds <- df1.name %in% df2.name
    df2.inds <- df2.name %in% df1.name
    
  cbind(df1[df1.inds],
        df2[df2.inds], "exp")
}


df1$IDs
```

    ## [1] "gene1" "gene2" "gene3"

``` r
df1[, "IDs"]
```

    ## [1] "gene1" "gene2" "gene3"

``` r
gene.colname = "IDs"
df1[, gene.colname]  
```

    ## [1] "gene1" "gene2" "gene3"

``` r
x <- df1$IDs
y <- df2$IDs
gene_intersect(x, y)
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

``` r
gene_intersect2(df1, df2)
```

    ##     IDs exp   IDs exp "exp"
    ## 2 gene2   1 gene2  -2   exp
    ## 3 gene3   1 gene3   1   exp
