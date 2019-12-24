setwd("~/git/topk-taupath/R")
source('fastbcs.r')
N <- 512

x = c(1, 2, 3, 4)
y = c(1, 2, 3, 4)

sub1 <- function(num) {
    return(num-1)
}

pi <- fastbcs(x,y)
cat(sapply(pi, sub1), "\n")
pi2 <- fastbcs(x,y)
cat(pi2, "\n")

cat("Equivalent?", all(pi == pi2), "\n")



