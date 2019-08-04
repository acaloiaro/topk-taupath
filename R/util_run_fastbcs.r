setwd("~/git/topk-taupath/R")
source('fastbcs.r')

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Please provide two equal-length vectors as input:\ne.g. Rscript --vanilla util_run_fastbcs.r 1,2,3,4 4,3,2,1", call.=FALSE)
}

x <- sapply(strsplit(args[1], ","), as.numeric)
y <- sapply(strsplit(args[2], ","), as.numeric)

sub1 <- function(num) {
    return(num-1)
}

pi <- get.pi3(x,y)
pi.sub1 <- sapply(pi, sub1)
cat(paste(pi.sub1, collapse = ", "), "\n")

