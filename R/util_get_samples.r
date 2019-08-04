# Generate two n-length vectors for validating various implementations

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Please provide N as a parameter, e.g. Rscript gen_samples.r 256", call.=FALSE)
}

N = as.numeric(args[1])

cat(toString(sample(N)))
