# The canonical FastBCS algorithm
fastbcs <- function(x, y) {
  n <- length(x)
  id <- 1:n
  pie <- id
  tie.list = vector("list", n)

  # create the initial concordance matrix of x,y
  concord <- get.c(x, y)
  i <- n

  permute <- TRUE
  while (permute) {
    # get the permuted concordance matrix for the current π
    c.ia = concord[pie, pie]

    # calculate column sums for the current stage and find the minimum
    tail.i <- pie[1:i]
    colsums.i <- colSums(concord[tail.i, tail.i])

    # find any stage k prior to stage i which contain column i as a member of its tie set.
    ties.k = ties.above(i, pie[i], tie.list)

    # l = the permuted index of the column with the minimum sum in stage i
    min.i <- tail.i[order(colsums.i)][1]
    l = match(min.i, pie)

    # check how many other column sums tie with the minimum
    minColsum.i <- min(colsums.i)
    n.tie <- sum(colsums.i == minColsum.i)
    if (n.tie > 1) {
      # set TRUE in the row of the ith column for each column_j that ties with the minimum
      tie.id <- tail.i[colsums.i == minColsum.i]
      tie.list[[i]] = tie.id
    }

    # transpose the minimum column with pip[i]
    l <- id[pie == min.i]
    pie[pie == min.i] <- pie[i]
    pie[i] <- min.i

    switch <- FALSE

    # determine if a forward step is necessary
    if (length(ties.k) > 0) {
      c.i <- concord[pie, pie]

      for (k in ties.k) {
        if (k > i) {
          old <- cumsum.k(c.i[, i], i, k)
          new <- cumsum.k(c.i[, k], i, k)

          if (all(new >= old) & any(new > old)) {
            temp <- pie[k]
            pie[k] <- pie[i]
            pie[i] <- temp

            switch <- TRUE

            i <- k - 1
            pk = pie[1:k]
            tie.list[pk] = NULL

            break
          }
        }
      }
    }

    # take a backwoard step
    if (switch == FALSE) {
      i <- i - 1
    }

    # AC (07-20-2019) - I've forgotten a lot about this algorithm; I believe this means the remaining
    # columns are fully concordant and that further iterations would be futile.
    if (sum(concord[pie[1:i], pie[1:i]]) == i * (i - 1)) {
      c.ia = concord[pie, pie]
      permute <- FALSE
    }
  }

  return(pie = pie)
}

get.c <- function(x, y) {
  N <- length(x)
  c <- matrix(data = rep(0, N * N),
              nrow = N,
              ncol = N)

  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        product <- ((y[i] - y[j]) * (x[i] - x[j]))


        if (product > 0) {
          c[i, j] = 1
        } else if (product < 0) {
          c[i, j] = -1
        }
      }
    }
  }
  return(c)
}

# Calculate the cumulative sum
cumsum.k <- function(vec, lo, hi) {
  return(c(cumsum = cumsum(vec)[lo:hi]))
}

# currentStage: the stage being examined
# x: integer being tested for in list of lists
# ll: list of lists
ties.above = function(currentStage, x, ll) {
  ties = c()
  if (currentStage < length(ll)) {
    for (j in (currentStage + 1):length(ll)) {
      if (x %in% ll[[j]]) {
        ties = c(ties, j)
      }
    }
  }
  rev(ties)
}
