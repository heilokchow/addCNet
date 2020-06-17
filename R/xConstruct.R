xConstruct <- function(n, A) {
  x = array(0, dim = c((n*(n-1)), 2*(n - 1)))
  p = dim(A)[3]
  zij = matrix(0, nrow = n*(n-1), ncol = p)
  one = matrix(1, nrow = n*(n-1), ncol = 1)

  A = zProjection(A)

  co = 1
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {

        for (q in seq_len(p)) {
          zij[co, q] = A[i, j, q]
        }

        if (i == 1) {
          for (k in 1:(n-1)) {
            x[co, k] = -1
          }
        } else {
          x[co, i-1] = 1
        }
        if (j == 1) {
          for (k in 1:(n-1)) {
            x[co, k+(n-1)] = -1
          }
        } else {
          x[co, j-1+(n-1)] = 1
        }
        co = co + 1
      }
    }
  }

  list(intercept = matrix(1, nrow = n*(n-1), ncol = 1), x = x, zij = zij)
}
