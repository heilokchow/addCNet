#' Non Parametric Estimate
#'
#' @description
#' This function allows you to do non parametric estimation for additive network model
#'
#' @param data The whole dataset containing pairwise communications, the first column
#' is the person with outgoing communication and the second column is the person with
#' incoming communication. The last column is the communication time point (should be
#' in numeric format).
#' @param zij The transfromed non-homogeneous effect (array with dimension (n,n,p,t)). The
#' transformation can be done using zProjection() function.
#' @param n number of individuals in the whole network.
#' @param p number of non-homogeneuous effect.
#' @param tz Time points for non-homofily effects (ranging from (0,1)), the length of this
#' vector should equal to the 4th dimension of zij which reflect the non-homofily effect
#' on or before those specific time points.
#' @param h1 Bandwidth for calculate \eqn{\alpha(t)} and \eqn{\beta(t)}, default is 0.05.
#'
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @useDynLib addCNet, .registration = TRUE
#'
#' @return
#' \item{homo_coefficients}{List Individual's cumulative effect from scaled time range
#' [0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
#' \item{homo_sd}{List Individual's standard deviation for cumulative effect from scaled
#' time range [0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
#' \item{nonhomo_coefficients}{List nonhomogeneous effect from scaled time range
#' [0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
#' \item{nonhomo_sd}{List nonhomogeneous parameter's standard deviation for cumulative
#' effect from scaled time range [0,1]. The parameters are estimated at (0.01, 0.02, ..., 1)}
#'
#' @export


nonParametric <- function(data, zij, n, p, tz = 1, h1 = 0.05, test = 0) {

  t1 = Sys.time()

  t12 = 0
  t23 = 0
  t34 = 0

  dim_check = dim(zij)
  n1 = ncol(data)
  t = nrow(data)

  t_max <- (data[t, n1] - data[1, n1])
  t_norm <- (data[, n1] - data[1, n1][[1]]) / t_max[[1]]
  t_sep_t = seq(0.01, 1, 0.01)

  # if (length(tz) != dim_check[4]) {
  #   stop("dimension not match!")
  # }

  itz = c()
  for (i in seq_len(length(tz))) {
    itz = c(itz, which.max(tz[i] < t_sep_t))
  }
  cat(itz)

  B = matrix(0, nrow = (2*n-1)+p, ncol = 100)
  b = matrix(0, nrow = (2*n-1)+p, ncol = 100)
  V = matrix(0, nrow = (2*n-1)+p, ncol = (2*n-1)+p)
  B1 = matrix(0, nrow = 2, ncol = 100)
  b1 = matrix(0, nrow = 2, ncol = 100)
  SDo = matrix(0, nrow = n, ncol = 100)
  SDi = matrix(0, nrow = n, ncol = 100)
  SDoa = matrix(0, nrow = n, ncol = 100)
  SDib = matrix(0, nrow = n, ncol = 100)
  N_t = matrix(0, nrow = n*(n-1), ncol = 1)
  N_tall = matrix(0, nrow = n*(n-1), ncol = 100)
  N_tallM = array(0, c(n, n, 100))
  Nij = matrix(0, nrow = n, ncol = n)

  t2 = Sys.time()
  temp = xConstruct(n, zij)
  if (p == 1) {
    P = cbind(temp$intercept, temp$x, matrix(temp$zij[,,1]))
  } else {
    P = cbind(temp$intercept, temp$x, temp$zij[,,1])
  }
  t3 = Sys.time()
  if (test == 0) {
    VP = pvp(P)
    Va = VP$Va
    Pa = VP$Pa
  } else {
    Va = t(P) %*% P
    Pa = solve(Va) %*% t(P)
  }
  t4 = Sys.time()

  hash = rep(0, n^2)
  co = 1
  co1 = 1
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        hash[co1] = co
        co = co + 1
      }
      co1 = co1 + 1
    }
  }

  co = 1
  co1 = 1
  for (i in seq_len(t)) {

    z1 = data[i, 1]
    z2 = data[i, 2]
    z3 = hash[(z1-1)*n+z2]
    N_t[z3, 1] = N_t[z3, 1] + 1
    Nij[z1, z2] = Nij[z1, z2] + 1

    Patemp = Pa[, z3]

    splinecalc(Patemp, N_tall, N_tallM, b, n, nrow(b), z1, z2, z3, h1, t_norm[i], t_sep_t)

    if (t_norm[i] >= t_sep_t[co]) {
      if (co == 1) {
        B[, co] = Pa %*% N_t
      } else {
        B[, co] = B[, co - 1] + Pa %*% N_t
      }
      B1[1, co] = -sum(B[2:n, co])
      B1[2, co] = -sum(B[(n+1):(2*n-1), co])
      SDo[, co] = 1/n^2 * colSums(Nij)
      SDi[, co] = 1/n^2 * rowSums(Nij)
      co = co + 1
      V = V + Va

      if (co == itz[co1]) {
        co1 = co1 + 1
        if (p == 1) {
          P = cbind(temp$intercept, temp$x, matrix(temp$zij[,,co1]))
        } else {
          P = cbind(temp$intercept, temp$x, temp$zij[,,co1])
        }
        Va = t(P) %*% P
        Pa = solve(Va) %*% t(P)
      }

      N_t = N_t * 0
    }
    if (co == 101) {
      break
    }

  }

  t12 = t12 + t2 - t1
  t23 = t23 + t3 - t2
  t34 = t34 + t4 - t3

  b1[1, ] = - colSums(b[2:n,])
  b1[2, ] = - colSums(b[(n+1):(2*n-1),])

  for (j in 1:100) {
    SDoa[, j] = 1/n^2 * colSums(N_tallM[,,j])
    SDib[, j] = 1/n^2 * rowSums(N_tallM[,,j])
  }

  if (!is.null(zij)) {
    output <- list(homo_coefficients = list(baseline = matrix(B[1,], nrow = 1),
                                            outgoing = rbind(matrix(B1[1,], nrow = 1),
                                                             B[2:n,]),
                                            incoming = rbind(matrix(B1[2,], nrow = 1),
                                                             B[(n+1):(2*n-1),]),
                                            sdout = SDo,
                                            sdinc = SDi),
                   nonhomo_coefficients = B[(2*n):(2*n-1+p),],
                   ab = list(baseline = matrix(b[1,], nrow = 1),
                             outgoing = rbind(matrix(b1[1,], nrow = 1),
                                              b[2:n,]),
                             incoming = rbind(matrix(b1[2,], nrow = 1),
                                              b[(n+1):(2*n-1),]),
                             sdout = SDoa,
                             sdinc = SDib),
                   th = b[(2*n):(2*n-1+p),],
                   V = V/100,
                   ts = c(t12, t23, t34))
  } else {
    output <- list(homo_coefficients = list(baseline = matrix(B[1,], nrow = 1),
                                            outgoing = rbind(matrix(B1[1,], nrow = 1),
                                                             B[2:n,]),
                                            incoming = rbind(matrix(B1[2,], nrow = 1),
                                                             B[(n+1):(2*n-1),])),
                   V = V/100)
  }
  return(output)
}

