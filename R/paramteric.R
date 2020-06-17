#' Semi Parametric Estimate
#'
#' @description
#' This function allows you to do non parametric estimation for additive network model
#'
#' @param data The whole dataset containing pairwise communications, the first column
#' is the person with outgoing communication and the second column is the person with
#' incoming communication. The last column is the communication time point (should be
#' in numeric format).
#' @param zij The transfromed non-homogeneous effect (array with dimension (n,n,p)). The
#' transformation can be done using zProjection() function.
#' @param n number of individuals in the whole network.
#' @param p number of non-homogeneuous effect.
#' @param k number of orthogonal time dependent covariates
#'
#' @return
#' \item{homo_coef}{List Individual's parameters' estimation}
#' \item{homo_cov}{List Individual's parameters' covairance matrix}
#' \item{nonhomo_coef}{List nonhomogeneous parameters' estimation}
#' \item{nonhomo_sd}{List nonhomogeneous parameter's standard deviation}
#'
#' @export

parametric <- function(data, zij, n, p, k = 0) {

  k1 = k
  dim_check = dim(zij)
  n1 = ncol(data)
  t = nrow(data)

  if (n1 < 3 || t < 1) {
    stop("Please check the format of data")
  }

  if (n != dim_check[1] || n != dim_check[2] || p != dim_check[3]) {
    stop("Please check whether the format of zij is correct or entered n, p is correct")
  }

  if (min(data[,1:2]) != 1 || max(data[,1:2]) != n) {
    stop("Please check whether the format of data is correct or entered n is correct")
  }

  if (!is.numeric(sum(data[,n1]))) {
    stop("Please check whether the time stamp is of numeric format")
  }

  t_max <- (data[t, n1] - data[1, n1])
  t_norm <- (data[, n1] - data[1, n1][[1]]) / t_max[[1]]
  t_sep_t = seq(0.01, 1, 0.01)

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

  F_poly = function(k1, x, n) {

    Pns = function(x, n) {
      2^0.5*sin(2*pi*n*x)
    }

    Pnc = function(x, n) {
      2^0.5*cos(2*pi*n*x)
    }

    if (n == 1) {
      return(1)
    } else if (n <= k1 + 1) {
      return(Pns(x, n-1))
    } else {
      return(Pnc(x, n-k1-1))
    }
  }

  n2 = k1*2+1
  N_t = matrix(0, nrow = n*(n-1), ncol = k1*2+1)
  b1k = matrix(0, nrow = (2*n-2), ncol = k1*2+1)
  sd1k = matrix(0, nrow = 2*n*n2, ncol = 2*n*n2)
  theta = matrix(0, nrow = p, ncol = 1)
  theta_sd = matrix(0, nrow = p, ncol = 1)

  b1 = matrix(0, nrow = 2, ncol = k1*2+1)

  temp0 = xConstruct(n, zij)
  temp = temp0$x
  temp1 = temp0$zij
  P1a = solve(t(temp) %*% temp) %*% t(temp)
  P2a = solve(t(temp1) %*% temp1) %*% t(temp1)

  # Homogeneous Parameter Estimation
  for (i in seq_len(t)) {
    z1 = data[i, 1]
    z2 = data[i, 2]
    z3 = hash[(z1-1)*n+z2]

    for (i1 in seq_len(n2)) {
      N_t[z3, i1] = N_t[z3, i1] + F_poly(k1, t_norm[i], i1)
      if (t_norm[i] >= t_sep_t[100]) {
        b1k[, i1] = P1a %*% matrix(N_t[, i1])
        b1[1, i1] = -sum(b1k[1:(n-1), i1])
        b1[2, i1] = -sum(b1k[n:(2*n-2), i1])
      }
    }

    # Covariance Matrix Estimation
    for (j in seq_len(n2)) {
      for (z in seq_len(j)) {
        p1 = (z1-1)*n2
        q1 = (z2-1)*n2+n*n2

        sd1k[p1+j, p1+z] = sd1k[p1+j, p1+z] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        sd1k[p1+z, p1+j] = sd1k[p1+z, p1+j] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        sd1k[q1+j, q1+z] = sd1k[q1+j, q1+z] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        sd1k[q1+z, q1+j] = sd1k[q1+z, q1+j] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
      }
    }

    # Non-Homogeneous Parameters' Estimation and standard deviation
    if (t_norm[i] >= t_sep_t[100]) {
      theta = 1 / t_norm[i] * P2a %*% matrix(N_t[, 1])
      for (j in 1:p) {
        theta_sd[j, 1] = sum(P2a[j,]^2*N_t[, 1])
      }
    }
  }

  sd1k = sd1k / n

  output <- list(homo_coef = list(outgoing = rbind(matrix(b1[1,], nrow = 1),
                                                   as.matrix(b1k[1:(n-1),])),
                                  incoming = rbind(matrix(b1[2,], nrow = 1),
                                                   as.matrix(b1k[n:(2*n-2),]))),
                 homo_cov = sd1k,
                 nonhomo_coef = theta,
                 nonhomo_sd = theta_sd^0.5)

  output
}
