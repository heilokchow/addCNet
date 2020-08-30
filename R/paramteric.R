#' Semi Parametric Estimate
#'
#' @description
#' This function allows you to do non parametric estimation for additive network model.
#' For non-homofily effect, we only consider non-time dependent covariates since it is
#' the most common case in real world while time dependent cases are rare and difficult
#' to be specified. Also, since the projection should be done at any time t, the model
#' with time dependent covariates requires much higher computation power.
#'
#' @param data The whole dataset containing pairwise communications, the first column
#' is the person with outgoing communication and the second column is the person with
#' incoming communication. The last column is the communication time point (should be
#' in numeric format).
#' @param zij The transfromed non-homogeneous effect (array with dimension (n,n,p)). The
#' transformation can be done using zProjection() function.
#' @param n number of individuals in the whole network.
#' @param p number of non-homogeneuous effect.
#' @param tz Time points for non-homofily effects (ranging from (0,1)), the length of this
#' vector should equal to the 4th dimension of zij which reflect the non-homofily effect
#' on or before those specific time points.
#'
#' @param k number of orthogonal time dependent covariates
#'
#' @return
#' \item{homo_coef}{List Individual's parameters' estimation}
#' \item{homo_cov}{List Individual's parameters' covairance matrix}
#' \item{nonhomo_coef}{List nonhomogeneous parameters' estimation}
#' \item{nonhomo_cov}{List nonhomogeneous parameter's covariance estimation}
#'
#' @export

parametric <- function(data, zij, n, p, tz, k = 0) {

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

  F_poly_norm = function(k1, x, n) {

    Pns = function(x, n) {
      x - sin(4*pi*n*x) / (4*pi*n)
    }

    Pnc = function(x, n) {
      x + sin(4*pi*n*x) / (4*pi*n)
    }

    if (n == 1) {
      return(x)
    } else if (n <= k1 + 1) {
      return(Pns(x, n-1))
    } else {
      return(Pnc(x, n-k1-1))
    }
  }

  n2 = k1*2+1
  N_t = matrix(0, nrow = n*(n-1), ncol = n2)
  b1k = matrix(0, nrow = (2*n-2), ncol = n2)
  sd1k = matrix(0, nrow = 2*n*n2, ncol = 2*n*n2)

  b1 = matrix(0, nrow = 2, ncol = k1*2+1)

  temp0 = xConstruct(n, zij)
  temp = temp0$x
  P1a = solve(t(temp) %*% temp) %*% t(temp)

  N_tz1 = matrix(0, nrow = n*(n-1), ncol = n2)
  N_tz2 = matrix(0, nrow = n*(n-1), ncol = n2^2)

  count1 = array(0, c(p*n2, p*n2))
  count2 = matrix(0, nrow = p*n2, ncol = 1)
  count3 = array(0, c(p*n2, p*n2))

  theta = matrix(0, nrow = p, ncol = n2)
  theta_cov = array(0, c(p*n2, p*n2))

  t0 = 0
  co_t = 1
  t1 = tz[co_t]

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

    if (t_norm[i] < t1) {
      for (i1 in seq_len(n2)) {
        N_tz1[z3, i1] = N_tz1[z3, i1] + F_poly(k1, t_norm[i], i1)
        for (j1 in i1:n2) {
          idx = (i1-1)*n2+j1
          N_tz2[z3, idx] = N_tz2[z3, idx] + F_poly(k1, t_norm[i], i1) * F_poly(k1, t_norm[i], j1)
        }
      }
    } else {
      count2_temp = matrix(NA, nrow = 0, ncol = 1)
      temp = temp0$zij[,,co_t]
      temp1 = t(temp) %*% temp
      for (i1 in seq_len(n2)) {
        count2_temp = rbind(count2_temp, t(temp) %*% N_tz1[, i1])
        for (j1 in i1:n2) {
          idx0 = (i1-1)*n2+j1
          idx1 = (i1-1)*p+1
          idx2 = (j1-1)*p+1

          F_poly_cov = function(x) {
            if (i1 == 1 && j1 == 1) {
              return(x^0)
            } else {
              return(F_poly(k1, x, i1) * F_poly(k1, x, j1))
            }
          }

          if (i1 == j1) {
            count1[idx1:(idx1+p-1), idx2:(idx2+p-1)] = count1[idx1:(idx1+p-1), idx2:(idx2+p-1)] +
              temp1 * integrate(F_poly_cov, t0, t1)$value
            count3[idx1:(idx1+p-1), idx2:(idx2+p-1)] = count3[idx1:(idx1+p-1), idx2:(idx2+p-1)] +
              t(temp) %*% diag(N_tz2[, idx]) %*% temp
          } else {
            count1[idx1:(idx1+p-1), idx2:(idx2+p-1)] = count1[idx1:(idx1+p-1), idx2:(idx2+p-1)] +
              temp1 * integrate(F_poly_cov, t0, t1)$value
            count1[idx2:(idx2+p-1), idx1:(idx1+p-1)] = count1[idx2:(idx2+p-1), idx1:(idx1+p-1)] +
              temp1 * integrate(F_poly_cov, t0, t1)$value

            count3[idx1:(idx1+p-1), idx2:(idx2+p-1)] = count3[idx1:(idx1+p-1), idx2:(idx2+p-1)] +
              t(temp) %*% diag(N_tz2[, idx]) %*% temp
            count3[idx2:(idx2+p-1), idx1:(idx1+p-1)] = count3[idx2:(idx2+p-1), idx1:(idx1+p-1)] +
              t(temp) %*% diag(N_tz2[, idx]) %*% temp
          }
        }
      }
      count2[, ] = count2[, ] + count2_temp
      t0 = t1
      if (co_t < length(tz)) {
        co_t = co_t + 1
        t1 = tz[co_t]
      }
      N_tz1 = matrix(0, nrow = n*(n-1), ncol = n2)
      N_tz2 = matrix(0, nrow = n*(n-1), ncol = n2^2)
    }

    # Covariance Matrix Estimation
    for (j in seq_len(n2)) {
      for (z in seq_len(j)) {
        p1 = (z1-1)*n2
        q1 = (z2-1)*n2+n*n2

        sd1k[p1+j, p1+z] = sd1k[p1+j, p1+z] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        sd1k[q1+j, q1+z] = sd1k[q1+j, q1+z] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        if (j != z) {
          sd1k[p1+z, p1+j] = sd1k[p1+z, p1+j] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
          sd1k[q1+z, q1+j] = sd1k[q1+z, q1+j] + F_poly(k1, t_norm[i], j)*F_poly(k1, t_norm[i], z)
        }
      }
    }

    # Non-Homogeneous Parameters' Estimation and standard deviation
    if (t_norm[i] >= t_sep_t[100]) {
      inv_count = solve(count1, count2)
      for (z in seq_len(n2)) {
        theta[, z] = inv_count[((z-1)*p+1):(z*p),1]
      }
      inv_count1 = solve(count1)
      theta_cov = inv_count1 %*% count3 %*% inv_count1
    }
  }

  sd1k = sd1k / n^2

  output <- list(homo_coef = list(outgoing = rbind(matrix(b1[1,], nrow = 1),
                                                   as.matrix(b1k[1:(n-1),])),
                                  incoming = rbind(matrix(b1[2,], nrow = 1),
                                                   as.matrix(b1k[n:(2*n-2),]))),
                 homo_cov = sd1k,
                 nonhomo_coef = theta,
                 nonhomo_cov = theta_cov,
                 check = list(x = count1, y = count3))

  output
}
