#' Non Parametric Estimate
#' 
#' @description 
#' This function allows you to do non parametric estimation for additive network model
#' 
#' @param data The whole dataset containing pairwise communications, the first column
#' is the person with outgoing communication and the second column is the person with
#' incoming communication. The last column is the communication time point (should be 
#' in numeric format).
#' @param zij The transfromed non-homogeneous effect (array with dimension (n,n,p)). The 
#' transformation can be done using z_projection() function.
#' @param n number of individuals in the whole network.
#' @param p number of non-homogeneuous effect.
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

non_parametric <- function(data, zij, n, p) {
  
  dim_check = dim(zij)
  n1 = ncol(data)
  t = nrow(data)
  
  if (n1 < 3 || t < 1) {
    stop("Please check the format of data")
  }

  if (n != dim_check[1] || n != dim_check[2] || p != dim_check[3]) {
    stop("Please check whether the format of zij is correct or entered n, p is correct")
  }
  
  if (min(data[,1:2] != 1 || max(data[,1:2] != n))) {
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
  
  N_t = matrix(0, nrow = n*(n-1), ncol = 1)
  B = matrix(0, nrow = (2*n-1+p), ncol = 100)
  SD = matrix(0, nrow = (2*n-1+p), ncol = 100)
  
  B1 = matrix(0, nrow = 2, ncol = 100)
  SD1 = matrix(0, nrow = 2, ncol = 100)
  N1 = matrix(0, nrow = 2, ncol = 1)
  
  co = 1
  for (i in seq_len(t)) {
    z1 = whole_time_enron[i, 1]
    z2 = whole_time_enron[i, 2]
    z3 = hash[(z1-1)*n+z2]
    N_t[z3, 1] = N_t[z3, 1] + 1
    
    if (z1 == 1) {
      N1[1, 1] = N1[1, 1] + 1
    }
    
    if (z2 == 1) {
      N1[2, 1] = N1[2, 1] + 1
    }
    
    if (t_norm[i] >= t_sep_t[co]) {
      B[, co] = Pa %*% N_t
      B1[1, co] = -sum(B[2:n, co])
      B1[2, co] = -sum(B[(n+1):(2*n-1), co])
      
      for (j in 1:(2*n-1+p)) {
        SD[j, co] = sum(Pa[j,]^2*N_t)
      }
      SD1[1, co] = N1[1, 1] / n^2
      SD2[1, co] = N1[2, 1] / n^2
      
      co = co + 1
      if (co == 101) {
        break
      }
    }
  }
  
  output <- list(homo_coefficients = list(baseline = matrix(B[1,], nrow = 1),
                                outgoing = cbind(matrix(B1[1,], nrow = 1),
                                                 B[2:n,]),
                                incoming = cbind(matrix(B1[2,], nrow = 1),
                                                 B[(n+1):(2*n-1),])),
       homo_sd = list(baseline = matrix(SD[1,], nrow = 1),
                      outgoing = cbind(matrix(SD1[1,], nrow = 1)^0.5,
                                       SD[2:n,]),
                      incoming = cbind(matrix(SD1[2,], nrow = 1)^0.5,
                                       SD[(n+1):(2*n-1),])),
       nonhomo_coefficients = B[(2*n):(2*n-1+p),],
       nonhomo_sd = SD[(2*n):(2*n-1+p),]^0.5)
  
  output
}
