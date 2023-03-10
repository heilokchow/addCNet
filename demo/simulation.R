library(addCNet)
library(ggplot2)
library(latex2exp)

# Functions ---------------------------------------------------------------

bs <- function(t) {1}
fs <- function(t, i) {0.04*(i-(n+1)/2) + 0.04*(i-(n+1)/2)*sin(2*pi*t)} # cautious about negative hazard
fr <- function(t, i) {0.04*(i-(n+1)/2) + 0.04*(i-(n+1)/2)*cos(2*pi*t)}
# fs <- function(t, i) {0.02*(i-(n+1)/2)}
# fr <- function(t, i) {0.02*(i-(n+1)/2)}
fg <- function(t, k) {0.2*k}

fsI <- function(t, q) {
  if (t <= 0) {
    return(0)
  }
  integrate(fs, 0, t, i = q)[[1]]
}

fsI1 = Vectorize(fsI)

frI <- function(t, q) {
  if (t <= 0) {
    return(0)
  }
  integrate(fr, 0, t, i = q)[[1]]
}

frI1 = Vectorize(frI)

fgI <- function(t, q) {
  if (t <= 0) {
    return(0)
  }
  integrate(fg, 0, t, k = q)[[1]]
}

fgI1 = Vectorize(fgI)


# Simulation begin --------------------------------------------------------


n = 30
p = 1

t_sep_t = seq(0.01, 1, 0.01)

trueBo = matrix(0, nrow = n, ncol = 100)
trueBi = matrix(0, nrow = n, ncol = 100)

for (i in 1:n) {
  for (j in 1:100) {
    trueBo[i, j] = fsI1(t_sep_t[j], i)
    trueBi[i, j] = frI1(t_sep_t[j], i)
  }
}

# Constructing Non-homofily effect
zij <- array(rnorm(n*n), c(n, n, 1, 1))

all_result = list()
all_result_zval = list()
all_result_kh = list()

set.seed(1)
for (i in 1:100) {
  cat(i, '\n')
  result <- tGenerate(bs, fs, fr, fg, n, 1, array(zij, c(n,n,p,1)), tz = 1, maxit = 20)
  np0 = nonParametric(result$trail, array(zij, c(n,n,p,1)), n, p, h1 = 0.05)

  all_result[[i]] = rbind(np0$homo_coefficients$outgoing, np0$homo_coefficients$incoming, np0$nonhomo_coefficients)
  all_result_zval[[i]] = rbind((np0$homo_coefficients$outgoing - trueBo) / sqrt(np0$homo_coefficients$sdout),
                               (np0$homo_coefficients$incoming - trueBi) / sqrt(np0$homo_coefficients$sdinc))
  all_result_kh[[i]] = rbind(np0$ab$outgoing, np0$ab$incoming, np0$th)
}

all_mean = array(0, c(2*n + p, 100, 100))
all_mean_kh = array(0, c(2*n + p, 100, 100))
for (i in 1:100) {
  all_mean[,,i] = all_result[[i]]
  all_mean_kh[,,i] = all_result_kh[[i]]
}

resMean = apply(all_mean, c(1,2), mean)
reskhMean = apply(all_mean_kh, c(1,2), mean)
resl = apply(all_mean, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
resu = apply(all_mean, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
rekhsl = apply(all_mean_kh, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
rekhsu = apply(all_mean_kh, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)


# Plot --------------------------------------------------------------------


q = 30
p1 = data.frame(x = seq(0.01,1,0.01), y = resMean[q,], yl = resl[q,], yu = resu[q,])
ggplot(p1, aes(x = x, y = y)) + geom_line() + stat_function(fun = fsI1, args = list(q = q), color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{A}_{n}$(t)')) +
  theme_bw()

p1 = data.frame(x = seq(0.1,0.9,0.01), y = reskhMean[q,10:90], yl = rekhsl[q,10:90], yu = rekhsu[q,10:90])
ggplot(p1, aes(x = x, y = y)) + geom_line() + stat_function(fun = fs, args = list(i = q), color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\alpha}_{n}$(t)')) +
  theme_bw()


q = 60
p2 = data.frame(x = seq(0.01,1,0.01), y = resMean[q,], yl = resl[q,], yu = resu[q,])
ggplot(p2, aes(x = x, y = y)) + geom_line() + stat_function(fun = frI1, args = list(q = q-n), color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{B}_{n}$(t)')) +
  theme_bw()

p2 = data.frame(x = seq(0.1,0.9,0.01), y = reskhMean[q,10:90], yl = rekhsl[q,10:90], yu = rekhsu[q,10:90])
ggplot(p2, aes(x = x, y = y)) + geom_line() + stat_function(fun = fr, args = list(i = q-n), color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\beta}_{n}$(t)')) +
  theme_bw()


q = 1
p3 = data.frame(x = seq(0.01,1,0.01), y = resMean[2*n+q,], yl = resl[2*n+q,], yu = resu[2*n+q,])
ggplot(p3, aes(x = x, y = y)) + geom_line() + stat_function(fun = function(x) {0.2*x}, color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\Theta}_{1}$(t)')) +
  theme_bw()

p3 = data.frame(x = seq(0.1,0.9,0.01), y = reskhMean[2*n+q,10:90], yl = rekhsl[2*n+q,10:90], yu = rekhsu[2*n+q,10:90])
ggplot(p3, aes(x = x, y = y)) + geom_line() + stat_function(fun = function(x) {0.2}, color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\theta}_{1}$(t)')) +
  theme_bw()

zon = c()
zin = c()

for (i in 1:100) {
  zon = c(zon, all_result_zval[[i]][n,][which(abs(all_result_zval[[i]][n,]) < Inf)])
  zin = c(zin, all_result_zval[[i]][2*n,][which(abs(all_result_zval[[i]][2*n,]) < Inf)])
}

hist(zon, seq(-7, 7, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\A}_n$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

hist(zin, seq(-7, 7, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\B}_n$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)
