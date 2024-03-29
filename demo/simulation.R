library(addCNet)
library(ggplot2)
library(latex2exp)

# Functions ---------------------------------------------------------------

n = 30
p = 1

t_sep_t = seq(0.01, 1, 0.01)

# SIM1
bs <- function(t) {3}
fs <- function(t, i) {0.04*(i-(n+1)/2)*sin(2*pi*t)} # cautious about negative hazard
fr <- function(t, i) {0.04*(i-(n+1)/2)*cos(2*pi*t)}
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


# SIM1.5

n = 60
p = 1

bs <- function(t) {0.2}
fs <- function(t, i) {
  if (i <= n/3) return(0.1+0.1*sin(2*pi*t))
  return(-0.05-0.05*sin(2*pi*t))
} # cautious about negative hazard
fr <- function(t, i) {
  if (i <= n/3) return(0.1+0.1*cos(2*pi*t))
  return(-0.05-0.05*cos(2*pi*t))
}
fg <- function(t, k) {0.05+t*0}

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


# SIM1.6

n = 100
p = 1

bs <- function(t) {2/sqrt(n)}
fs <- function(t, i) {
  if (i <= n/2) return((0.5+0.5*sin(2*pi*t))/sqrt(n))
  return(-(0.5+0.5*sin(2*pi*t))/sqrt(n))
} # cautious about negative hazard
fr <- function(t, i) {
  if (i <= n/2) return((0.5+0.5*cos(2*pi*t))/sqrt(n))
  return(-(0.5+0.5*cos(2*pi*t))/sqrt(n))
}
fg <- function(t, k) {0.25/sqrt(n)+t*0}

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


# SIM2

n = 400
p = 1

bs <- function(t) {2}
fs <- function(t, i) {
  if (i < n/3) return(sin(2*pi*t)/2*0)
  if (i > 2*n/3) return(-sin(2*pi*t)/2*0)
  return(0)
}
fr <- function(t, i) {
  if (i < n/3) return(sin(2*pi*t)/2)
  if (i > 2*n/3) return(-sin(2*pi*t)/2)
  return(0)
}
fg <- function(t, k) {0.2*k}

fsI <- function(t, q) {
  if (t <= 0) {
    return(0)
  }
  integrate(Vectorize(fs), 0, t, i = q)[[1]]
}

fsI1 = Vectorize(fsI)

frI <- function(t, q) {
  if (t <= 0) {
    return(0)
  }
  integrate(Vectorize(fr), 0, t, i = q)[[1]]
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

trueBo = matrix(0, nrow = n, ncol = 100)
trueBi = matrix(0, nrow = n, ncol = 100)


for (i in 1:n) {
  for (j in 1:100) {
    trueBo[i, j] = fsI1(t_sep_t[j], i)
    trueBi[i, j] = frI1(t_sep_t[j], i)
  }
}
trueBoi = rbind(trueBo, trueBi)


truebo = matrix(0, nrow = n, ncol = 100)
truebi = matrix(0, nrow = n, ncol = 100)

for (i in 1:n) {
  for (j in 1:100) {
    truebo[i, j] = fs(t_sep_t[j], i)
    truebi[i, j] = fr(t_sep_t[j], i)
  }
}

trueboi = rbind(truebo, truebi)

# Constructing Non-homofily effect

all_result = list()
all_result_var = list()
all_result_zval = list()
all_result_kh = list()
all_result_kh_var = list()
all_result_kh_zval = list()

rep = 100
theta_sd = matrix(nrow = rep, ncol = 3)
# theta_sd = rbind(theta_sd, matrix(0, nrow = rep-100, ncol = 3))
theta_ab_sd = matrix(nrow = rep, ncol = 3)
# theta_ab_sd = rbind(theta_ab_sd, matrix(0, nrow = rep-100, ncol = 3))

set.seed(100)
zij <- array(runif(n*n), c(n, n, 1, 1))
for (i in 1:rep) {
  cat(i, '\n')
  set.seed(i)

  # Generate Data
  result <- tGenerateC(n, shift1 = 1, shift2 = 0, Zij = array(zij, c(n,n,p)))

  # Model Run
  t1 = Sys.time()
  np0 = nonParametric(result$trail, array(zij, c(n,n,p,1)), n, p, h1 = 0.07, test = 0, pvpval = 0)
  t2 = Sys.time()

  all_result[[i]] = rbind(np0$homo_coefficients$outgoing, np0$homo_coefficients$incoming, np0$nonhomo_coefficients)
  all_result_var[[i]] = rbind(np0$homo_coefficients$sdout, np0$homo_coefficients$sdinc)
  all_result_zval[[i]] = rbind((np0$homo_coefficients$outgoing - trueBo) / sqrt(np0$homo_coefficients$sdout),
                               (np0$homo_coefficients$incoming - trueBi) / sqrt(np0$homo_coefficients$sdinc))

  all_result_kh[[i]] = rbind(np0$ab$outgoing, np0$ab$incoming, np0$th)
  all_result_kh_var[[i]] = rbind(np0$ab$sdout, np0$ab$sdinc)
  all_result_kh_zval[[i]] = rbind((np0$ab$outgoing - truebo) / sqrt(np0$ab$sdout),
                                  (np0$ab$incoming - truebi) / sqrt(np0$ab$sdinc))

  # Additional output for tests

  tk = c(30, 50, 70)
  PaT = np0$Pa[2*n,]

  for (j in 1:3) {
    NTs = np0$NTs[, tk[j]]
    theta_ab_sd[i, j] = sum(NTs * PaT * PaT)

    NT = np0$NT[, tk[j]]
    theta_sd[i, j] = sum(NT * PaT * PaT)
  }
}

all_mean = array(0, c(2*n + p, 100, rep))
all_mean_sd = array(0, c(2*n, 100, rep))
all_mean_kh = array(0, c(2*n + p, 100, rep))
all_mean_sd_kh = array(0, c(2*n, 100, rep))
for (i in 1:rep) {
  all_mean[,,i] = all_result[[i]]
  all_mean_sd[,,i] = sqrt(all_result_var[[i]])
  all_mean_kh[,,i] = all_result_kh[[i]]
  all_mean_sd_kh[,,i] = sqrt(all_result_kh_var[[i]])
}

resMean = apply(all_mean, c(1,2), mean)
reskhMean = apply(all_mean_kh, c(1,2), mean)
resl = apply(all_mean, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
resu = apply(all_mean, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
rekhsl = apply(all_mean_kh, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
rekhsu = apply(all_mean_kh, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)


# Plot --------------------------------------------------------------------

Nij = matrix(0, n, n)
for (i in 1:nrow(result$trail)) {
  p1 = result$trail[i, 1]
  q1 = result$trail[i, 2]
  Nij[p1, q1] = Nij[p1, q1] + 1
}

t_seq = c(20, 40, 80, 160, 320, 640)
t_seq_1 = c(0.035, 0.109, 0.567, 3.329, 25.01, 270.62)
fit1 = lm(log(t_seq_1)~log(t_seq))

plot(log(t_seq), log(t_seq_1), ylab = "Log(Runtime in seconds)", xlab = "Log(Sample size)")
plot(t_seq, t_seq_1, ylab = "Log(Runtime in seconds)", xlab = "Log(Sample size)", log = c('xy'))
abline(a = -13.5, b = 3, col = "red")


# ABT plot ----------------------------------------------------------------

q = 10
q = n
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


q = n*2
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
ggplot(p3, aes(x = x, y = y)) + geom_line() + stat_function(fun = function(x) {0.25/sqrt(n)*x}, color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\Theta}_{1}$(t)')) +
  theme_bw()

p3 = data.frame(x = seq(0.1,0.9,0.01), y = reskhMean[2*n+q,10:90], yl = rekhsl[2*n+q,10:90], yu = rekhsu[2*n+q,10:90])
ggplot(p3, aes(x = x, y = y)) + geom_line() + stat_function(fun = function(x) {0.25/sqrt(n)}, color = "red") +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\theta}_{1}$(t)')) +
  theme_bw()

# Not use
p4 = data.frame(x = rep(seq(0.01,0.99,0.01), n),
                y = c(t(reskhMean[1:n,1:99])),
                id = sort(rep(seq(1,n,1), 99)),
                group = c(rep("1", 99*16), rep("2", 99*17), rep("3", 99*17)))

ggplot(p4, aes(x = x, y = y, group = id)) + geom_line()

# Not use
p5 = data.frame(x = rep(seq(0.01,0.99,0.01), 3),
                yl = c(apply(rekhsl[1:15,1:99], 2, mean),
                       apply(rekhsl[16:33,1:99], 2, mean),
                       apply(rekhsl[34:50,1:99], 2, mean)),
                yu = c(apply(rekhsu[1:15,1:99], 2, mean),
                       apply(rekhsu[16:33,1:99], 2, mean),
                       apply(rekhsu[34:50,1:99], 2, mean)),
                group = c(rep("1", 99), rep("2", 99), rep("3", 99)))

ggplot(p5, aes(x = x, group = group, fill = group)) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "", values=c("#619CFF", "red", "green")) +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\alpha}$(t)'))

p6 = data.frame(x = rep(seq(0.01,0.99,0.01), 3),
                yl = c(apply(reskhMean[1:(n/3-1),1:99], 2, min),
                       apply(reskhMean[(n/3):(n/3*2),1:99], 2, min),
                       apply(reskhMean[(n/3*2+1):n,1:99], 2, min)),
                yu = c(apply(reskhMean[1:(n/3-1),1:99], 2, max),
                       apply(reskhMean[(n/3):(n/3*2),1:99], 2, max),
                       apply(reskhMean[(n/3*2+1):n,1:99], 2, max)),
                group = c(rep("1", 99), rep("2", 99), rep("3", 99)))

ggplot(p6, aes(x = x, group = group, fill = group)) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "", values=c("#619CFF", "red", "green")) +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\alpha}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


p7 = data.frame(x = rep(seq(0.01,0.99,0.01), 3),
                yl = c(apply(reskhMean[n+1:(n/3-1),1:99], 2, min),
                       apply(reskhMean[n+(n/3):(n/3*2),1:99], 2, min),
                       apply(reskhMean[n+(n/3*2+1):n,1:99], 2, min)),
                yu = c(apply(reskhMean[n+1:(n/3-1),1:99], 2, max),
                       apply(reskhMean[n+(n/3):(n/3*2),1:99], 2, max),
                       apply(reskhMean[n+(n/3*2+1):n,1:99], 2, max)),
                group = c(rep("1", 99), rep("2", 99), rep("3", 99)))

ggplot(p7, aes(x = x, group = group, fill = group)) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5, linetype = "dashed") +
  scale_fill_manual(name = "", values=c("#619CFF", "red", "green")) +
  xlab(expression(italic("t"))) +
  ylab(TeX('$\\widehat{\\beta}$(t)')) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# Coverage frequency ------------------------------------------------------


pk = c(n/3, 2*n/3, n+n/3, n+2*n/3)
pk = c(n/2, n, n+n/2, n+n)
tk = c(30, 50, 70)

e1 = 0
e1t = 0
for (i in 1:rep) {
  temp = abs(all_result_kh[[i]][1:(2*n), ] - trueboi)/1.96
  e1 = e1 + (temp[pk, tk] > sqrt(all_result_kh_var[[i]][pk, tk]))

  temp1 = abs(all_result_kh[[i]][2*n+p, ] - 0.25/sqrt(n))/1.96
  e1t = e1t + (temp1[tk] > sqrt(theta_ab_sd[i,]))
}
e1 = e1 / rep
e1t = e1t / rep
1-e1
1-e1t

xkk_sum_sd = apply(all_mean_sd_kh, c(1, 2), mean, na.rm = TRUE)
xkk_sum_sd_th = apply(theta_ab_sd, 2, mean, na.rm = TRUE)
round(xkk_sum_sd[pk, tk]*1.96*2, 2)
round(sqrt(xkk_sum_sd_th)*1.96*2, 2)


e2 = 0
e2t = 0

for (i in 1:rep) {
  temp = abs(all_result[[i]][1:(2*n), ] - trueBoi)/1.96
  e2 = e2 + (temp[pk, tk] > sqrt(all_result_var[[i]][pk, tk]))

  temp1 = abs(all_result[[i]][2*n+p, ] - 0.25/sqrt(n)*t_sep_t)/1.96
  e2t = e2t + (temp1[tk] > sqrt(theta_sd[i,]))
}
e2 = e2 / rep
e2t = e2t / rep
1-e2
1-e2t

xkk_sum_sd = apply(all_mean_sd, c(1, 2), mean, na.rm = TRUE)
xkk_sum_sd_th = apply(theta_sd, 2, mean, na.rm = TRUE)
round(xkk_sum_sd[pk, tk]*1.96*2, 2)
round(sqrt(xkk_sum_sd_th)*1.96*2, 2)



# T statistic plot --------------------------------------------------------


zon = c()
zin = c()

for (i in 1:100) {
  zon = c(zon, all_result_zval[[i]][n,][which(abs(all_result_zval[[i]][n,]) < Inf)])
  zin = c(zin, all_result_zval[[i]][n*3/2,][which(abs(all_result_zval[[i]][n*3/2,]) < Inf)])
}
#
# for (i in 1:100) {
#   zon = c(zon, all_result_zval[[i]][n,][which(abs(all_result_zval[[i]][n,]) < Inf)])
#   zin = c(zin, all_result_zval[[i]][2*n,][which(abs(all_result_zval[[i]][2*n,]) < Inf)])
# }

hist(zon, seq(-7, 7, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\A}_{n/2}$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

hist(zin, seq(-7, 7, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\B}_{n/2}$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)


zons = c()
zins = c()

for (i in 1:100) {
  zons = c(zons, all_result_kh_zval[[i]][n/3, 30:70])
  zins = c(zins, all_result_kh_zval[[i]][n*4/3, 30:70])
}

# zons = zons[which(zons > -7)]
hist(zons, seq(-7, 7, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\alpha}_{n/2}$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)

# zins = zins[which(zins > -7)]
hist(zins, seq(-8, 8, 0.5), ylim = c(0, 0.45), freq = FALSE, main = "", xlab = TeX('Standardized $\\widehat{\\beta}_{n/2}$'))
curve(dnorm(x), xlab = "", ylab = "", add = T, lwd = 2.0)


# Test of Directed effect and degree heterogeneity ------------------------

ap = 0
sall = c()
sall1 = c()

k1 = c((1):(n/3-1), (n/3*2):(n))

for (i in 1:100) {
  # out = directTest(all_result_kh[[i]][1:n,], all_result_kh_var[[i]][1:n,],
  #                  all_result_kh[[i]][(n+1):(2*n),], all_result_kh_var[[i]][(n+1):(2*n),], seq(10,90,20))
  # out = degreeHTest(all_result_kh[[i]][(n+1):(n+n/3),], all_result_kh_var[[i]][(n+1):(n+n/3),], 90)
  out = degreeHTest(all_result_kh[[i]][k1,], all_result_kh_var[[i]][k1,], 25)
  # out = degreeHTest(all_result_kh[[i]][1:n,], all_result_kh_var[[i]][1:n,], 75)

  sall = c(sall, out[[1]])
  sall1 = c(sall1, sum(all_result_kh[[i]][(1:n),60]^2/all_result_kh_var[[i]][(1:n),60]-1)/sqrt(2*n))
  ap = ap + out[[2]]
}

hist(all_result_kh[[i]][1:n,50]/sqrt(all_result_kh_var[[i]][1:n,50]),
     seq(-5, 5, 0.5), ylim = c(0, 0.8), freq = FALSE, main = "", xlab = TeX('$T_d$'))
curve(dnorm, xlab = "", ylab = "", add = T, lwd = 2.0)

mean((all_result_kh[[i]][1:n,]/sqrt(all_result_kh_var[[i]][1:n,]))^2)

test = c()
for (i in 1:100) {
  test = c(test, mean(all_result_kh_zval[[i]][(n+1):(2*n),]^2))
}

# Under Null
g = function(x) {9*pnorm(x)^8*dnorm(x)}
g = function(x) {5*pnorm(x)^4*dnorm(x)}
h = Vectorize(g)

hist(sall, seq(-6, 10, 0.5), ylim = c(0, 0.8), freq = FALSE, main = "", xlab = TeX('$T_d$'))
hist(sall, seq(-2, 6, 0.5), ylim = c(0, 0.8), freq = FALSE, main = "", xlab = TeX('$T_{dh}$'))
curve(h, xlab = "", ylab = "", add = T, lwd = 2.0)
curve(dnorm, xlab = "", ylab = "", add = T, lwd = 2.0)

base <- list(null = list(estimate = all_result_kh, var = all_result_kh_var))

# TD 3-3-3 200
# 5 9 40 90 100
plot(seq(0,0.5,0.125), c(0.05, 0.09, 0.40, 0.90, 1), ylim = c(0, 1), ylab = "Size and Power", xlab = expression(italic(C[1])))
lines(seq(0,0.5,0.125), c(0.05, 0.09, 0.40, 0.90, 1))

# TDH 3-3-3 200
# 14 24 40 86 100
plot(seq(0,1,0.25), c(0.12, 0.24, 0.40, 0.86, 1), ylim = c(0, 1), ylab = "Size and Power", xlab = expression(italic(C[2])))
lines(seq(0,1,0.25), c(0.12, 0.24, 0.40, 0.86, 1))

# TDH three
# 0.03 0.92 1
#      0.03 0.96
#           0.08
