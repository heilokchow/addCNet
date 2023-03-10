npara = proxyData$modelfit
trail = proxyData$trail
zij = proxyData$zij
n = proxyData$n
temp_date = proxyData$tz


# Correlation plot --------------------------------------------------------

k = seq(5, 95, 5)

corr1 = rep(0, length(k))
corrZ = matrix(NA, nrow = 10, ncol = length(k))

for (z in k) {

  ts = z
  te = z + 5

  Nij = matrix(0, nrow = n, ncol = n)

  for (i in 1:nrow(trail)) {
    if (trail[i,3] > t_sep_t[ts] && trail[i,3] <  t_sep_t[te] ) {
      p = trail[i,1]
      q = trail[i,2]
      Nij[p, q] = Nij[p, q] + 1
    }
  }

  for (z1 in 1:10) {

    Lij = matrix(0, nrow = n, ncol = n)

    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          k2 = which.max(z/100 < temp_date)
          Lij[i, j] = Lij[i, j] + npara$nonhomo_coefficients[z1, te] * zij[i,j,z1,k2]
          Lij[i, j] = Lij[i, j] - npara$nonhomo_coefficients[z1, ts] * zij[i,j,z1,k2]
        }
      }
    }

    corrZ[z1, z/5] = cor(c(Nij), c(Lij))
  }

  Lij = matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        Lij[i, j] = npara$homo_coefficients$outgoing[i, te] + npara$homo_coefficients$incoming[j, te]
        Lij[i, j] = Lij[i, j] - npara$homo_coefficients$outgoing[i, ts] - npara$homo_coefficients$incoming[j, ts]
      }
    }
  }

  corr1[z/5] = cor(c(Nij), c(Lij))

  cat(z, "\n")
}

corrAB = corr1
corrSy = colSums(corrZ[1:5,])
corrSf = corrZ[6,]
corrsj = colSums(corrZ[7:10,])

plot_t1 = seq(as.Date("2008/09/08"),as.Date("2009/06/25"),length = 19)

df1 = data.frame(x = rep(plot_t1, 4), y = c(corrAB, corrSy, corrSf, corrsj), group = c(rep("heterogeneous", 19),
                                                                                       rep("Survey data", 19),
                                                                                       rep("Same floor", 19),
                                                                                       rep("Same year", 19)))


ggplot(df1, aes(x = x, y = y, group = factor(group, levels = c("Survey data", "Same year", "Same floor", "heterogeneous")), fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  scale_y_continuous(limits = c(-0.05, 1)) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("Time") +
  ylab("Correlation")
