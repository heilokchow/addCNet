npara = proxyData$modelfit
trail = proxyData$trail
zij = proxyData$zij
n = proxyData$n
p = dim(proxyData$zij)[[3]]
temp_date = proxyData$tz
sect = proxyData$sector


# Refit Data --------------------------------------------------------------

npara = nonParametric(trail, zij, n, p, tz = temp_date, h1 = 0.05)

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



# motivation --------------------------------------------------------------

(as.Date("2009/06/25") - as.Date("2008/09/08")) * 0.3294 + as.Date("2008/09/08")
(as.Date("2009/06/25") - as.Date("2008/09/08")) * 0.6123 + as.Date("2008/09/08")

library(scales)
library(igraph)

Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nrow(trail)) {
  if (trail[z, 3] > 0.6123) {
    break
  }
  p = trail[z,1]
  q = trail[z,2]
  Nij[p, q] = Nij[p, q] + 1
}

rownames(Nij) = 1:n
colnames(Nij) = 1:n

colrs <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
yg = unique(sect$year_school)

cont_plot = sect$year_school
for (i in 1:5) {
  cont_plot[which(cont_plot == yg[i])] = i
}

cont_plot[which(cont_plot == "")] = 1
arrange = order(cont_plot)

net2 <- graph_from_adjacency_matrix(zij[arrange, arrange, 1, 3], weighted = TRUE, mode = "directed")
cont_plot = as.numeric(cont_plot)

cont_plot1 = cont_plot[arrange]
vsize = (log(rowSums(Nij) + colSums(Nij)+1) + 1)[arrange]
vsize = (sqrt(rowSums(Nij) + colSums(Nij)+1)/10)[arrange]

V(net2)$color = colrs[cont_plot1]
V(net2)$size = 15*(vsize - min(vsize))/ (max(vsize) - min(vsize)) + 2
V(net2)$label <- NA

plot(net2, layout=layout.circle, edge.width = log(E(net2)$weight), edge.color = "gray80", edge.arrow.size=.5)


k = seq(5, 95, 5)

samefloor = matrix(NA, nrow = 2, ncol = length(k))
t_sep_t = seq(1, 100, 1)/100

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

  samefloor[1, z/5] = sum(Nij[zij[,,6,1] == 0])
  samefloor[2, z/5] = sum(Nij[zij[,,6,1] == 1])

  cat(z, "\n")
}


plot_t1 = seq(as.Date("2008/09/08"),as.Date("2009/06/25"),length = 19)
df1 = data.frame(x = rep(plot_t1, 2), y = c(samefloor[1,], samefloor[2,]), group = c(rep("Not same floor", 19), rep("Same floor", 19)))

ggplot(df1, aes(x = x, y = y, group = group, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "")  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("Time") +
  ylab("No. communication")



Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nrow(trail)) {
  if (trail[z, 3] > 1) {
    break
  }
  p = trail[z,1]
  q = trail[z,2]
  Nij[p, q] = Nij[p, q] + 1
}

rownames(Nij) = 1:n
colnames(Nij) = 1:n

colrs <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
yg = unique(sect$year_school)

cont_plot = sect$year_school
cont_plot[which(cont_plot == "")] = "GRT / Other"
arrange = order(rowSums(Nij) + colSums(Nij), decreasing = TRUE)

df2 = data.frame(x = seq(1, n, 1), y = (rowSums(Nij) + colSums(Nij))[arrange], group = cont_plot[arrange])

ggplot(df2, aes(x = x, y = y, group = group, color = group)) +
  geom_point() +
  scale_color_discrete(name = "")  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("Individual") +
  ylab("No. communication")

# Group / Ind Het plot ----------------------------------------------------


plot_t1 = seq(as.Date("2008/09/08"),as.Date("2009/06/25"),length = 100)

# Baseline

Nij = matrix(0, nrow = n, ncol = n)
sdB = rep(0, 100)
co = 1
for (i in 1:nrow(trail)) {

  p = trail[i,1]
  q = trail[i,2]
  t = trail[i,3]
  Nij[p, q] = Nij[p, q] + 1

  if (t >= t_sep_t[co]) {
    sdB[co] =  sqrt(sum(Nij) / n^4)
    co = co + 1
  }
}

bsestimate = npara$homo_coefficients$baseline

df = data.frame(x = plot_t1, y = c(bsestimate),
                yl = c(bsestimate - 1.96 * sdB),
                yu = c(bsestimate + 1.96 * sdB))

ggplot(df, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2008/12/13"), linetype="dotted") +
  geom_vline(xintercept = as.Date("2009/01/05"), linetype="dotted") +
  geom_vline(xintercept = as.Date("2009/03/20"), linetype="dotted") +
  geom_vline(xintercept = as.Date("2009/03/30"), linetype="dotted") +
  geom_vline(xintercept = as.Date("2009/05/30"), linetype="dotted") +
  xlab(expression(italic("t"))) +
  ylab("Cumulative baseline hazard")  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")

# Group by floor

uni_sect = unique(sect[,3])
uni_sect = uni_sect[-which(uni_sect == "")]


# alp

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$floor == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$ab$outgoing[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$ab$sdout[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$ab$outgoing[floorid, ]
    floorsd = sqrt(npara$ab$sdout[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\alpha}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# A

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$floor == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$homo_coefficients$outgoing[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$homo_coefficients$sdout[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$homo_coefficients$outgoing[floorid, ]
    floorsd = sqrt(npara$homo_coefficients$sdout[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{A}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# bet

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$floor == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$ab$incoming[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$ab$sdinc[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$ab$incoming[floorid, ]
    floorsd = sqrt(npara$ab$sdinc[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\beta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# B

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$floor == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$homo_coefficients$incoming[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$homo_coefficients$sdinc[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$homo_coefficients$incoming[floorid, ]
    floorsd = sqrt(npara$homo_coefficients$sdinc[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{B}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# Group by Year

uni_sect = unique(sect[,2])
uni_sect = uni_sect[-which(uni_sect == "")]


# alp

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$year_school == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$ab$outgoing[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$ab$sdout[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$ab$outgoing[floorid, ]
    floorsd = sqrt(npara$ab$sdout[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\alpha}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# A

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$year_school == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$homo_coefficients$outgoing[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$homo_coefficients$sdout[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$homo_coefficients$outgoing[floorid, ]
    floorsd = sqrt(npara$homo_coefficients$sdout[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{A}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# bet

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$year_school == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$ab$incoming[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$ab$sdinc[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$ab$incoming[floorid, ]
    floorsd = sqrt(npara$ab$sdinc[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\beta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# B

floor = data.frame()

for (i in seq_len(length(uni_sect))) {
  floorid = which(sect$year_school == uni_sect[i])
  if (length(floorid) > 1) {
    floormean = apply(npara$homo_coefficients$incoming[floorid, ], 2, mean)
    floorsd = sqrt(apply(npara$homo_coefficients$sdinc[floorid, ], 2, mean) / length(floorid))
  } else {
    floormean = npara$homo_coefficients$incoming[floorid, ]
    floorsd = sqrt(npara$homo_coefficients$sdinc[floorid, ])
  }
  floordf = data.frame(x = plot_t1,
                       y = floormean,
                       yl = floormean - 1.96 * floorsd,
                       yu = floormean + 1.96 * floorsd,
                       group = rep(uni_sect[i], 100))
  floor = rbind(floor, floordf)

}


ggplot(floor, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{B}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# Thetas

tk = seq(1, 100, 1)
PaT = npara$Pa[(2*n):(2*n+10-1),]

rd_theta_sd = matrix(0, nrow = 10, ncol = 100)
rd_theta_ab_sd = matrix(0, nrow = 10, ncol = 100)

for (j in 1:100) {
  NTs = npara$NTs[, tk[j]]
  rd_theta_ab_sd[,j] = diag(PaT %*% diag(NTs) %*% t(PaT))


  NT = npara$NT[, tk[j]]
  rd_theta_sd[,j] = diag(PaT %*% diag(NT) %*% t(PaT))

  cat(j, "\n")
}

groupName = c("Close friend", "Socialize twice per week", "Political discussant",
              "Facebook all tagged photos", "Blog livejournal twitter", "Same floor",
              "Same year", "1 year difference", "2 year difference", "3 year difference")

# tht

tht = data.frame()
for (i in seq_len(10)) {
  thtdf = data.frame(x = plot_t1,
                     y = npara$th[i,],
                     yl = npara$th[i,] - 1.96 * sqrt(rd_theta_ab_sd[i,]),
                     yu = npara$th[i,] + 1.96 * sqrt(rd_theta_ab_sd[i,]),
                     group = rep(groupName[i], 100))
  tht = rbind(tht, thtdf)

}

ggplot(tht, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\theta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# T

tht = data.frame()
for (i in seq_len(10)) {
  thtdf = data.frame(x = plot_t1,
                     y = npara$nonhomo_coefficients[i,],
                     yl = npara$nonhomo_coefficients[i,] - 1.96 * sqrt(rd_theta_sd[i,]),
                     yu = npara$nonhomo_coefficients[i,] + 1.96 * sqrt(rd_theta_sd[i,]),
                     group = rep(groupName[i], 100))
  tht = rbind(tht, thtdf)

}

ggplot(tht, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\Theta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# Tests -------------------------------------------------------------------

floor_comp = matrix(0, nrow = 8, ncol = 8)
for (j in seq_len(length(uni_sect))) {
  for (i in seq_len(length(uni_sect))) {
    floorid1 = which(sect$floor == uni_sect[j])
    floorid2 = which(sect$floor == uni_sect[i])
    floorid = c(floorid1, floorid2)

    floor_comp[j, i] = degreeHTest(npara$ab$outgoing[floorid, ], npara$ab$sdinc[floorid, ], 25)$tstat

  }
}

