# Load Dataset ------------------------------------------------------------


library(ggplot2)
library(ggpubr)
library(purrr)
library(igraph)
library(addCNet)
library(latex2exp)

setwd("Z:/new data/Friends&Family/DataRelease")

sms = read.csv("SMSLog.csv.bz2")
call = read.csv("CallLog.csv.bz2")
friend = read.csv("SurveyFriendship.csv")
copule = read.csv("SurveyCouples.csv")



# Data Preprocessing ------------------------------------------------------

unique(sms[,1])
unique(friend[,2])


rm1 = which(substring(friend[,2], 1, 2) == "fa")
rm2 = which(substring(friend[,3], 1, 2) == "fa")
rmall = union(rm1, rm2)

friend1 = friend[-rmall,]

rm1 = which(substring(sms[,1], 1, 2) == "fa")
rm2 = which(substring(sms[,2], 1, 2) == "fa")
rmall = union(rm1, rm2)
rm3 = which(sms[,2] == "")
rmall = union(rmall, rm3)

sms1 = sms[-rmall,]

ids = union(unique(sms1[,1]), unique(sms1[,2]))

mapid = list()

for (i in 1:length(ids)) {
  mapid[[ids[i]]] = i
}
n = length(ids)

trail = sms1[,1:3]

nn = nrow(trail)

for (i in 1:nn) {
  p1 = trail[i, 1]
  q1 = trail[i, 2]

  if (sms1[i, 4] == "incoming") {
    trail[i, 1] = mapid[[q1]]
    trail[i, 2] = mapid[[p1]]
  } else {
    trail[i, 1] = mapid[[p1]]
    trail[i, 2] = mapid[[q1]]
  }
}


temp = as.POSIXct(sms1$local_time)
temp2 = temp
temp = as.numeric(temp)
temp1 = temp[order(temp)]
temp2 = temp2[order(temp)]
plot(temp1)
sms2 = sms1[order(temp),]

tz0 = as.POSIXct(unique(friend1[,5]))
tz0 = as.numeric(tz0)
tz = (tz0 - min(temp)) / (max(temp) - min(temp))
tz[4] = 1

trail = trail[order(temp),]
temp = temp1

trail$local_time = (temp-temp[1])/(temp[length(temp)]-temp[1])

zij = array(0, dim = c(n, n, 3, 4))

# Couple
# Have Child

copule[,1] = tolower(copule[,1])
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      p1 = which(copule[,1] == ids[i])
      q1 = which(copule[,1] == ids[j])
      if (copule[p1, 2] == copule[q1, 2]) {
        if (copule[p1, 5] == 1) {
          zij[i,j,2,] = 1
        } else {
          zij[i,j,1,] = 1
        }
      }
    }
  }
}


# Friend

friend_tz = unique(friend1[,5])

for (i in 1:nrow(friend1)) {
  p1 = friend1[i, 2]
  q1 = friend1[i, 3]
  weight = friend1[i, 4]

  if (is.null(mapid[[p1]]) || is.null(mapid[[q1]])) {

  } else {

    pp1 = mapid[[p1]]
    qq1 = mapid[[q1]]
    if (pp1 != qq1) {
      tf = friend1[i, 5]
      kk = which(friend_tz == tf)
      if (zij[pp1,qq1,3,kk] == 0 && weight > 0) {
        zij[pp1,qq1,3,kk] = 1
      }
      if (zij[qq1,pp1,3,kk] == 0 && weight > 0) {
        zij[qq1,pp1,3,kk] = 1
      }
      # zij[pp1,qq1,3,kk] = zij[pp1,qq1,3,kk] + weight
      # zij[qq1,pp1,3,kk] = zij[qq1,pp1,3,kk] + weight
    }
  }
}

trail[,1] = as.numeric(trail[,1])
trail[,2] = as.numeric(trail[,2])

table(trail[,2])

t_sep_t = seq(1, 100, 1)/100

# Run Real Data -----------------------------------------------------------

p = 3
npara = nonParametric(trail, zij, n, 3, tz = tz)
npara_homo = nonParametric_homo(trail, zij, n, 3, tz = tz)
npara_del_1 = nonParametric(trail, zij[,,-1,], n, 2, tz = tz)
npara_del_2 = nonParametric(trail, zij[,,-2,], n, 2, tz = tz)
npara_del_3 = nonParametric(trail, zij[,,-3,], n, 2, tz = tz)



# Partial R^2 -------------------------------------------------------------


k = seq(5, 95, 5)

p = 3
RSS = rep(0, length(k))
RSS_homo = rep(0, length(k))
RSS_del_1 = rep(0, length(k))
RSS_del_2 = rep(0, length(k))
RSS_del_3 = rep(0, length(k))

for (z in k) {

  ts = z
  te = z + 5

  Nij = matrix(0, nrow = n, ncol = n)

  for (i in 1:nrow(trail)) {
    if (trail[i,3] > t_sep_t[ts] && trail[i,3] <  t_sep_t[te] ) {
      p1 = trail[i,1]
      q1 = trail[i,2]
      Nij[p1, q1] = Nij[p1, q1] + 1
    }
  }

  Lij = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        k2 = which.max(z/100 < tz)
        Lij[i, j] = Lij[i, j] + sum(npara$nonhomo_coefficients[, te] * zij[i,j,,k2])
        Lij[i, j] = Lij[i, j] - sum(npara$nonhomo_coefficients[, ts] * zij[i,j,,k2])
        Lij[i, j] = Lij[i, j] + npara$homo_coefficients$outgoing[i, te] + npara$homo_coefficients$incoming[j, te]
        Lij[i, j] = Lij[i, j] - npara$homo_coefficients$outgoing[i, ts] - npara$homo_coefficients$incoming[j, ts]
        Lij[i, j] = Lij[i, j] + npara$homo_coefficients$baseline[te] - npara$homo_coefficients$baseline[ts]
      }
    }
  }

  Lij0 = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        k2 = which.max(z/100 < tz)
        Lij0[i, j] = Lij0[i, j] + sum(npara_homo$nonhomo_coefficients[, te] * zij[i,j,,k2])
        Lij0[i, j] = Lij0[i, j] - sum(npara_homo$nonhomo_coefficients[, ts] * zij[i,j,,k2])
        Lij0[i, j] = Lij0[i, j] + npara_homo$homo_coefficients$baseline[te] - npara_homo$homo_coefficients$baseline[ts]
      }
    }
  }

  Lij1 = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        k2 = which.max(z/100 < tz)
        Lij1[i, j] = Lij1[i, j] + sum(npara_del_1$nonhomo_coefficients[, te] * zij[i,j,-1,k2])
        Lij1[i, j] = Lij1[i, j] - sum(npara_del_1$nonhomo_coefficients[, ts] * zij[i,j,-1,k2])
        Lij1[i, j] = Lij1[i, j] + npara_del_1$homo_coefficients$outgoing[i, te] + npara_del_1$homo_coefficients$incoming[j, te]
        Lij1[i, j] = Lij1[i, j] - npara_del_1$homo_coefficients$outgoing[i, ts] - npara_del_1$homo_coefficients$incoming[j, ts]
        Lij1[i, j] = Lij1[i, j] + npara_del_1$homo_coefficients$baseline[te] - npara_del_1$homo_coefficients$baseline[ts]
      }
    }
  }

  Lij2 = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        k2 = which.max(z/100 < tz)
        Lij2[i, j] = Lij2[i, j] + sum(npara_del_2$nonhomo_coefficients[, te] * zij[i,j,-2,k2])
        Lij2[i, j] = Lij2[i, j] - sum(npara_del_2$nonhomo_coefficients[, ts] * zij[i,j,-2,k2])
        Lij2[i, j] = Lij2[i, j] + npara_del_2$homo_coefficients$outgoing[i, te] + npara_del_2$homo_coefficients$incoming[j, te]
        Lij2[i, j] = Lij2[i, j] - npara_del_2$homo_coefficients$outgoing[i, ts] - npara_del_2$homo_coefficients$incoming[j, ts]
        Lij2[i, j] = Lij2[i, j] + npara_del_2$homo_coefficients$baseline[te] - npara_del_2$homo_coefficients$baseline[ts]
      }
    }
  }

  Lij3 = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        k2 = which.max(z/100 < tz)
        Lij3[i, j] = Lij3[i, j] + sum(npara_del_3$nonhomo_coefficients[, te] * zij[i,j,-3,k2])
        Lij3[i, j] = Lij3[i, j] - sum(npara_del_3$nonhomo_coefficients[, ts] * zij[i,j,-3,k2])
        Lij3[i, j] = Lij3[i, j] + npara_del_3$homo_coefficients$outgoing[i, te] + npara_del_3$homo_coefficients$incoming[j, te]
        Lij3[i, j] = Lij3[i, j] - npara_del_3$homo_coefficients$outgoing[i, ts] - npara_del_3$homo_coefficients$incoming[j, ts]
        Lij3[i, j] = Lij3[i, j] + npara_del_3$homo_coefficients$baseline[te] - npara_del_3$homo_coefficients$baseline[ts]
      }
    }
  }

  RSS[z/5] = sum((Nij - Lij)^2)
  RSS_homo[z/5] = sum((Nij - Lij0)^2)
  RSS_del_1[z/5] = sum((Nij - Lij1)^2)
  RSS_del_2[z/5] = sum((Nij - Lij2)^2)
  RSS_del_3[z/5] = sum((Nij - Lij3)^2)


  cat(z, "\n")
}

1-RSS/RSS_homo
1-RSS/RSS_del_1
1-RSS/RSS_del_2
1-RSS/RSS_del_3

plot_t1 = seq(as.Date("2010/03/18"),as.Date("2011/07/15"),length = 19)

df1 = data.frame(x = rep(plot_t1), y = c(1-RSS/RSS_homo), group = c(rep("heterogeneous", 19)))

q1 = ggplot(df1, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#E69F00") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.4)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab(TeX('Partial $R^2')) +
  ggtitle("heterogeneous")

df2 = data.frame(x = rep(plot_t1), y = c(1-RSS/RSS_del_1), group = c(rep("couple", 19)))

q2 = ggplot(df2, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.4)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab(TeX('Partial $R^2')) + ggtitle("couple without children")


df3 = data.frame(x = rep(plot_t1), y = c(1-RSS/RSS_del_2))

q3 = ggplot(df3, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#009E73") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.4)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab(TeX('Partial $R^2')) +
  ggtitle("couple with children")

df4 = data.frame(x = rep(plot_t1), y = c(1-RSS/RSS_del_3))

q4 = ggplot(df4, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#F0E442") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.4)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab(TeX('Partial $R^2')) +
  ggtitle("friend")

ggarrange(q1, q4, q2, q3, ncol = 2, nrow = 2)

# Correlation plot --------------------------------------------------------


k = seq(5, 95, 5)

p = 3
corr1 = rep(0, length(k))
corrZ = matrix(NA, nrow = p, ncol = length(k))

for (z in k) {

  ts = z
  te = z + 5

  Nij = matrix(0, nrow = n, ncol = n)

  for (i in 1:nrow(trail)) {
    if (trail[i,3] > t_sep_t[ts] && trail[i,3] <  t_sep_t[te] ) {
      p1 = trail[i,1]
      q1 = trail[i,2]
      Nij[p1, q1] = Nij[p1, q1] + 1
    }
  }

  for (z1 in 1:p) {

    Lij = matrix(0, nrow = n, ncol = n)

    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          k2 = which.max(z/100 < tz)
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
corrCP = corrZ[1,]
corrCh = corrZ[2,]
corrFr = corrZ[3,]

plot_t1 = seq(as.Date("2010/03/18"),as.Date("2011/07/15"),length = 19)

df1 = data.frame(x = rep(plot_t1), y = c(corrAB), group = c(rep("heterogeneous", 19)))

q1 = ggplot(df1, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#E69F00") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.7)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab("Correlation") +
  ggtitle("heterogeneous")

df2 = data.frame(x = rep(plot_t1), y = c(corrCP), group = c(rep("couple", 19)))

q2 = ggplot(df2, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.7)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab("Correlation") + ggtitle("couple without children")


df3 = data.frame(x = rep(plot_t1), y = c(corrCh))

q3 = ggplot(df3, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#009E73") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.05, 0.7)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab("Correlation") +
  ggtitle("couple with children")

df4 = data.frame(x = rep(plot_t1), y = c(corrFr))

q4 = ggplot(df4, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "#F0E442") +
  scale_fill_discrete(name = "") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  scale_y_continuous(limits = c(-0.3, 0.7)) +
  theme(legend.position="top", panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  xlab("Time") +
  ylab("Correlation") +
  ggtitle("friend")

ggarrange(q1, q4, q2, q3, ncol = 2, nrow = 2)


# motivation --------------------------------------------------------------


library(scales)
library(igraph)
library(plot.matrix)

# Plot matrix

zij1 = array(0, dim = c(n, n, 1, 4))

for (i in 1:nrow(friend1)) {
  p1 = friend1[i, 2]
  q1 = friend1[i, 3]
  weight = friend1[i, 4]

  if (is.null(mapid[[p1]]) || is.null(mapid[[q1]])) {

  } else {

    pp1 = mapid[[p1]]
    qq1 = mapid[[q1]]
    if (pp1 != qq1) {
      tf = friend1[i, 5]
      kk = which(friend_tz == tf)
      # if (zij[pp1,qq1,3,kk] == 0 && weight > 0) {
      #   zij[pp1,qq1,3,kk] = 1
      # }
      # if (zij[qq1,pp1,3,kk] == 0 && weight > 0) {
      #   zij[qq1,pp1,3,kk] = 1
      # }
      zij1[pp1,qq1,1,kk] = zij1[pp1,qq1,1,kk] + weight
      zij1[qq1,pp1,1,kk] = zij1[qq1,pp1,1,kk] + weight
    }
  }
}

zij1[,,,] = zij1[,,,]/2

par(mar=c(5.1, 4.1, 4.1, -0.2))
plot(zij1[,,,1], border = NA, xlab = "Individual", ylab = "Individual", main = "2010-09-01")
plot(zij1[,,,4], border = NA, xlab = "Individual", ylab = "Individual", main = "2011-05-01")

plot(zij[,,3,1]==1, col = c("white", "#CCCCCC"), border = NA, xlab = "Individual", ylab = "Individual", main = "2010-09-01")
plot(zij[,,3,4]==1, col = c("white", "#CCCCCC"), border = NA, xlab = "Individual", ylab = "Individual", main = "2011-05-01")


# Plot in-out degree

t_sep_t = seq(0, 1, 1/16)
in_ave = matrix(0, nrow = n, ncol = length(t_sep_t))
out_ave = matrix(0, nrow = n, ncol = length(t_sep_t))

k = 1
for (i in 1:nrow(trail)) {
  t0 = t_sep_t[k]
  p1 = trail[i, 1]
  p2 = trail[i, 2]
  tn = trail[i, 3]
  if (tn > t0) {
    k = k + 1
  }

  if (k > length(t_sep_t)) {
    break
  }

  in_ave[p2, k] = in_ave[p2, k]+1
  out_ave[p1, k] = out_ave[p1, k]+1

}

in_ave_dev = in_ave[, -1]
out_ave_dev = out_ave[, -1]

for (k in 2:length(t_sep_t)) {
  in_ave[, k] = in_ave[, k] + in_ave[, k-1]
  out_ave[, k] = out_ave[, k] + out_ave[, k-1]
}

(as.Date("2011/07/15") - as.Date("2010/03/18")) * 0.0625

plot_t1 = seq(as.Date("2010/03/18"),as.Date("2011/07/15"),length = 17)

p4 = data.frame(t = rep(plot_t1[-1], 4),
                y = c(t(in_ave_dev[c(1,2,3,4),])),
                group = c(rep("sp10-01-07", length(t_sep_t) - 1),
                          rep("sp10-01-08", length(t_sep_t) - 1),
                          rep("sp10-01-48", length(t_sep_t) - 1),
                          rep("sp10-01-38", length(t_sep_t) - 1)))

ggplot(p4, aes(x = t, y = y)) +
  geom_line(aes(color = group), size = 0.75) +
  scale_y_continuous(breaks=seq(0, 500, 100), limits = c(0, 550)) +
  scale_color_discrete(name="")+
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  ylab("Monthly In-degree") +
  theme(panel.background = element_rect(fill = "white"), legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

p5 = data.frame(t = rep(plot_t1[-1], 4),
                y = c(t(out_ave_dev[c(1,2,3,4),])),
                group = c(rep("sp10-01-07", length(t_sep_t) - 1),
                          rep("sp10-01-08", length(t_sep_t) - 1),
                          rep("sp10-01-48", length(t_sep_t) - 1),
                          rep("sp10-01-38", length(t_sep_t) - 1)))


ggplot(p5, aes(x = t, y = y)) +
  geom_line(aes(color = group), size = 0.75) +
  scale_y_continuous(breaks=seq(0, 500, 100), limits = c(0, 550)) +
  scale_color_discrete(name="")+
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  ylab("Monthly Out-degree") +
  theme(panel.background = element_rect(fill = "white"), legend.position = c(0.15, 0.8)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

in_ave_dev_person = apply(in_ave_dev, 2, mean)
in_ave_dev_sd = apply(in_ave_dev, 2, sd)

df_degree = data.frame(x = plot_t1[-1], y = c(in_ave_dev_person),
                yl = c(pmax(in_ave_dev_person - 1.96 * in_ave_dev_sd, 0)),
                yu = c(in_ave_dev_person + 1.96 * in_ave_dev_sd))

ggplot(df_degree, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  ylab("Monthly average In-degree")  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")

# Homo

homo_effect = matrix(0, nrow = 2, ncol = length(t_sep_t))

k = 1
for (i in 1:nrow(trail)) {
  t0 = t_sep_t[k]
  p1 = trail[i, 1]
  p2 = trail[i, 2]
  tn = trail[i, 3]
  if (tn > t0) {
    k = k + 1
  }

  if (k > length(t_sep_t)) {
    break
  }

  if (zij[p1, p2, 1, 1] == 1 || zij[p1, p2, 2, 1]) {
    homo_effect[1, k] = homo_effect[1, k] + 1
  } else {
    homo_effect[2, k] = homo_effect[2, k] + 1
  }
}

homo_effect = homo_effect[,-1]
sf = sum(zij[,,1,1]+zij[,,2,1])
nsf = n*(n-1) - sf

p3 = data.frame(t = rep(plot_t1[-1], 2),
                y = c(homo_effect[1,], homo_effect[2,]),
                group = c(rep("copule", length(t_sep_t) - 1),
                          rep("non couple", length(t_sep_t) - 1)))


ggplot(p3, aes(x = t, y = y)) +
  geom_line(aes(color = group), size = 0.75) +
  # scale_y_continuous(breaks=seq(0, 15, 5), limits = c(0, 15)) +
  scale_color_discrete(name="", labels = c("couple", "non couple"))+
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  ylab("Monthly number of SMS") +
  theme(panel.background = element_rect(fill = "white"), legend.position = c(0.2, 0.85)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# Out-degree sp10-01-07

library(tidyr)
library(dplyr)

o1 = which(trail[,1] == 1)
dfo1 = data.frame(t = as.Date(sms1[o1, 3]), id = trail[o1, 2])
dfevent = matrix(0, nrow = nrow(dfo1), ncol = n-1)
for (i in 2:n) {
  dfevent[,i-1] = dfo1$id == i
}
dfevent = as.data.frame(dfevent)
colnames(dfevent) = sprintf("%02d", 2:28)
dfevent$t = dfo1$t

data <- dfevent %>% pivot_longer(-c(t), names_to = 'id')

ggplot(data) + geom_point(aes(x = t, y = 2))+
  geom_vline(aes(xintercept  = t), data = filter(data,value == 1)) +
  xlab(expression(italic("t"))) +
  scale_x_date(limits = c(as.Date("2010-03-18"), as.Date("2011/07/15"))) +
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~id, ncol = 1, strip.position = 'left') +
  theme(axis.title.y       = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())



# In-degree sp10-01-07

o2 = which(trail[,2] == 1)
dfo2 = data.frame(t = as.Date(sms1[o2, 3]), id = trail[o2, 1])
dfevent2 = matrix(0, nrow = nrow(dfo2), ncol = n-1)
for (i in 2:n) {
  dfevent2[,i-1] = dfo2$id == i
}
dfevent2 = as.data.frame(dfevent2)
colnames(dfevent2) = sprintf("%02d", 2:28)
dfevent2$t = dfo2$t

data2 <- dfevent2 %>% pivot_longer(-c(t), names_to = 'id')

ggplot(data2) + geom_point(aes(x = t, y = 2))+
  geom_vline(aes(xintercept  = t), data = filter(data2,value == 1)) +
  scale_x_date(limits = c(as.Date("2010-03-18"), as.Date("2011/07/15"))) +
  xlab(expression(italic("t"))) +
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~id, ncol = 1, strip.position = 'left') +
  theme(axis.title.y       = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())




# Baseline plot -----------------------------------------------------------


plot_t1 = seq(as.Date("2010/03/18"),as.Date("2011/07/15"),length = 100)

# Baseline

Nij = matrix(0, nrow = n, ncol = n)
sdB = rep(0, 100)
sdA = rep(0, 100)
co = 1
for (i in 1:nrow(trail)) {

  p = trail[i,1]
  q = trail[i,2]
  t = trail[i,3]
  Nij[p, q] = Nij[p, q] + 1

  if (t >= t_sep_t[co]) {
    sdB[co] =  sqrt(sum(Nij) / n^4)
    sdA[co] =  sum(Nij)
    co = co + 1
  }
}

bsestimate = npara$homo_coefficients$baseline

df = data.frame(x = plot_t1, y = c(bsestimate),
                yl = c(bsestimate - 1.96 * sdB),
                yu = c(bsestimate + 1.96 * sdB))

ggplot(df, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.2, linetype = "dashed") +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  ylab("Cumulative baseline hazard")  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")


# Plot alp bet ------------------------------------------------------------


# alp bet id 1 07

k = 1
roledf = data.frame(x = plot_t1,
                    y = npara$ab$outgoing[k, ],
                    yl = npara$ab$outgoing[k, ] - 1.96 * npara$ab$sdout[k, ],
                    yu = npara$ab$outgoing[k, ] + 1.96 * npara$ab$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-20, 130, 20), limits = c(-20, 130)) +
  ylab(TeX('$\\widehat{\\alpha_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 1
roledf = data.frame(x = plot_t1,
                    y = npara$ab$incoming[k, ],
                    yl = npara$ab$incoming[k, ] - 1.96 * npara$ab$sdinc[k, ],
                    yu = npara$ab$incoming[k, ] + 1.96 * npara$ab$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-20, 130, 20), limits = c(-20, 130)) +
  ylab(TeX('$\\widehat{\\beta_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# alp bet id 2 08

k = 2
roledf = data.frame(x = plot_t1,
                    y = npara$ab$outgoing[k, ],
                    yl = npara$ab$outgoing[k, ] - 1.96 * npara$ab$sdout[k, ],
                    yu = npara$ab$outgoing[k, ] + 1.96 * npara$ab$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-10, 170, 20), limits = c(-10, 170)) +
  ylab(TeX('$\\widehat{\\alpha_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 2
roledf = data.frame(x = plot_t1,
                    y = npara$ab$incoming[k, ],
                    yl = npara$ab$incoming[k, ] - 1.96 * npara$ab$sdinc[k, ],
                    yu = npara$ab$incoming[k, ] + 1.96 * npara$ab$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-10, 170, 20), limits = c(-10, 170)) +
  ylab(TeX('$\\widehat{\\beta_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# alp bet id 11 46

k = 11
roledf = data.frame(x = plot_t1,
                    y = npara$ab$outgoing[k, ],
                    yl = npara$ab$outgoing[k, ] - 1.96 * npara$ab$sdout[k, ],
                    yu = npara$ab$outgoing[k, ] + 1.96 * npara$ab$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-10, 90, 20), limits = c(-10, 90)) +
  ylab(TeX('$\\widehat{\\alpha_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 11
roledf = data.frame(x = plot_t1,
                    y = npara$ab$incoming[k, ],
                    yl = npara$ab$incoming[k, ] - 1.96 * npara$ab$sdinc[k, ],
                    yu = npara$ab$incoming[k, ] + 1.96 * npara$ab$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-10, 90, 20), limits = c(-10, 90)) +
  ylab(TeX('$\\widehat{\\beta_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# Plot A B ----------------------------------------------------------------


# A B id 1 07

k = 1
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$outgoing[k, ],
                    yl = npara$homo_coefficients$outgoing[k, ] - 1.96 * npara$homo_coefficients$sdout[k, ],
                    yu = npara$homo_coefficients$outgoing[k, ] + 1.96 * npara$homo_coefficients$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 22, 4), limits = c(-2, 22)) +
  ylab(TeX('$\\widehat{A_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 1
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$incoming[k, ],
                    yl = npara$homo_coefficients$incoming[k, ] - 1.96 * npara$homo_coefficients$sdinc[k, ],
                    yu = npara$homo_coefficients$incoming[k, ] + 1.96 * npara$homo_coefficients$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 22, 4), limits = c(-2, 22)) +
  ylab(TeX('$\\widehat{B_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# A B id 2 08

k = 2
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$outgoing[k, ],
                    yl = npara$homo_coefficients$outgoing[k, ] - 1.96 * npara$homo_coefficients$sdout[k, ],
                    yu = npara$homo_coefficients$outgoing[k, ] + 1.96 * npara$homo_coefficients$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 50, 13), limits = c(-2, 50)) +
  ylab(TeX('$\\widehat{A_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 2
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$incoming[k, ],
                    yl = npara$homo_coefficients$incoming[k, ] - 1.96 * npara$homo_coefficients$sdinc[k, ],
                    yu = npara$homo_coefficients$incoming[k, ] + 1.96 * npara$homo_coefficients$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 50, 13), limits = c(-2, 50)) +
  ylab(TeX('$\\widehat{B_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# A B id 11 46

k = 11
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$outgoing[k, ],
                    yl = npara$homo_coefficients$outgoing[k, ] - 1.96 * npara$homo_coefficients$sdout[k, ],
                    yu = npara$homo_coefficients$outgoing[k, ] + 1.96 * npara$homo_coefficients$sdout[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 22, 8), limits = c(-2, 22)) +
  ylab(TeX('$\\widehat{A_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


k = 11
roledf = data.frame(x = plot_t1,
                    y = npara$homo_coefficients$incoming[k, ],
                    yl = npara$homo_coefficients$incoming[k, ] - 1.96 * npara$homo_coefficients$sdinc[k, ],
                    yu = npara$homo_coefficients$incoming[k, ] + 1.96 * npara$homo_coefficients$sdinc[k, ])

ggplot(roledf, aes(x = x, y = y)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.3) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(breaks=seq(-2, 22, 8), limits = c(-2, 22)) +
  ylab(TeX('$\\widehat{B_i}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# Plot theta --------------------------------------------------------------

plot_t1 = seq(as.Date("2010/03/18"),as.Date("2011/07/15"),length = 100)


# Thetas

p = dim(zij)[3]
tk = seq(1, 100, 1)
PaT = npara$Pa[(2*n):(2*n+p-1),]

rd_theta_sd = matrix(0, nrow = p, ncol = 100)
rd_theta_ab_sd = matrix(0, nrow = p, ncol = 100)

for (j in 1:100) {
  NTs = npara$NTs[, tk[j]]
  rd_theta_ab_sd[,j] = diag(PaT %*% diag(NTs) %*% t(PaT))


  NT = npara$NT[, tk[j]]
  rd_theta_sd[,j] = diag(PaT %*% diag(NT) %*% t(PaT))

  cat(j, "\n")
}

groupName = c("copule without children", "copule with children", "friend")

# tht

tht = data.frame()
for (i in seq_len(p)) {
  thtdf = data.frame(x = plot_t1,
                     y = npara$th[i,],
                     yl = npara$th[i,] - 1.96 * sqrt(rd_theta_ab_sd[i,]),
                     yu = npara$th[i,] + 1.96 * sqrt(rd_theta_ab_sd[i,]),
                     group = rep(groupName[i], 100))
  tht = rbind(tht, thtdf)

}

ggplot(tht, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\theta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), legend.position="top")

# T

tht = data.frame()
for (i in seq_len(p)) {
  thtdf = data.frame(x = plot_t1,
                     y = npara$nonhomo_coefficients[i,],
                     yl = npara$nonhomo_coefficients[i,] - 1.96 * sqrt(rd_theta_sd[i,]),
                     yu = npara$nonhomo_coefficients[i,] + 1.96 * sqrt(rd_theta_sd[i,]),
                     group = rep(groupName[i], 100))
  tht = rbind(tht, thtdf)

}

ggplot(tht, aes(x = x, y = y, fill = group)) + geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0.5) +
  scale_x_date(date_labels = "%d-%b-%y", breaks = "12 week") +
  xlab(expression(italic("t"))) +
  scale_fill_discrete(name = "") +
  ylab(TeX('$\\widehat{\\Theta}$(t)'))  +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), legend.position="top")

# Tests -------------------------------------------------------------------

directTest <- function(inc, incvar, out, outvar, tsep) {
  n = nrow(inc)
  s0 = -100
  for (i in tsep) {
    s = sum((inc[, i] - out[, i])^2/(incvar[, i] + outvar[, i]) - 1) / sqrt(n) / sqrt(2)
    s0 = max(s0, s)
  }
  return(list(tstat = s0,
              ind = s0 > qnorm(0.95^(1/length(tsep)))))
}

directTest(npara$ab$outgoing, npara$ab$sdout, npara$ab$incoming, npara$ab$sdinc, seq(10,90,10))

