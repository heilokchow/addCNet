np0 = nonParametric(result$trail, array(zij, c(n,n,p,1)), n, p, h1 = 0.05, test = 0)

Pa = np0$Pa
NTs = diag(np0$NTs[,50])
test = Pa %*% NTs %*% t(Pa)

NT = diag(np0$NT[,50])
test = Pa %*% NT %*% t(Pa)


diag(test)[2:n]
np0$homo_coefficients$sdout[-1,50]
np0$ab$sdout[-1,50]

plot(np0$homo_coefficients$sdout[-1,50], diag(test)[2:n])

test = Pa %*% np0$NTs[,50]
np0$ab$outgoing[,50]
