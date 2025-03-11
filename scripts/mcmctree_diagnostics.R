# MCMCTree diagnostics
library(ape)
setwd('~/workspace/Pterocarya/species_phylogeny/MCMCTree_calibration/')

# plot tree and get relevant node labels 
tre <- read.tree("NodeLabels.tree")
plot(tre)
nodelabels(text=tre$node.label)

# Root - 21
# Pterocarya crown age - 27
# Juglans crown age - 31
# Carya crown age - 37
# Pterocarya - Cyclocarya divergence - 26


# read in MCMC trace files
mcmc1 <- read.table("mcmc_rep1/mcmc1.txt", head=TRUE)
mcmc2 <- read.table("mcmc_rep2/mcmc2.txt", head=TRUE)

# each data frame contains 23 columns:
# MCMC generation number, node ages (divergence times), mean mutation rate,
# rate drift coefficients, and sample log-likelihood values

names(mcmc1)

# to check for convergence of the MCMC runs, we calculate the posterior
# means of times for each run, and plot them against each other
t.mean1 <- apply(mcmc1[,2:20], 2, mean) * 100
t.mean2 <- apply(mcmc2[,2:20], 2, mean) * 100
# good convergence is indicated when the points fall on the y = x line.
# posterior times for run 1 vs run 2:
plot(t.mean1, t.mean2, main="a) Posterior times, chain 1 vs. chain 2"); abline(0, 1)



# ###############################################
# PRIOR VS POSTERIOR
# ###############################################
# prior is determined by fossil calibrations (user-input) 
# and the birth-death model in mcmctree. 
# prior distributions may be distorted slightly from original
# specifications due to the constraints of the birth-death model
mcmc1.p <- read.table("prior1/mcmc_prior1.txt", head=TRUE)

dev.off()
par(mfrow=c(4,2), mar = c(4,4,2,1), mgp = c(2,0.5,0))

# ------ Pterocarya crown age -----

dpr <- density(mcmc1.p[,'t_n27'], adj=.1) # prior of specific node
dPr1 <- density(mcmc1[,'t_n27'], adj=.1)   # Posterior of the same node
dPr2 <- density(mcmc2[,'t_n27'], adj=.1)   # Posterior of the same node

# plotting range
xl <- range(c(dpr$x, dPr2$x))
yl <- range(c(dpr$y, dPr2$y))

# check that it resembles the specified gamma distribution 
plot(dpr, main=expression(bolditalic("Pterocarya")~bold("crown age")), 
    xlim = xl, ylim = yl, col = 'gray', xlab = '100 Myr')

#axis(1, mgp = c(0,1,0), cex = 0.7, line = -1)
#axis(2, las = 2, mgp = c(0,1,0), cex = 0.7)
#mtext("100 Myr", side = 1, line = 2, cex = 0.7, tck = 0.5)
#mtext("Density", side = 2, line = 2, cex = 0.7)

lines(dPr1, col="black")
lines(dPr2, col="blue")

# trace plot
plot(mcmc1$t_n27, ty='l', main="", xlab = 'Generation', ylab = '100 Myr')
#mtext("Generation", side = 1, line = 2, cex = 0.7, tck = 0.5)


# ---- Walnut crown age prior vs posterior----

dpr <- density(mcmc1.p[,'t_n31'], adj=.1) # prior of specific node
dPr1 <- density(mcmc1[,'t_n31'], adj=.1)   # Posterior of the same node
dPr2 <- density(mcmc2[,'t_n31'], adj=.1)   # Posterior of the same node

# plotting range
xl <- range(c(dpr$x, dPr2$x))
yl <- range(c(dpr$y, dPr2$y))

# prior resembles the specified gamma distribution 
plot(dpr, main=expression(bolditalic("Juglans")~bold("crown age")), xlab="100 myr", ylab="Density", las=1, xlim=xl, ylim=yl, col="darkgrey")
lines(dPr1, col="black")
lines(dPr2, col="blue")

# trace plot
plot(mcmc1$t_n31, ty='l', main="", xlab = 'Generation', ylab = '100 Myr')


# ---- Carya crown age prior vs posterior----

dpr <- density(mcmc1.p[,'t_n37'], adj=.1) # prior of specific node
dPr1 <- density(mcmc1[,'t_n37'], adj=.1)   # Posterior of the same node
dPr2 <- density(mcmc2[,'t_n37'], adj=.1)   # Posterior of the same node

# plotting range
xl <- range(c(dpr$x, dPr2$x))
yl <- range(c(dpr$y, dPr2$y))

# prior resembles the specified gamma distribution 
plot(dpr, main=expression(bolditalic("Carya")~bold("crown age")), xlab="100 Myr", ylab="Density", las=1, xlim=xl, ylim=yl, col="darkgrey")
lines(dPr1, col="black")
lines(dPr2, col="blue")

# trace plot
plot(mcmc1$t_n37, ty='l', main="", xlab = 'Generation', ylab = '100 Myr')




# ---- Cyclocarya, Pterocarya divergence -----

dpr <- density(mcmc1.p[,'t_n26'], adj=.1) # prior of specific node
dPr1 <- density(mcmc1[,'t_n26'], adj=.1)   # Posterior of the same node
dPr2 <- density(mcmc2[,'t_n26'], adj=.1)   # Posterior of the same node

# plotting range
xl <- range(c(dpr$x, dPr2$x))
yl <- range(c(dpr$y, dPr2$y))

# check that it resembles the specified gamma distribution 
plot(dpr, main=expression(bolditalic("Pterocarya Cyclocarya")~bold("divergence")), xlab="100 Myr", ylab="Density", las=1, xlim=xl, ylim=yl, col="darkgrey")
lines(dPr1, col="black")
lines(dPr2, col="blue")


# trace plot
plot(mcmc1$t_n26, ty='l', main="", xlab = 'Generation', ylab = '100 Myr')

mcmc2



# # ---- root prior ----
# 
# dpr <- density(mcmc1.p[,'t_n21'], adj=.1) # prior of specific node
# dPr1 <- density(mcmc1[,'t_n21'], adj=.1)   # Posterior of the same node
# dPr2 <- density(mcmc2[,'t_n21'], adj=.1)   # Posterior of the same node
# 
# # plotting range
# xl <- range(c(dpr$x, dPr2$x))
# yl <- range(c(dpr$y, dPr2$y))
# 
# 
# # check that it resembles the specified gamma distribution 
# plot(dpr, main="t_n21", xlab="", ylab="", las=1, xlim=xl, ylim=yl, col="darkgrey")
# lines(dPr1, col="black")
# lines(dPr2, col="blue")

