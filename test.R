V_0 <- matrix(0, nrow = 6, ncol = 6)
diag(V_0) <- 1
u_0 <- 7

setwd("~/Dropbox/Projects/RA Biomarker/")

m <- jags.model("code/jags/test.jags", data = list(V_0 = V_0, u_0 = u_0))

jags.samples(m, 100000, variable.names = "Omega_0")
