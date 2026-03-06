

library(distrEx)

x <- rnorm(100)

TotalVarDist(Norm(), x)


M1  <- UnivarMixingDistribution(Norm(3,.25), Norm(1,.5))

plot(Norm())

plot(M1)

TotalVarDist(e1 = M1, e2 = x)
