# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------
library(copula)
library(scatterplot3d)
library(ggplot2)
library(grid)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(plot3D)
library(rglwidget)
library(knitr)
library(rgl)
library(VineCopula)
library(VC2copula)
library(asbio)
library(lcopula)
library(goftest)

# ------------------------------------------------------------------------------
# DATA IMPORT
# ------------------------------------------------------------------------------

# setwd(dir = "C:/Users/gwend/Documents/ISFA/3A/SEMESTRE 2/ESTIMATION DE COPULES/PROJET") # Florian
# setwd("~/Desktop/Estimation de copules") # Tina

# dat <- read.csv(file = "chinesetyph.csv") # Florian
# dat <- read.csv(file = "ChineseTyphoon.csv") # Tina

# ------------------------------------------------------------------------------
# VARIABLES DESCRIPTION
# ------------------------------------------------------------------------------

# YYYYNNMMDDHH YYYY-year NN-code MM-month DD-day HH-hour
# Latitude
# Longitude
# Lowest pressure
# Max wind speed
# Radius in 34n mile/h (tropical storms +)
# Satellite Sensor

head(dat)

# ------------------------------------------------------------------------------
# DATA STUDY AND FIRST DEFINITIONS
# ------------------------------------------------------------------------------
dat <- dat[, c(-1, -7)]
n <- nrow(dat)

pressu <- dat[, 3]
windspeed <- dat[, 4]
sir <- dat[, 5]

par(mar = c(4, 4, .1, .1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
df_ranked = as.data.frame(lapply(dat, rank))/n

pairs(dat, pch = '.', col = "deepskyblue4", upper.panel = NULL, main = "Diagrammes de dispersion")
pairs(df_ranked, pch = '.', col = "deepskyblue4", upper.panel = NULL, main = "Rank-Rank Plot")

# Latitude and longitude do not greatly affect other variables. Thus, we will stick to the relations between the PRS, WND and SiR34 variables

dat <- dat[, c(-1, -2)]
df_ranked = as.data.frame(lapply(dat, rank))/n

pairs(dat, pch = '.', col = "deepskyblue4", upper.panel = NULL, main = "Diagrammes de dispersion")
pairs(df_ranked, pch = '.', col = "deepskyblue4", upper.panel = NULL, main = "Rank-Rank Plot")

# Representation of the correlation through a 3D histogram
cols <- colorRampPalette(c("white", "deepskyblue4"))(100)

par(mfrow=c(1,1))

# PRS & WND
x_c <- cut(rank(dat$PRS)/length(dat), 30)
y_c <- cut(rank(dat$WND)/length(dat), 30)
z <- table(x_c, y_c)/sum(table(x_c, y_c))
m = max(z)
hist3D(z = z, border = "black", xlab = "PRS", ylab = "WND", zlab = "", theta = 130, phi = 25, 
       col = cols, #zlim=c(0,m),
       colkey = list(labels = F, axis = F), nticks = 2, ticktype = "detailed")
filled.contour(z)

# PRS & SiR34
x_c <- cut(rank(dat$PRS)/length(dat), 30)
y_c <- cut(rank(dat$SiR34)/length(dat), 30)
z <- table(x_c, y_c)/sum(table(x_c, y_c))
max(z)
hist3D(z = z, border = "black", xlab = "PRS", ylab = "SiR34", zlab = "", theta = 130, phi = 25, 
       col = cols, # zlim=c(0,m),
       colkey = list(labels = F, axis = F), nticks = 2, ticktype = "detailed")
filled.contour(z)

# WND & SiR34
x_c <- cut(rank(dat$WND)/length(dat), 30)
y_c <- cut(rank(dat$SiR34)/length(dat), 30)
z <- table(x_c, y_c)/sum(table(x_c, y_c))

hist3D(z = z, border = "black", xlab = "WND", ylab = "SiR34", zlab = "", theta = 220, phi = 25, 
       col = cols, # zlim=c(0,m),
       colkey = list(labels = F, axis = F), nticks = 2, ticktype = "detailed")

filled.contour(z)

# Chi-plot
par(mfrow = c(1, 1))

Chi_plot_PRS_WND <- chi.plot(windspeed, pressu, pch = '.', main = "Chi-plot de la relation\n pression/vitesse du vent", col ="skyblue3", 
                             xlim = c(-1, 1))
Chi_plot_PRS_Radius <- chi.plot(sir, pressu, pch = '.', main = "Chi-plot de la relation\n pression/rayon", col = "skyblue3")

# K-plot
par(mfrow = c(1, 1))

K.plot(data = dat[, c(1, 2)]) # Other arguments not working
K.plot(data = dat[, c(1, 3)])

# ------------------------------------------------------------------------------
# PARAMETRIC APPROACH
# ------------------------------------------------------------------------------

# Histograms of the marginals -----

### Pressure
hist(x = pressu, breaks = 25, freq = FALSE, col = "skyblue", main = "Histogramme de la pression", xlab = "Pression") # Gamma modelling looks appropriate
lines(density(x = pressu), col = 17, lwd = 4)

xmin <- min(pressu)
xmax <- max(pressu)

plot(dgamma(x = seq(xmin, xmax, 0.01), shape = 967, rate = 1), type = 'l') # Manual test

### Wind speed
hist(x = windspeed, breaks = 25, freq = FALSE, col = "salmon", main = "Histogramme de la vitesse du vent", xlab = "Vit. du vent") # Gamma modelling here too
lines(density(x = windspeed), col = 17, lwd = 4)

xmin <- min(windspeed)
xmax <- max(windspeed)

plot(dgamma(x = seq(xmin, xmax, 0.01), shape = 12, rate = 0.3), type = 'l')

### Radius
hist(x = sort(sir, decreasing = TRUE)[-1:-200], breaks = 25, freq = FALSE, col = "pink3", main = "Histogramme du rayon", xlab = "Rayon") # Normal distribution here
lines(density(x = sir), col = 17, lwd = 4)

xmin <- min(sir)
xmax <- max(sir)

plot(dnorm(x = seq(xmin, xmax, 0.01), mean = 200, sd = 50), type = 'l')

# Parameters fitting -----

# Empirical mean
mean_pressu <- mean(pressu)
mean_windspeed <- mean(windspeed)
mean_sir <- mean(sir)

# Empirical variance
var_pressu <- var(pressu)
var_windspeed <- var(windspeed)
var_sir <- var(sir)

# FITTING OF THE PRESSURE DISTRIBUTION PARAMETERS

# First attempt with the optim function
# Gamma parameters
params_pressu <- c(mean_pressu^2/var_pressu, var_pressu/mean_pressu)
loglikelihood.margins <- function(b, x) sum(log((x^(b[1]-1)*b[2]^(b[1])*exp(-b[2]*x))/gamma(b[1]))) # Gamma log likelihood
ctrl <- list(fnscale = -1) # Turns the problem into a maximization problem
paramhat_pressu <- optim(par = params_pressu, fn = loglikelihood.margins, x = pressu, control = ctrl)$par # Doesn't work due to a too high initial shape value

# Gamma fitting (second attempt)
paramhat_pressu_gamma <- fitdistr(x = pressu, densfun = "gamma", start = list(shape = params_pressu[1], scale = params_pressu[2]))

# Representing the simulated histogram
hist(rgamma(n = n, shape = paramhat_pressu_gamma$estimate[1], scale = paramhat_pressu_gamma$estimate[2]), freq = FALSE, col = "skyblue")
lines(density(pressu))

# params_pressu_lognorm <- c(exp(mean_sir + var_sir/2), (exp(var_sir)-1)*exp(2*mean_sir+var_sir))
paramhat_pressu_lognorm <- fitdistr(x = pressu, densfun = "lognormal")
hist(rlnorm(n = n, meanlog = paramhat_pressu_lognorm$estimate[1], sdlog = paramhat_pressu_lognorm$estimate[2]), freq = FALSE, col = "skyblue")
lines(density(pressu))

# Conclusion : Gamma and lognormal distributions finally don't work very well
# It's not working due to a negative skewness. As a result, it tends to reach a symmetrical distribution but can't go further

paramhat_pressu_weibull <- fitdistr(x = pressu, densfun = "weibull")

hist(rweibull(n = n*1000, shape = paramhat_pressu_weibull$estimate[1], scale = paramhat_pressu_weibull$estimate[2]), freq = FALSE, col = "skyblue", 
     xlab = "Pression", ylab = "Densité", main = "")   
lines(density(pressu), col = "red", lwd = 3)
  
# AIC/BIC criterion to choose the model
AIC(paramhat_pressu_gamma, paramhat_pressu_lognorm, paramhat_pressu_weibull) # We choose the Weibull distribution
BIC(paramhat_pressu_gamma, paramhat_pressu_lognorm, paramhat_pressu_weibull)

# FITTING OF THE WINDSPEED DISTRIBUTION PARAMETERS

# First attempt with the optim function
params_windspeed <- c(mean_windspeed^2/var_windspeed, var_windspeed/mean_windspeed)
paramhat_windspeed <- optim(par = params_windspeed, fn = loglikelihood.margins, x = windspeed, control = ctrl, )$par
paramhat_windspeed_gamma <- fitdistr(x = windspeed, densfun = "gamma", start = list(shape = params_windspeed[1], scale = params_windspeed[2]))

# Other attempts
paramhat_windspeed_lognorm <- fitdistr(x = windspeed, densfun = "lognormal")
paramhat_windspeed_weibull <- fitdistr(x = windspeed, densfun = "weibull")

# Representing the simulated histogram
hist(rgamma(n = n, shape = paramhat_windspeed[1], rate = paramhat_windspeed[2]), freq = FALSE, col = "salmon") # First method
hist(rgamma(n = n*1000, shape = paramhat_windspeed_gamma$estimate[1], scale = paramhat_windspeed_gamma$estimate[2]), freq = FALSE, col = "salmon") # Second method
lines(density(windspeed))

hist(rlnorm(n = n*1000, meanlog = paramhat_windspeed_lognorm$estimate[1], sdlog = paramhat_windspeed_lognorm$estimate[2]), freq = FALSE, col = "skyblue",
     xlab = "Vitesse du vent", ylab = "Densité", main = "")
lines(density(windspeed), col = "red", lwd = 3)

hist(rweibull(n = n*1000, shape = paramhat_windspeed_weibull$estimate[1], scale = paramhat_windspeed_weibull$estimate[2]), freq = FALSE, col = "salmon")
lines(density(windspeed))

# AIC/BIC criterion to choose the model
AIC(paramhat_windspeed_gamma, paramhat_windspeed_lognorm, paramhat_windspeed_weibull) # lognormal distribution
BIC(paramhat_windspeed_gamma, paramhat_windspeed_lognorm, paramhat_windspeed_weibull) 

# FITTING OF THE RADIUS DISTRIBUTION PARAMETERS
loglikelihood.radius <- function(b, x) sum((1/(sqrt(2*pi)*b[2]))*exp(-(((x-b[1])^2)/(2*(b[2]^2)))))
paramhat_sir <- optim(par = params_sir, fn = loglikelihood.radius, x = sir, control = ctrl)$par
paramhat_sir_normal <- fitdistr(x = sir, densfun = "normal")

paramhat_sir_gamma <- fitdistr(x = sir, densfun = "gamma")
paramhat_sir_lognorm <- fitdistr(x = sir, densfun = "log-normal")
paramhat_sir_weibull <- fitdistr(x = sir, densfun = "weibull")

# Representing the simulated histogram
hist(rnorm(n = n, mean = paramhat_sir_normal$estimate[1], sd = paramhat_sir_normal$estimate[2]), freq = FALSE, col = "pink3")
lines(density(sir))

hist(rlnorm(n = n, meanlog = paramhat_sir_lognorm$estimate[1], sdlog = paramhat_sir_lognorm$estimate[2]), freq = FALSE, col = "pink3")
lines(density(sir))

hist(rgamma(n = n*1000, shape = paramhat_sir_gamma$estimate[1], rate = paramhat_sir_gamma$estimate[2]), freq = FALSE, col = "skyblue",
     xlab = "Rayon", ylab = "Densité", main = "")
lines(density(sir), col = "red", lwd = 3)

hist(rweibull(n = n, shape = paramhat_sir_weibull$estimate[1], scale = paramhat_sir_weibull$estimate[2]), freq = FALSE, col = "pink3")
lines(density(sir))

# AIC/BIC validation

AIC(paramhat_sir_gamma, paramhat_sir_normal, paramhat_sir_lognorm, paramhat_sir_weibull)
BIC(paramhat_sir_gamma, paramhat_sir_normal, paramhat_sir_lognorm, paramhat_sir_weibull)



# ------------------------------------------------------------------------------
# ESTIMATIONS
# ------------------------------------------------------------------------------
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))
dexponential <- function(x, l) l*exp(-l*x)
pexponential <- function(p, l) 1-exp(-l*(1-p))
# pressu

fit.gamma <- fitdist(pressu, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)
fit.gamma
fit.weibull <- fitdist(pressu, distr = "weibull", method = "mle")
fit.lnorm <- fitdist(pressu, distr = "lnorm", method = "mle")
  # fit.gumbel <- fitdist(pressu, "gumbel", start=list(a=10, b=10))
par(mfrow = c(1, 1))

plot.legend <- c("Weibull", "Gamma", "Lognorm")
denscomp(list(fit.weibull, fit.gamma, fit.lnorm), fitcol = c("red", "red", "red"), legendtext = plot.legend, main = "Répartition de la pression",
         datacol = "deepskyblue3", ylim = c(0,0.033), xlab="", ylab="")

points(density(pressu, kernel = "gaussian"),type = "l", col= "chartreuse")

# windspeed

fit.gamma <- fitdist(windspeed, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)
fit.gamma
fit.weibull <- fitdist(windspeed, distr = "weibull", method = "mle")
fit.lnorm <- fitdist(windspeed, distr = "lnorm", method = "mle")

par(mfrow = c(1, 1))
plot.legend <- c("Weibull", "Gamma", "Lognorm")
denscomp(list(fit.weibull, fit.gamma, fit.lnorm), fitcol = c("red", "red", "red"), legendtext = plot.legend, main = "Répartition de la vitesse du vent",
         datacol = "deepskyblue3", ylim = c(0,0.05), xlab="", ylab="")

points(density(windspeed, kernel = "gaussian"),type = "l", col= "chartreuse")

# sir

hist(sir, breaks =200)
plot(density(sir))

fit.gamma <- fitdist(sir, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)
fit.gamma
fit.weibull <- fitdist(sir, distr = "weibull", method = "mle")
fit.lnorm <- fitdist(sir, distr = "lnorm", method = "mle")
par(mfrow = c(1, 1))
plot.legend <- c("Weibull", "Gamma", "Lognorm")
denscomp(list(fit.weibull, fit.gamma, fit.lnorm), fitcol = c("red", "red", "red"), legendtext = plot.legend, main = "Répartition du rayon",
         datacol = "deepskyblue3", ylim = c(0,0.01), xlab="", ylab="")

points(density(sir, kernel = "gaussian"),type = "l", col= "chartreuse")


# ------------------------------------------------------------------------------
# SEMI PARAMETRIC METHOD - CML (Canonical Maximum Likelihood)
# ------------------------------------------------------------------------------

  # Pressu and windspeed : 2 candidates

selectedCopula <- BiCopSelect(pobs(windspeed), pobs(pressu), family = 1:6, selectioncrit = "logLik", method = "mle")
summary(selectedCopula)

selectedCopula <- BiCopSelect(pobs(windspeed), pobs(pressu), family = 1:10, selectioncrit = "logLik", method = "mle")
summary(selectedCopula)

    ## Compute the MPLE and its standard error

rotGcop_PRS_WND = fitCopula(rotCopula(gumbelCopula(), flip=c(FALSE,TRUE)),
                            data = pobs(dat[,c(2,1)]),
                            method = "mpl")
summary(rotGcop_PRS_WND)

rotGcop_PRS_WND = fitCopula(rotCopula(gumbelCopula(), flip=c(FALSE,TRUE)),
                            data = pobs(dat[,c(2,1)]),
                            method = "itau")
summary(rotGcop_PRS_WND)

rotGcop_PRS_WND = fitCopula(r270GumbelCopula(),
                            #rotCopula(gumbelCopula(), flip=c(FALSE,TRUE)),
                            data = pobs(dat[,c(2,1)]),
                            method = "irho")
summary(rotGcop_PRS_WND)


rotBcop_PRS_WND = fitCopula(r270BB6Copula(),
                            data = pobs(dat[,c(2,1)]),
                            method = "ml")
summary(rotBcop_PRS_WND)

## Perspective presentation of selected copulas
par(mfrow=c(1,1))
persp(rotCopula(gumbelCopula(dim=2, 7.0), flip=c(FALSE,TRUE)),dCopula,col='deepskyblue4', theta = 310, phi = 25)
;persp(rotCopula(BB6Copula(param = c(1.6, 5.2)), flip=c(FALSE,TRUE)),dCopula,col='deepskyblue4', theta = 310, phi = 25)

fc1 <- r270GumbelCopula(-7) 
wireframe2(fc1, FUN = dCopula, delta = 0.00001, zlim = 0:90, col ="deepskyblue4", n.grid =40, xlab="WND", ylab="PRS", zlab ="") # wireframe plot (density)

fc2 <- r270BB6Copula(param = c(-1.6,-5.2)) 
wireframe2(fc2, FUN = dCopula, delta = 0.00001, zlim = 0:90, col ="deepskyblue4", n.grid =40, xlab="WND", ylab="PRS", zlab ="") # wireframe plot (density)

contourplot2(fc1, FUN = dCopula, n.grid = 72, # contour plot (density)
             lwd = 1/2, xlab ="WND", ylab ="PRS",
             col.regions = colorRampPalette(c("white", "deepskyblue4")), col = "deepskyblue4")

contourplot2(fc2, FUN = dCopula, n.grid = 72, # contour plot (density)
             lwd = 1/2, xlab ="WND", ylab ="PRS",
             col.regions = colorRampPalette(c("white", "deepskyblue4")), col = "deepskyblue4")


  # Pressu and radius: 2 candidates
selectedCopula <- BiCopSelect(pobs(sir), pobs(pressu), family = 1:6, selectioncrit = "logLik")
summary(selectedCopula)
selectedCopula <- BiCopSelect(pobs(sir), pobs(pressu), family = 1:10, selectioncrit = "logLik")
summary(selectedCopula)


    ## Compute the MPLE and its standard error

cop_PRS_SiR34 = fitCopula(frankCopula(), data = pobs(dat[,c(3,1)]), method = "mpl")
summary(cop_PRS_SiR34)

cop_PRS_SiR34 = fitCopula(frankCopula(), data = pobs(dat[,c(3,1)]), method = "itau")
summary(cop_PRS_SiR34)

cop_PRS_SiR34 = fitCopula(frankCopula(), data = pobs(dat[,c(3,1)]), method = "irho")
summary(cop_PRS_SiR34)



rotBcop_PRS_SiR34 = fitCopula(r270BB8Copula(),
                            data = pobs(dat[,c(3,1)]),
                            method = "ml")
summary(rotBcop_PRS_SiR34)


  ## Perspective presentation of selected copulas

persp(rotCopula(BB8Copula(param = c(6, 0.62)), flip=c(FALSE,TRUE)),dCopula,col='deepskyblue4', theta = 310, phi = 25)
persp(frankCopula(dim=2, 5.3), dCopula,col='deepskyblue4', theta = 40, phi = 25)

fc1 <- frankCopula(-5.4) # define a Frank copula
wireframe2(fc1, FUN = dCopula, delta = 0.00001, zlim = 0:8, col ="deepskyblue4", n.grid =40, xlab="SiR34", ylab="PRS", zlab ="") # wireframe plot (density)
contourplot2(fc1, FUN = dCopula, n.grid = 72, # contour plot (density)
             lwd = 1/2, xlab ="SiR34", ylab ="PRS",
             col.regions = colorRampPalette(c("white", "deepskyblue4")), col = "deepskyblue4")

fc2 <- r270BB8Copula(param = c(-6,-0.62)) # define a Frank copula
wireframe2(fc2, FUN = dCopula, delta = 0.00001, zlim = 0:8, col ="deepskyblue4", n.grid =40, xlab="SiR34", ylab="PRS", zlab ="") # wireframe plot (density)
contourplot2(fc2, FUN = dCopula, n.grid = 72, # contour plot (density)
             lwd = 1/2, xlab ="SiR34", ylab ="PRS",
             col.regions = colorRampPalette(c("white", "deepskyblue4")), col = "deepskyblue4")

# ------------------------------------------------------------------------------
# TESTS D'INDEPENDANCE
# ------------------------------------------------------------------------------

R <- pobs(windspeed)
S <- pobs(sir)
P <- pobs(pressu)

  # ------------------------------------------------------------------------------
  # Test parametrique: pearson (pas adapte ici)
  # ------------------------------------------------------------------------------

cor(dat, method = "pearson") # Pas adapte car donnees ne sont pas distribuees comme une normale
cor.test(dat$WND, dat$PRS,  method="pearson")
  
  # ------------------------------------------------------------------------------
  # Test non parametriques: kendall, spearman
  # ------------------------------------------------------------------------------

# Correlation
cor(dat, method = "kendall")
cor(dat, method = "spearman")

cor.test(dat$WND, dat$PRS,  method="kendall")
cor.test(dat$SiR34, dat$PRS,  method="kendall")

cor.test(dat$WND, dat$PRS,  method = "spearman", exact = FALSE)
cor.test(dat$SiR34, dat$PRS,  method = "spearman", exact = FALSE)

  # ------------------------------------------------------------------------------
  # Independence tests (Kendall's test)
  # ------------------------------------------------------------------------------

BiCopIndTest(R, P) # p.value equals 0 => rejecting the H0 independence hypothesis
BiCopIndTest(S, P) # same
# BiCopIndTest(R, S)

# ------------------------------------------------------------------------------
# TESTS D' ADEQUATION DES COPULES
# ------------------------------------------------------------------------------

# Goodness-of-fit tests
  # Test Cramer Von Mises et Kolmogorov Smirnov
  # Avec bootstrap

BiCopGofTest(R, P, family = 34, method = "Kendall", B = 10) # 34 rotated Gumbel
BiCopGofTest(R, P, family = 38, method = "Kendall", B = 10) # 38 rotated BB6

P2 <- 1-P # Flipping for rotated Frank ...
BiCopGofTest(S, P2, family = 5, method = "Kendall", B = 10) # Frank on rotated data
BiCopGofTest(S, P, family = 40, method = "Kendall", B = 10) # rotated BB8
