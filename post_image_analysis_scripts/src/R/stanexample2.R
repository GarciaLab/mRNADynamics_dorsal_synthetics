datTemp <- read.csv('C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\dugongdata.csv')

dat = list("N" = length(datTemp$x), "x" = datTemp$x, "Y"=datTemp$Y)

nlm <- nls(Y ~ alpha - beta * lambda^x, data=dat,
           start=list(alpha=1, beta=1, lambda=0.9))
summary(nlm)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit <- stan(file = "C:\\Users\\owner\\Desktop\\LivemRNA\\mRNADynamics\\test.stan", 
            model_name = "GrowthCurve", 
            data = dat)
print(fit, c("alpha", "beta", "lambda", "sigma"))


Y_mean <- extract(fit, "Y_mean")
Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)

Y_pred <- extract(fit, "Y_pred")
Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)

plot(dat$Y ~ dat$x, xlab="x", ylab="Y", 
     ylim=c(1.6, 2.8), main="Non-linear Growth Curve")
lines(dat$x, Y_mean_mean)
points(dat$x, Y_pred_mean, pch=19)
lines(dat$x, Y_mean_cred[1,], col=4)
lines(dat$x, Y_mean_cred[2,], col=4)
lines(dat$x, Y_pred_cred[1,], col=2)
lines(dat$x, Y_pred_cred[2,], col=2)
legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
       legend=c("observation", "prediction", "mean prediction",
                "90% mean cred. interval", "90% pred. cred. interval"),
       col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))
# 
# library(bayesplot)
# library(tidyverse)
# 
# fit%>%
#   mcmc_trace()
# 
# fit %>%
#   rhat() %>%
#   mcmc_rhat() +
#   yaxis_text()
