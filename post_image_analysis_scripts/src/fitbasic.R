library(pracma)
getFrac <- function(x, c, kd, tcycle) {
  n <- 5
  out = c();
  for (i in 1:length(x)) {
  arg = (c*x[i]*tcycle /(x[i]+kd))
  out[i] <- 1 - (incgam(arg, n)/gamma(n));
  }
  return(out)
}

getOnset <- function(x, c, kd, tcycle) {
  n <- 5
  arg = (c*x*tcycle /(x+kd))
  out <- ((x+kd)/(c*x))*(gamma(n+1) - incgam(n+1, arg))/(gamma(n) - incgam(n, arg));
  return(out)
}


datTemp <-
  read.csv(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\basic.csv")

binned <-
  read.csv(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\binned.csv")

x = datTemp$x;
fraction = datTemp$fraction;
onset = datTemp$onset;

theta0 = list("p0"=c(1, 1E3, 8), "lb"=c(1E-2, 1E2, 4), "ub"=c(1E2, 1E5, 9))
  
dat = list(
  "N" = length(x),
  "x" = x,
  "fraction" = fraction,
  "onset" = onset,
  "p0" = theta0$p0,
  "lb" = theta0$lb,
  "ub" = theta0$ub,
  "K" = length(theta0$p0)
)

df <- data.frame(x = x, fraction = fraction, onset=onset)
#nlm <- nls(Y ~ (w*x/KD)/(1+(x/KD)+(w*x/KD)), data=dat,
#           start=list(w=1, KD=500))
#summary(nlm)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit <-
  stan(file = "C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\src\\fitbasic_ode.stan",
       model_name = "fitbasic_ode",
       data = dat)
print(fit, c("c", "kd", "tcycle","sigma_fraction", "sigma_onset"))


fraction_mean <- extract(fit, "fraction_mean")
fraction_mean_cred <- apply(fraction_mean$fraction_mean, 2, quantile, c(0.05, 0.95))
fraction_mean_mean <- apply(fraction_mean$fraction_mean, 2, mean)


onset_mean <- extract(fit, "onset_mean")
onset_mean_cred <- apply(onset_mean$onset_mean, 2, quantile, c(0.05, 0.95))
onset_mean_mean <- apply(onset_mean$onset_mean, 2, mean)


fraction_pred <- extract(fit, "fraction_pred")
fraction_pred_cred <- apply(fraction_pred$fraction_pred, 2, quantile, c(0.05, 0.95))
fraction_pred_mean <- apply(fraction_pred$fraction_pred, 2, mean)


onset_pred <- extract(fit, "onset_pred")
onset_pred_cred <- apply(onset_pred$onset_pred, 2, quantile, c(0.05, 0.95))
onset_pred_mean <- apply(onset_pred$onset_pred, 2, mean)


predicted_df <- data.frame(
  fraction_mean_mean = fraction_mean_mean,
  fraction_pred_mean = fraction_pred_mean,
  fraction_mean_cred_lower = fraction_mean_cred[1,],
  fraction_mean_cred_upper = fraction_mean_cred[2,],
  fraction_pred_cred_lower = fraction_pred_cred[1,],
  fraction_pred_cred_upper = fraction_pred_cred[2,],
  onset_mean_mean = onset_mean_mean,
  onset_pred_mean = onset_pred_mean,
  onset_mean_cred_lower = onset_mean_cred[1,],
  onset_mean_cred_upper = onset_mean_cred[2,],
  onset_pred_cred_lower = onset_pred_cred[1,],
  onset_pred_cred_upper = onset_pred_cred[2,],
  x = x
)

z = binned$binMidValues
yy = mean(fraction)
zz = getFrac(z, 1, 595.06, 8.36)
plot(z, zz)
lines(z, binned$frac_mean)

plot(x, fraction_mean_mean)
plot(x, onset_mean_mean)


plot(dat$fraction ~ dat$x, xlab="x", ylab="fraction", 
     ylim=c(0, 1), main="fraction")
lines(dat$x, fraction_mean_mean)
points(dat$x, fraction_pred_mean, pch=19)
lines(dat$x, fraction_mean_cred[1,], col=4)
lines(dat$x, fraction_mean_cred[2,], col=4)
lines(dat$x, fraction_pred_cred[1,], col=2)
lines(dat$x, fraction_pred_cred[2,], col=2)
legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
       legend=c("observation", "prediction", "mean prediction",
                "90% mean cred. interval", "90% pred. cred. interval"),
       col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))
# 



library(extrafont)
library(viridis)
library(hrbrthemes)
library(ggthemes)
#font_import()
#loadfonts(device = "win")
#windowsFonts("Arial" = windowsFont("Arial"))

dfnames = c()

namesList = c("1Dg", "1DgS", "1DgW", "1DgAW3", "1DgSVW2", "1DgVVW", "1DgVW")

for (j in 1:nx) {
  dfnames[j] = namesList[dsid[j]]
}

df$dfnames = dfnames

predicted_df$dfnames = dfnames


windows()

df %>%
  ggplot(aes(x = x, y = Y, group = dfnames)) +
  # geom_point() +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(size = 8, family = "ArialMT"),
    plot.title = element_text(size = 13, family = "ArialMT"),
    text = element_text(family = "ArialMT")
  ) +
  xlab("[Dorsal] (au)") +
  ylab("max fluorescence (au)") +
  facet_wrap( ~ dfnames) +
  geom_line(data = predicted_df,
            aes(x = x, y = Y_mean_mean, group = dfnames)) +
  geom_line(
    data = predicted_df,
    aes(x = x, y = Y_pred_cred_lower, group = dfnames),
    linetype = "dashed",
    color = "red"
  ) +
  geom_line(
    data = predicted_df,
    aes(x = x, y = Y_pred_cred_upper, group = dfnames),
    linetype = "dashed",
    color = "red"
  ) +
  geom_ribbon(
    data = predicted_df,
    aes(ymin = Y_mean_cred_lower, ymax = Y_mean_cred_upper, group = dfnames),
    fill = "blue",
    alpha = 0.3
  )
#legend=c("observation", "prediction", "mean prediction",
#        "90% mean cred. interval", "90% pred. cred. interval")

ggsave('C:\\Users\\owner\\Dropbox\\DorsalSyntheticsDropbox\\manuscript\\ggplot.pdf')


library(bayesplot)
posteriormat <- as.matrix(fit)
posteriorarr <- as.array(fit)

windows()

mcmc_pairs(posteriorarr,
           pars = c("w", "R", "KD1", "KD3", "KD4", "sigma"))




plot(
  dat$Y ~ dat$x,
  xlab = "x",
  ylab = "Y",
  main = "simple weak fraction active",
  col = "green",
  bg = "green"
)
lines(dat$x, Y_mean_mean)
points(dat$x, Y_pred_mean, pch = 19)
lines(dat$x, Y_mean_cred[1, ], col = 4)
lines(dat$x, Y_mean_cred[2, ], col = 4)
lines(dat$x, Y_pred_cred[1, ], col = 2)
lines(dat$x, Y_pred_cred[2, ], col = 2)
legend(
  x = "bottomright",
  bty = "n",
  lwd = 2,
  lty = c(NA, NA, 1, 1, 1),
  legend = c(
    "observation",
    "prediction",
    "mean prediction",
    "90% mean cred. interval",
    "90% pred. cred. interval"
  ),
  col = c(1, 1, 1, 4, 2),
  pch = c(1, 19, NA, NA, NA)
)
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

np <- nuts_params(fit)
color_scheme_set("red")
mcmc_intervals(posteriormat,
               pars = c("w", "KD1", "KD2",
                        "KD3", "KD4", "KD5",
                        "KD6", "KD7", "R", "sigma"))
windows()

(p <- mcmc_pairs(
  posteriorarr,
  pars = c("w", "KD1", "KD2",
           "KD3", "KD4", "KD5",
           "KD6", "KD7", "R", "sigma")
))




mcmc_scatter(posteriorarr, pars = c("KD7", "R"))

windows()
mcmc_trace(posteriorarr,
           pars = c("w", "KD1", "KD2",
                    "KD3", "KD4", "KD5",
                    "KD6", "KD7", "R", "sigma"))


xaxis_title(size = 13, family = "sans") +
  xaxis_ticks(on = FALSE) +
  yaxis_ticks(on = FALSE)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
windows()
mcmc_parcoord(
  posteriormat,
  transformations = "log",
  pars = c("w", "KD1", "KD2",
           "KD3", "KD4", "KD5",
           "KD6", "KD7", "R", "sigma")
)



windows()
color_scheme_set("darkgray")
mcmc_parcoord(
  posteriormat,
  transformations = "log",
  pars = c("w", "R"),
  np = np
)

mcmc_areas(
  posterior,
  pars = c("w", "KD1", "KD2",
           "KD3", "KD4", "KD5",
           "KD6", "KD7", "R", "sigma"),
  prob = 0.8,
  transformations = "log"
) + plot_title

# mcmc_areas(posterior,
#            prob = 0.8) + plot_title
