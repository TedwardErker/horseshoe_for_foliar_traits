library(rstanarm)
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

options(mc.cores = parallel::detectCores())

lma <- read_excel("data/FFT_LMA.xlsx")
colnames(lma) <- c("sample", "lma") #renames columns to be same as "x" names and simplify them
x <- read_csv("data/FFT_Spectra_and_NIT_broadleaf.csv")
cn <- colnames(x)
colnames(x) <- ifelse(!grepl("^[0-9].*",cn), cn, paste0("X",cn)) # rename columns so they don't start with numbers
x <- select(x, -nitrogen_percent) # drop nitrogen, not our target
d <- left_join(lma, x)# join data

dc <- dplyr::filter(d, complete.cases(d))

head(dc)

dcr <- dc %>% sample_n(., 30) %>% gather(key = wv, value = refl, -sample, -lma)
dcr <- dcr %>% rowwise() %>% mutate(wv = as.numeric(str_sub(wv, 2, nchar(wv))))

plt <- ggplot(dcr, aes(x = wv, y = refl, color = lma, group = sample)) + geom_line() + theme(legend.position = c(.85,.85))
print(plt)

#library(plotly)
#ggplotly(p, dynamicTicks = T)

library(corrplot)
a <- spread(dcr, key = wv, value = refl)
colnames(a)[-2] <- ifelse(nchar(colnames(a)[-2]) == 3, paste0("0",colnames(a)[-2]), colnames(a)[-2])  # to make alphabetic ordering work right
a[,c(2,seq(3,ncol(a),20))] %>% cor %>% corrplot(order = "alphabet", type = "upper", method = "color")

dim(dc)

n <- 50 # number of observations to sample
p <- 80 # number of wavelengths to sample
rows <- sample(1:nrow(dc),n)
cols <- sample(3:ncol(dc),p)
dcs <- dc[rows, c(2,cols)]  # make sure I keep column 2, which contains the target response, lma.

dim(dcs)

lma_center <- mean(dcs$lma)
dcs$lma <- dcs$lma - lma_center

p0 <- 5 # prior guess for the number of non zero coefficients
  sigmaguess <- 5  # guess for sigma of gaussian response
  tau0 <- p0 / (p - p0) * sigmaguess/sqrt(n)
tau0

fit <- stan_glm(lma ~ ., data = dcs, gaussian(), prior = hs(global_scale=tau0), prior_intercept = normal())

plot(fit)

unsampledrows <- (1:nrow(dc))[!1:nrow(dc) %in% rows]
newrows <- sample(unsampledrows, 50)
newdata <- dc[newrows, c(2,cols)]
y_rep <- posterior_predict(fit, newdata)
predictions_mean <- apply(y_rep, 2, mean) + lma_center

sqrt(mean((newdata$lma - predictions_mean)^2))

fit_gaus <- stan_glm(lma ~ ., data = dcs, gaussian(), prior = normal())

plot(fit_gaus)

y_rep_gaus <- posterior_predict(fit_gaus, newdata)
predictions_mean_gaus <- apply(y_rep_gaus, 2, mean) + lma_center

sqrt(mean((newdata$lma - predictions_mean_gaus)^2))

dcss <- dcs[,1:5]
m <- lm(lma ~ ., dcss)

summary(m)

lm_newdata <- dc[newrows,colnames(dcss)]
lm_predict <- predict(m, lm_newdata)

sqrt(mean((lm_newdata$lma - lm_predict)^2))
