#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2022-07-11 19:24
# @Description: appendix table 2: model fitting and parameter selection for growth models
#--------------------------------------------------------------------------#
# clean cp memory
cat("\014")
rm(list = ls())
gc()

# create a work dictionary
# replace the file path with the one you used
setwd("D:/我的坚果云/002项目/王老师合作项目/202207_NEW_RESULT")
getwd()

# load packages
library(dplyr)
library(readr)

# import the data of exp. 1
Mdata <- read_csv("./data_Nmodel_2022.csv")
Mdata <- data.frame(Mdata)

# an overview of the data
View(Mdata)
str(Mdata)
glimpse(Mdata)

# --------------------------------------------------------------------------
# step 1: model preparation
# --------------------------------------------------------------------------
# six candidate growth models
types <- c("liner", "power", "expontial", "logistic.3p", "logistic.4p", "gompertz")

# initial mass of parental plants
M0 <- 22.38

# ==========================================================================
# function of model fitting
# ==========================================================================
model = function(types, time, M0, r, N, gamma, delta, sigma, L){
					if(missing(N)) N = 0
					if(missing(gamma)) gamma = 0
					if(missing(delta)) delta = 0
					if(missing(sigma)) sigma = 0
					if(missing(L)) L = 0

					K <- gamma * N + delta
					if(missing(K)) K = 0

					try ({
						if(types == "liner"){out <- M0 + r * time}
						if(types == "expontial"){out <- M0 * exp(r * time)}
						if(types == "power"){out <- (M0^(1 - sigma) + r * time * (1 - sigma))^(1/(1 - sigma))}
						if(types == "monomolecular"){out <- K - (K - M0) * exp(-r * time)}
						if(types == "logistic.3p"){out <- M0 * K / (M0 + (K - M0) * exp (-r * time))}
						if(types == "logistic.4p"){out <- L + M0 * (K - L)/ (M0 + (K - L - M0) * exp (-r * time))}
						if(types == "gompertz"){out <- K * (M0 / K)^exp(- r * time)}
					}, silent = TRUE)

					return(out)}

# ==========================================================================
# function of model fitting
# ==========================================================================
# --------------------------------------------------------------------------
# step 2: parameter setting
# --------------------------------------------------------------------------
# "liner" model, starting at r = 91.33
model_out1 <- nls(tmass ~ model(types = "liner", time = htime, M0 = M0, r, N = nlevel),
                  data = Mdata,
                  start = list(r = 91.33),
                  algorithm = "port",
                  lower = c(0))

(model_out1)

# "expontial" model, starting at r = 0.08132
model_out2 <- nls(tmass ~ model(types = "expontial", time = htime, M0 = M0, r, N = nlevel),
                 data = Mdata,
                 start = list(r = 0.08132),
                 algorithm = "port",
                 lower = c(0))
(model_out2)

# "power" model, starting at r = 0.000000005652916657 and sigma = 75
# some unsolved problems in "power" mdoel
model_out3 <- nls(tmass ~ model(types = "power", time = htime, M0 = M0, r, N = nlevel),
                 data = Mdata,
                 start = list(r = 0.000000005652916657, sigma = 75),
                 algorithm = "port",
                 lower = c(0),
                 control = list(maxiter = 1000, trace = TRUE, warnOnly = TRUE))
(model_out3)

# "monomolecular" model, starting at r = 0.000000005652916657, gamma = 171538365.926158, and delta = 7576293277.62156
# some unsolved problems in "monomolecular" model
model_out4 <- nls(tmass ~ model(types = "monomolecular", time = htime, M0 = M0, r, N = nlevel, gamma, delta),
                 data = Mdata,
                 start = list(r = 0.000000005652916657, gamma = 171538365.926158, delta = 7576293277.62156),
                 algorithm = "port",
                 lower = c(0),
                 control = list(maxiter = 1000, trace = TRUE, warnOnly = TRUE))
(model_out4)

# "logistic.3p" model, starting at r = 0.1061, gamma = 163.4468, and delta = 2946.6494
model_out5 <- nls(tmass ~ model(types = "logistic.3p", time = htime, M0 = M0, r, N = nlevel, gamma, delta),
                 data = Mdata,
                 start = list(r = 0.1061, gamma = 163.4468, delta = 2946.6494),
                 algorithm = "port",
                 lower = c(0, 0, 0))
(model_out5)

# "logistic.4p" model, starting at r = 0.1061, gamma = 163.4468, delta = 2946.6493, and L = 0
model_out6 <- nls(tmass ~ model(types = "logistic.4p", time = htime, M0 = M0, r, N = nlevel, gamma, delta, sigma, L),
                 data = Mdata,
                 start = list(r = 0.1061, gamma = 163.4468, delta = 2946.6493, L = 0),
                 algorithm = "port",
                 lower = c(0, 0, 0, 0))
(model_out6)

# "gompertz" model, starting at r = 0.01976, gamma = 844.5, and delta = 14890
model_out7 <- nls(tmass ~ model(types = "gompertz", time = htime, M0 = M0, r, N = nlevel, gamma, delta),
                 data = Mdata,
                 start = list(r = 0.01976, gamma = 844.5, delta = 14890),
                 algorithm = "port",
                 lower = c(0, 0, 0))
(model_out7)

# --------------------------------------------------------------------------
# step 3: calculate AICs and R squared
# --------------------------------------------------------------------------
# "liner" model
fit.model_out1 <- predict(model_out1)
cor_value1     <- cor(Mdata$tmass, fit.model_out1)
SSE1           <- sum((Mdata$tmass - fit.model_out1)^2)
SST1           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_1     <- 1 - SSE1/SST1
aic1           <- AIC(model_out1)

# "expontial" model
fit.model_out2 <- predict(model_out2)
cor_value2     <- cor(Mdata$tmass, fit.model_out2)
SSE2           <- sum((Mdata$tmass - fit.model_out1)^2)
SST2           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_2     <- 1 - SSE2/SST2
aic2           <- AIC(model_out2)

# "power" model
fit.model_out3 <- predict(model_out3)
cor_value3     <- cor(Mdata$tmass, fit.model_out3)
SSE3           <- sum((Mdata$tmass - fit.model_out3)^2)
SST3           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_3     <- 1 - SSE3/SST3
aic3           <- AIC(model_out3)

# "monomolecular" model
fit.model_out4 <- predict(model_out4)
cor_value4     <- cor(Mdata$tmass, fit.model_out4)
SSE4           <- sum((Mdata$tmass - fit.model_out4)^2)
SST4           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_4     <- 1 - SSE4/SST4
aic4           <- AIC(model_out4)

# "logistic.3p" model
fit.model_out5 <- predict(model_out5)
cor_value5     <- cor(Mdata$tmass, fit.model_out5)
SSE5           <- sum((Mdata$tmass - fit.model_out5)^2)
SST5           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_5     <- 1 - SSE5/SST5
aic5           <- AIC(model_out5)

# "logistic.4p" model
fit.model_out6 <- predict(model_out6)
cor_value6     <- cor(Mdata$tmass, fit.model_out6)
SSE6           <- sum((Mdata$tmass - fit.model_out6)^2)
SST6           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_6     <- 1-SSE6/SST6
aic6           <- AIC(model_out6)

# "gompertz" model
fit.model_out7 <- predict(model_out7)
cor_value7     <- cor(Mdata$tmass, fit.model_out7)
SSE7           <- sum((Mdata$tmass - fit.model_out7)^2)
SST7           <- sum((Mdata$tmass - mean(Mdata$tmass))^2)
Rsquared_7     <- 1 - SSE7/SST7
aic7           <- AIC(model_out7)

# combine all AICs and r squaredd
cor_values <- c(cor_value1, cor_value2, cor_value4, cor_value5, cor_value6, cor_value7)
aics       <- c(aic1, aic2, aic4, aic5, aic6, aic7)
Rsquared   <- c(Rsquared_1, Rsquared_2, Rsquared_4, Rsquared_5, Rsquared_6, Rsquared_7)
gofit      <- rbind(cor_values, aics, Rsquared)

# --------------------------------------------------------------------------
# step 4: obtaining model coefficients
# --------------------------------------------------------------------------
out_pa1 <- summary(model_out1)$coefficients[, 1]
out_pa2 <- summary(model_out2)$coefficients[, 1]
out_pa3 <- summary(model_out3)$coefficients[, 1]
out_pa4 <- summary(model_out4)$coefficients[, 1]
out_pa5 <- summary(model_out5)$coefficients[, 1]
out_pa6 <- summary(model_out6)$coefficients[, 1]
out_pa7 <- summary(model_out7)$coefficients[, 1]

# list of model coefficients
a_list <- list(out_pa1, out_pa2, out_pa4, out_pa5, out_pa6, out_pa7)

# combine model coefficients
out1 <- do.call(cbind, lapply(lapply(a_list, unlist), `length<-`, max(lengths(a_list))))
out2 <- rbind(out1, gofit)

# rename the rows
par_name <- c("r", "gamma", "delta", "L", "cor_values", "AICs", "Rsquared")
out      <- cbind(par_name, out2)

# --------------------------------------------------------------------------
# step 5: exporting model coefficients
# --------------------------------------------------------------------------
out <- as.data.frame(out, stringsAsFactors = F)
colnames(out) <- c("model_parameter", "liner", "expontial", "monomolecular", "logistic_3p", "logistic_4p", "gompertz")
head(out)

# exporting results - Appendix Table 2
write_excel_csv(out,"./appendix.table.2_growth_model_output.csv")


#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2022-07-11 19:24
# @Description:  appendix fig. 2: correlation between number (and mean mass) of propagules and total mass of paretnal plants
#--------------------------------------------------------------------------#
# clean cp memory
cat("\014")
rm(list = ls())
gc()

# create a work dictionary
# replace the file path with the one you used
setwd("D:/我的坚果云/002项目/王老师合作项目/202207_NEW_RESULT")
getwd()

# load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(readr)

# import the data of exp. 1
Mdata <- read_csv("./data_Nmodel_2022.csv")
Mdata <- data.frame(Mdata)

# an overview of the data
View(Mdata)
str(Mdata)
glimpse(Mdata)

# --------------------------------------------------------------------------
# step 1: nonlinear regression
# --------------------------------------------------------------------------
# ==========================================================================
# function for calculating R^2, COR, and AIC of non-linear mdoels
# ==========================================================================
fit.out = function(fit_model, x, y){
									# summary values
										summary.fit <- summary(fit_model)

									# predicted values
										y_pred <- predict(fit_model)

									# create plot using real and simulated data
										plot(x, y)
										lines(x, y_pred, lty = 2, col = "red", lwd = 3)

									# caculate the "cor" of y
										cor_y <- cor(y, y_pred)

									# caculate the "pesudo R squared" of y
									# caculauted R^2 "return 1 - ((y_test - y_true)**2).sum() / ((y_true - y_true.mean())**2).sum()"
										SSE <- sum((y-y_pred)^2)
										SST <- sum((y-mean(y))^2)
										pesudo_R2 <- 1-SSE/SST

									# caculate the "AIC
										AIC_y <- AIC(fit_model)

									# output of results

										cat("\n", "[1] User-defined fomulate:", "\n")
										print(summary.fit$formula)

										cat("\n", "[2] Coefficients and P values as listed below:", "\n")
										print(summary.fit$coefficients)

										cat("\n", "[3] Test goodness of fit:", "\n")
										cat(paste("pesudo_R2 = ", pesudo_R2), "\n", paste("Cor_y = ", cor_y), "\n", paste("AIC_y = ", AIC_y))
									}

# ==========================================================================
# function for calculating R2,COR, and AIC of non-linear mdoels
# ==========================================================================

# regression between number of propagule and total mass of parental plants
model_out0 <- nls(propagule_num ~ gamma*(tmass)^sigma,
                  data = Mdata, start = list(gamma = 0.4814, sigma = 0.6093))

fit.out(model_out0, Mdata$tmass, Mdata$propagule_num)

# regression between mean mass of propagule and total mass of parental plants
model_out1 <- nls(mean_propagule_mass ~ gamma*(tmass)^sigma,
                 data = Mdata, start = list(gamma = 2.727, sigma = 0.355))

fit.out(model_out1, Mdata$tmass, Mdata$mean_propagule_mass)


# --------------------------------------------------------------------------
# step 2: separate plots for two regressions
# --------------------------------------------------------------------------
# correlation between number of propagules and total mass of parental plants
plot_propagule_num <- ggplot(Mdata, aes(tmass, propagule_num)) +
				             geom_point(col = 'grey', size = 4) +
				             geom_line(aes(tmass, fitted(model_out0)), col = 'blue', size = 1) +
							 theme_classic() +
							 theme(legend.position = "none",
								   axis.line = element_line(size = 1, linetype = "solid"),
			                       axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
			                       axis.text = element_text(family = "serif", colour = "black", size = 14),
			                       axis.title = element_text(family = "serif", colour = "black", size = 14)) +
				             scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50)) +
				             scale_x_continuous(limits = c(0, 18000), breaks = seq(0, 15000, by = 5000)) +
				             labs(y = "Number of clonal propagules", x = "Total mass of one parental plant (mg)")

plot_propagule_num

eq_propagule_num <- substitute(plain(y) == gamma %.% plain(x)^sigma *", "~plain(R)^2~" = "~0.91*", "~italic(P)~"< "~0.05, list(gamma = 0.253, sigma = 0.684))
plot_propagule_num <- plot_propagule_num + geom_text(aes(x = 10000, y = 25, family = "serif", size = 14, label = as.character(as.expression(eq_propagule_num))), parse = TRUE)
plot_propagule_num

# correlation between mean mass of propagules and total mass of parental plants
plot_propagule_mean.mass <- ggplot(Mdata, aes(tmass, mean_propagule_mass)) +
								   geom_point(col = 'grey', size = 4) +
								   geom_line(aes(tmass, fitted(model_out1)), col = 'blue', size = 1) +
								   theme_classic() +
								   theme(legend.position = "none",
										 axis.line = element_line(size = 1, linetype = "solid"),
				                         axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
				                         axis.text = element_text(family = "serif", colour = "black", size = 14),
				                         axis.title = element_text(family = "serif", colour = "black", size = 14)) +
					               scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
					               scale_x_continuous(limits = c(0, 18000), breaks = seq(0, 15000, by = 5000)) +
					               labs(y = "Mean mass of clonal propagules (mg)", x = "Total mass of one parental plant (mg)")

plot_propagule_mean.mass

eq_propagule_mean.mass <- substitute(plain(y) == gamma %.% plain(x)^sigma*", "~plain(R)^2~" = "~0.77*", "~italic(P)~"< "~0.05, list(gamma = 2.727, sigma = 0.354))
plot_propagule_mean.mass <- plot_propagule_mean.mass + geom_text(aes(x = 10000, y = 10, family = "serif", size = 14, label = as.character(as.expression(eq_propagule_mean.mass))), parse = TRUE)
plot_propagule_mean.mass

# --------------------------------------------------------------------------
# step 3: combine two separate plots
# --------------------------------------------------------------------------
# combined figs
plots_correlation <- plot_propagule_num + plot_propagule_mean.mass + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')

# export the combined plot
ggexport(plots_correlation, filename = "./appendix.fig.2_correlation.png",
         width = 3000,
         height = 1500,
         pointsize = 12,
         res = 300)

#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2022-07-11 19:24
# @Description: model fitting and parameter selection for size distribution models
#--------------------------------------------------------------------------#
# clean cp memory
cat("\014")
rm(list = ls())
gc()

# create a work dictionary
# replace the file path with the one you used
setwd("D:/我的坚果云/002项目/王老师合作项目/202207_NEW_RESULT")
getwd()

# load packages
library(MASS)
library(readr)
library(fitdistrplus)

# import the data of exp. 2
data_chisq <- read_csv("./Data_Off_distri_2022.csv")
data_chisq <- data.frame(data_chisq)

# an overview of the data
View(data_chisq)

# N levels in parental plants
Nlevel <- c(10, 20, 40, 60)

# --------------------------------------------------------------------------
# step 2: model fitting
# --------------------------------------------------------------------------
# create a null vector
out <- vector()

# loop for size distribution of each parental plants
# i: replicates of plants
# j: N levels

for (i in 1:7){
for(j in 1:4){

# select a fixed N level
Nlevel[j]

# extract data with one fixed N level
data_chisq_test  <- data_chisq[which(data_chisq$nlevels == Nlevel[j]), ]
data_temp        <- data_chisq_test[data_chisq_test$reps == paste(i), ]$propagule_mass
data_temp_mean   <- rep(mean(data_temp), 4)
data_temp_length <- rep(length(data_temp), 4)

# fitting the individual mass of clonal propagules, using "gamma", "lnorm", "weibull", and "norm" distributions
fit.ga   <- fitdist(data_chisq_test[data_chisq_test$reps == paste(i), ]$propagule_mass, "gamma")
fit.ln   <- fitdist(data_chisq_test[data_chisq_test$reps == paste(i), ]$propagule_mass, "lnorm")
fit.wb   <- fitdist(data_chisq_test[data_chisq_test$reps == paste(i), ]$propagule_mass, "weibull")
fit.norm <- fitdist(data_chisq_test[data_chisq_test$reps == paste(i), ]$propagule_mass, "norm")
out_p    <- gofstat(list(fit.ga, fit.ln, fit.wb, fit.norm), fitnames = c("gamma", "lnorm", "weibull", "norm"))

# calculate ks
out_ks     <- out_p$ks
out_kstest <- out_p$kstest

# calculate cvm
out_cvm     <- out_p$cvm
out_cvmtest <- out_p$cvmtest

# calculate chisq
out_chisq       <- out_p$chisq
out_chisqpvalue <- out_p$chisqpvalue

# calculate aic and bic
out_aic <- out_p$aic
out_bic <- out_p$bic

# calculate model estimate
fit.ga_es   <- fit.ga$estimate
fit.ln_es   <- fit.ln$estimate
fit.wb_es   <- fit.wb$estimate
fit.norm_es <- fit.norm$estimate
out_es      <- rbind(fit.ga_es, fit.ln_es, fit.wb_es, fit.norm_es)

# calculate model sd
fit.ga_sd   <- fit.ga$sd
fit.ln_sd   <- fit.ln$sd
fit.wb_sd   <- fit.wb$sd
fit.norm_sd <- fit.norm$sd
out_sd      <- rbind(fit.ga_sd, fit.ln_sd, fit.wb_sd, fit.norm_sd)

# seting model names
fit_name <- c("gamma", "lnorm", "weibull", "norm")

# combine all results
out <- rbind(out, cbind(fit_name,
                        i,
                        Nlevel[j],
                        data_temp_mean,
                        data_temp_length,
                        out_es,
                        out_sd,
                        out_chisq,
                        out_chisqpvalue,
                        out_cvm,
                        out_cvmtest,
                        out_ks,
                        out_kstest,
                        out_aic,
                        out_bic))
}
}

# --------------------------------------------------------------------------
# step 3: ouput of model results
# --------------------------------------------------------------------------
out <- as.data.frame(out,stringsAsFactors=F)

colnames(out)<-c("model_name",
                 "reps",
                 "nlevels",
                 "propagule_mean.mass",
                 "propagule_num",
                 "mean/meanlog/shape_mean",
                 "sd/sdlog/rate/scale_mean",
                 "mean/meanlog/shape_se",
                 "sd/sdlog/rate/scale_se",
                 "out_chiq",
                 "out_chiqpvalue",
                 "out_cvm",
                 "out_cvmtest",
                 "out_ks",
                 "out_kestest",
                 "out_aic",
                 "out_bic")
head(out)

# export results
write.csv(out,"./20220710_size_distribution.csv", row.names = FALSE)

#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2022-07-11 19:24
# @Description: appendix fig.3: mean mass of clonal propagules and shape and scale of weibull distribution model
#--------------------------------------------------------------------------# 
# clean cp memory
cat("\014")
rm(list = ls())
gc()

# create a work dictionary
# replace the file path with the one you used
setwd("D:/我的坚果云/002项目/王老师合作项目/202207_NEW_RESULT")
getwd()


# load packages
library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(MASS)
library(readr)

# import the data of model fitting of size distribution
data_propagule_dist <- read_csv("./20220710_size_distribution.csv")
data_propagule_dist <- data.frame(data_propagule_dist)

# an overview of the data
View(data_propagule_dist)
str(data_propagule_dist)

# --------------------------------------------------------------------------
# step 1: linear regression
# --------------------------------------------------------------------------

# use out_cvmtest to select the best-fitting weibull parameters
data_propagule_dist_01 <- data_propagule_dist %>% filter(model_name == "weibull" & out_cvmtest == "not rejected")
View(data_propagule_dist_01)

# correlation between shape value and mean mass of propagules
model_out0  <-  lm(mean.meanlog.shape_mean ~ propagule_mean.mass, data = data_propagule_dist_01)
summary(model_out0)

# correlation between scale value and mean mass of propagules
model_out1  <-  lm(sd.sdlog.rate.scale_mean ~ propagule_mean.mass, data = data_propagule_dist_01)
summary(model_out1)

# --------------------------------------------------------------------------
# step 2: regression plots of size distribution
# --------------------------------------------------------------------------
# correlation between shape value and mean mass of propagules
plot_shape <- ggplot(data_propagule_dist_01, aes(propagule_mean.mass, mean.meanlog.shape_mean)) +
		             geom_point(col = 'grey', size = 4)+
		             geom_line(aes(propagule_mean.mass, fitted(model_out0)), col = 'blue', size = 1) +
		             theme_classic() +
					 theme(legend.position = "none",
						   axis.line = element_line(size = 1, linetype = "solid"),
		                   axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
		                   axis.text = element_text(family = "serif", colour = "black", size = 14),
		                   axis.title = element_text(family = "serif", colour = "black", size = 14)) +
		             scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5)) +
		             scale_x_continuous(limits = c(30, 80), breaks = seq(30, 80, by = 10)) +
		             labs(y = "Shape", x = "Mean mass of clonal propagules (mg)")


plot_shape

# add text in plot
eq_shape <- substitute(plain(y) == gamma %.% plain(x) + delta*","~~plain(R)^2~"="~0.16*","~~italic(P)~"="~0.088, list(gamma = 0.015, delta = 0.674))
plot_shape <- plot_shape + geom_text(aes(x = 60, y = 0.3, family = "serif", size = 14, label = as.character(as.expression(eq_shape))), parse = TRUE)
plot_shape


# correlation between scale value and mean mass of propagules
plot_scale <- ggplot(data_propagule_dist_01, aes(propagule_mean.mass, sd.sdlog.rate.scale_mean)) +
                     geom_point(col = 'grey', size = 4)+
	                 geom_line(aes(propagule_mean.mass, fitted(model_out1)), col = 'blue', size = 1) +
		             theme_classic() +
					 theme(legend.position = "none",
						   axis.line = element_line(size = 1, linetype = "solid"),
		                   axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
		                   axis.text = element_text(family = "serif", colour = "black", size = 14),
		                   axis.title = element_text(family = "serif", colour = "black", size = 14)) +
					 scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
		             scale_x_continuous(limits = c(30, 80), breaks = seq(30, 80, by = 10)) +
		             labs(y = "Scale", x = "Mean mass of clonal propagules (mg)")

plot_scale
eq_scale <- substitute(plain(y) == gamma %.% plain(x) - delta*","~~plain(R)^2~"="~0.98*","~~italic(P)~"<"~0.001, list(gamma = 1.183, delta = 5.083))
plot_scale  <-  plot_scale + geom_text(aes(x = 60, y = 10, family = "serif", size = 14, label = as.character(as.expression(eq_scale))), parse = TRUE)
plot_scale

# --------------------------------------------------------------------------
# step 3: ouput plots of size distribution
# --------------------------------------------------------------------------
# combined figs
plots_distribution <- plot_shape + plot_scale + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')

# save the combined plot
ggexport(plots_distribution, filename = "./appendix.fig.3_size_distribution.png",
         width = 3000,
         height = 1500,
         pointsize = 12,
         res = 300)
