#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2021-03-09 11:24:52
# @Last Modified by:
# @Last Modified time: 2019-10-30 23:32:49
# @Description: appendix.fig.1 3D plot
#--------------------------------------------------------------------------#
# clean cp memory
cat("\014")
rm(list = ls())
gc()

# ==========================================================================
# fucntion 1: 3D plot
# ==========================================================================
# details can be found on the following webpages
# https://blog.csdn.net/huozi07/article/details/45365105
# https://www.rdocumentation.org/packages/rgl/versions/0.100.30/topics/par3d

# predict grid
predictgrid <- function(data, model, xvar, yvar, zvar, res = 16, type = NULL){
                        xrange <- range(data[[xvar]])
                        yrange <- range(data[[yvar]])
                        newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                                               y = seq(yrange[1], yrange[2], length.out = res))
                        names(newdata) <- c(xvar, yvar)
                        newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
                        newdata}

# transition matrix of x, y, z
df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL){
                  if(is.null(xvar)) xvar <- names(p)[1]
                  if(is.null(yvar)) yvar <- names(p)[2]
                  if(is.null(zvar)) zvar <- names(p)[3]
                  x <- unique(p[[xvar]])
                  y <- unique(p[[yvar]])
                  z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))
                  m <- list(x, y, z)
                  names(m)<-c(xvar, yvar, zvar)
                  m}


# two vector elements appear interleaved
interleave <- function(v1, v2) as.vector(rbind(v1, v2))

# ==========================================================================
# fucntion 2: growth model
# ==========================================================================

model = function (types, time, M0, r, N, gamma, delta, sigma, L){
                  if(missing(N)) N <- 0
                  if(missing(gamma)) gamma <- 0
                  if(missing(delta)) delta <- 0
                  if(missing(sigma)) sigma <- 0
                  if(missing(L)) L <- 0
                  K <- gamma * N + delta
                  if(missing(K)) K <- 0

                  try ({
                    if (types == "liner"){out <- M0 + r * time}
                    if (types == "expontial"){out <- M0 * exp(r * time)}
                    if (types == "power"){out <- (M0^(1-sigma) + r * time *(1-sigma))^(1/(1-sigma))}
                    if (types == "monomolecular"){out <- K - (K - M0)*exp(-r * time)}
                    if (types == "logistic.3p"){out <- M0 * K / (M0 + (K - M0) * exp (-r * time))}
                    if (types == "logistic.4p"){out <- L + M0 * (K - L)/ (M0 + (K - L - M0) * exp (-r * time))}
                    if (types == "gompertz"){out <- K * (M0 / K)^exp(- r * time)}
                  }, silent = TRUE)
                  return(out)}



# create a work dictionary
# replace the file path with the one you used
setwd("D:/我的坚果云/002项目/王老师合作项目/202207_NEW_RESULT")
getwd()

# load packages
library(dplyr)
library(readr)
library(rgl)


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

# "logistic.3p" model, starting at r = 0.1061, gamma = 163.4468, and delta = 2946.6494
model_out5 <- nls(tmass ~ model(types = "logistic.3p", time = htime, M0 = M0, r, N = nlevel, gamma, delta),
                 data = Mdata,
                 start = list(r = 0.1061, gamma = 163.4468, delta = 2946.6494),
                 algorithm = "port",
                 lower = c(0, 0, 0))
(model_out5)

# model prediction
Mdata$pred_tmass <- predict(model_out5)
mpgrid_df <- predictgrid(data = Mdata, model = model_out5, xvar = 'htime', yvar = 'nlevel', zvar = 'tmass')
mpgrid_list <- df2mat(mpgrid_df)

# --------------------------------------------------------------------------
# step 2: 3D plot
# --------------------------------------------------------------------------
# scatter plots
plot3d(Mdata$htime, Mdata$nlevel, Mdata$tmass, xlab = '', ylab = '', zlab = '', axes = FALSE, size=.5, type = 's', lit = FALSE)
spheres3d(Mdata$htime, Mdata$nlevel, Mdata$pred_tmass, gamma = 0.4, type = 's', size = 0.5, lit = FALSE)

# prediction grid
surface3d(mpgrid_list$htime, mpgrid_list$nlevel, mpgrid_list$tmass, gamma = .4, front = 'lines', back = 'lines', color = "blue")

# more sets for 3D plot
rgl.bbox(color = 'grey50', emission = 'grey50', xlen = 0, ylen = 0, zlen = 0)
rgl.material(color = 'black')
axes3d(edges = c('x--', 'y--', 'z-+'), cex = 1.5)
mtext3d('time (days)', edge = 'x--', line = 2.5, cex = 1.5)
mtext3d('N levels (mg/l)', edge = 'y--', line = 2.5, cex = 1.5)
mtext3d('Total mass (mg)', edge = 'z-+', line = 3.5, cex = 1.5)

# rotate the fig
par3d(userMatrix = rotationMatrix(0 * pi/180, 1, 0, 0));
par3d(userMatrix = rotationMatrix(-90 * pi/180, 1, -25 * pi/180, -30 * pi/180))

# export 3d plot
# rgl.snapshot("growth_model_3D.png", fmt = "png", top = TRUE)
rgl.postscript("growth_model_3D.svg", fmt = "svg")