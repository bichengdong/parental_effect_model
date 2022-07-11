#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2021-01-16 10:14
# @Description: simulation of parental and offspring generations
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
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(readr)
library(Rmisc)
library(stringr)
library(tidyr)
library(dplyr)

# ==========================================================================
#
#
# growth simulation of parental generation
#
#
# ==========================================================================
# --------------------------------------------------------------------------
# step 1: preparation
# --------------------------------------------------------------------------
# initial mass (mg) of parental plants
M0 <- 22.08

# innate growth rate for three-parameter logistic model in Appendix Table 2
r <- 0.103

# parameters gamma and delta for three-parameter logistic model in Appendix Table 2
gamma <- 219.642
delta <- 1965.396

# growth time of parental plants
time <- 75

# N levels in parental generation (mg/l)
Nlevels <- seq(10, 60, 10)

# replicates for simulation of growth in parental generation
replicate <- 5

# create null vector
out_pa_gen <- vector()

# first loop for each replicate of plant in parental generation
# --------------------------------------------------------------------------
for(j in 1:replicate){

# print replicate #j
print(j)

# set random seed
    seed <- j
    set.seed(seed)

# second loop for each level of N in parental generation
# --------------------------------------------------------------------------
for (i in 1:length(Nlevels)){

# N level in each loop
    N <- Nlevels[i]

# environmental capacity K, determined by external N level
    K <- gamma * N + delta

# --------------------------------------------------------------------------
# step 2: growth simulation for parental generation
# --------------------------------------------------------------------------
# total mass of parental generation
    parent_tmass <- M0 * K / (M0 + (K - M0) * exp (-r * time))

# number of clonal propagules produced by parental plants
    alpha1 <- 0.253
    beta1 <- 0.684
    propagule_num <- alpha1*(parent_tmass)^beta1
    propagule_num <- floor(propagule_num)

# mean mass of clonal propagules produced by parental plants
    alpha2 <- 2.727
    beta2 <- 0.354
    propagule_mean.mass <- alpha2*(parent_tmass)^beta2

# size distribution of clonal propagules produced by parental plants
# weibull distributon
    alpha3 <- 0.015
    beta3  <- 0.674
    shape  <- alpha3 * propagule_mean.mass + beta3

    alpha4 <- 1.183
    beta4  <- -5.083
    scale  <- alpha4 * propagule_mean.mass + beta4

# individual size (mass) of clonal propagules produced by parental plants
    pesudo_propagule.indiv.mass <- rweibull(propagule_num, shape = shape, scale = scale)
    propagule.indiv.mass <- propagule_mean.mass / mean(pesudo_propagule.indiv.mass) * pesudo_propagule.indiv.mass

# output results
    row_seed                   <- rep(seed, propagule_num)
    row_parent_i               <- rep(i, propagule_num)
    row_parent_M0              <- rep(M0, propagule_num)
    row_parent_r               <- rep(r, propagule_num)
    row_parent_N               <- rep(N, propagule_num)
    row_parent_K               <- rep(K, propagule_num)
    row_parent_time            <- rep(time, propagule_num)
    row_parent_tmass           <- rep(parent_tmass, propagule_num)
    row_propagule_ID           <- seq(propagule_num)
    row_propagule_num          <- rep(propagule_num, propagule_num)
    row_propagule_mean.mass    <- rep(propagule_mean.mass, propagule_num)
    row_propagule.indiv.mass   <- propagule.indiv.mass

out_pa_gen <- rbind(out_pa_gen, cbind(row_seed,
                                      row_parent_i,
                                      row_parent_M0,
                                      row_parent_r,
                                      row_parent_N,
                                      row_parent_K,
                                      row_parent_time,
                                      row_parent_tmass,
                                      row_propagule_ID,
                                      row_propagule_num,
                                      row_propagule_mean.mass,
                                      row_propagule.indiv.mass
                                      ))}}

# --------------------------------------------------------------------------
# step 3: output of growth simulation for parental generation
# --------------------------------------------------------------------------
out_pa_gen <- data.frame(out_pa_gen)

# change the colnames
colnames(out_pa_gen)<-c("seed",
                        "parent_i",
                        "parent_M0",
                        "parent_r",
                        "parent_N",
                        "parent_K",
                        "parent_time",
                        "parent_tmass",
                        "propagule_ID",
                        "propagule_num",
                        "propagule_mean.mass",
                        "propagule.indiv.mass"
                        )

# output the growth simulation of paretal generation
write_excel_csv(out_pa_gen,"./20220710_simulation_pa_gen.csv")

# ==========================================================================
#
#
# survival simulation of clonal propagules
#
#
# ==========================================================================
# refs
# https://www.jianshu.com/p/0f2b7a571743
# https://www.cda.cn/discuss/post/details/5c8f4a7a5c9f1749e4f8b61b

# --------------------------------------------------------------------------
# step 1: preparation
# --------------------------------------------------------------------------
# create null vector for simulation of survival rate
out_survival <- vector()

# four fixed survival rate
survival_rate <- c(25, 50, 75, 100)

# --------------------------------------------------------------------------
# step 2: survival simulation
# --------------------------------------------------------------------------
# loop for each level of survival rate in clonal propagules
# --------------------------------------------------------------------------
for (i in 1: length(survival_rate)){

print(i)

# set random seed
  set.seed(i)

# choose one level of survival rate
  s_rate <- survival_rate[i]

# add a new column called s_rate
  out_pa_gen_temp_01 <- out_pa_gen %>% mutate(survival_rate = s_rate)

# rows were grouped by parent_i and randomly selected based on s_rate
  out_pa_gen_temp_02 <- out_pa_gen_temp_01 %>% group_by(seed, parent_i) %>% sample_n(floor(n() * s_rate/100))
  out_pa_gen_temp_02 <- out_pa_gen_temp_02 %>% arrange(seed, parent_i, propagule_ID)

# 将随机选择的数据，临时
  out_pa_gen_temp_02 <- data.frame(out_pa_gen_temp_02)
  out_survival <- rbind(out_survival, out_pa_gen_temp_02)
}

# --------------------------------------------------------------------------
# step 3: output of surivial simulation for clonal propagules
# --------------------------------------------------------------------------
out_pa_gen_survival <- out_survival

# output the survival simulation
write_excel_csv(out_pa_gen_survival, "./20220709_simulation_pa_gen_surivial.csv")

# ==========================================================================
#
#
# growth simulation of offspring generation
#
#
# ==========================================================================
# --------------------------------------------------------------------------
# step 1: preparation
# --------------------------------------------------------------------------
out_off_gen_survival <- out_pa_gen_survival

# innate growth rate for three-parameter Logistic model in Appendix Table 2
off_r <- 0.103

# parameters gamma and delta for three-parameter Logistic model in Appendix Table 2
off_gamma <- 219.642
off_delta <- 1965.396

# growth time for offspring plants
off_growth_periods <- c(30, 45, 60, 75, 150, 300)

# N levels in offspring generation (mg/l)
off_Nlevels <- seq(10, 60, 10)

# --------------------------------------------------------------------------
# step 2: growth simulation of clonal offspring
# --------------------------------------------------------------------------
# the first loop for each growth time in offspring generation
# --------------------------------------------------------------------------
for (off_t in 1: length(off_growth_periods)){

print(off_t)

# select one growth time of offspring generation for each loop
  off_tm <- off_growth_periods[off_t]

# the second loop for each level of N in offspring generation
# --------------------------------------------------------------------------
  for (off_n in 1:length(off_Nlevels)){

# N level in offspring generation for each loop
    off_N <- off_Nlevels[off_n]

# environmental capacity K, determined by external N level
    off_K <- off_gamma * off_N + off_delta

# create null vector for offspring generation
    out_off_gen <- vector()

# the third loop for each survival rate in offspring generation
# --------------------------------------------------------------------------
    for(off_i in 1:nrow(out_off_gen_survival)){

# initial mass of clonal propagule
      off_M0 <- out_off_gen_survival$propagule.indiv.mass[off_i]

# growth simulation of clonal offspring
      off_tmass <- off_M0 * off_K / (off_M0 + (off_K - off_M0) * exp (- off_r * off_tm))

# combine rows from out_off_gen and off_tmass
      out_off_gen <- rbind(out_off_gen, off_tmass)
      colnames(out_off_gen) <- paste("off_tmass_day_", off_tm, "_nlevel_", off_N, sep = "")

    }

    out_off_gen_survival <- cbind(out_off_gen_survival, out_off_gen)
  }
}

# --------------------------------------------------------------------------
# step 3: output of surivial simulation for clonal propagules
# --------------------------------------------------------------------------
# output the survival simulation
write_excel_csv(out_off_gen_survival, "./20220711_simulation_pa_gen_survival_off_gen.csv")

#--------------------------------------------------------------------------#
# -*- coding: utf-8 -*-
# @Author: Bicheng Dong
# @Date: 2020-11-10 20:27:25
# @Last Modified by:
# @Last Modified time: 2021-01-16 10:14
# @Description: analyses of parental and offspring generations
#--------------------------------------------------------------------------#
# ==========================================================================
#
#
# analyses of clonal propagules
#
#
# ==========================================================================
# load the data of simulation of parental and offspring generations
off_gen <- read_csv("./20220711_simulation_pa_gen_survival_off_gen.csv")

# a overview of data
View(off_gen)
glimpse(off_gen)

# > colnames(off_gen)
#  [1] "seed"                        "parent_i"
#  [3] "parent_M0"                   "parent_r"
#  [5] "parent_N"                    "parent_K"
#  [7] "parent_time"                 "parent_tmass"
#  [9] "propagule_ID"                "propagule_num"
# [11] "propagule_mean.mass"         "propagule.indiv.mass"
# [13] "survival_rate"               "off_tmass_day_30_nlevel_10"
# [15] "off_tmass_day_30_nlevel_20"  "off_tmass_day_30_nlevel_30"
# [17] "off_tmass_day_30_nlevel_40"  "off_tmass_day_30_nlevel_50"
# [19] "off_tmass_day_30_nlevel_60"  "off_tmass_day_45_nlevel_10"
# [21] "off_tmass_day_45_nlevel_20"  "off_tmass_day_45_nlevel_30"
# [23] "off_tmass_day_45_nlevel_40"  "off_tmass_day_45_nlevel_50"
# [25] "off_tmass_day_45_nlevel_60"  "off_tmass_day_60_nlevel_10"
# [27] "off_tmass_day_60_nlevel_20"  "off_tmass_day_60_nlevel_30"
# [29] "off_tmass_day_60_nlevel_40"  "off_tmass_day_60_nlevel_50"
# [31] "off_tmass_day_60_nlevel_60"  "off_tmass_day_75_nlevel_10"
# [33] "off_tmass_day_75_nlevel_20"  "off_tmass_day_75_nlevel_30"
# [35] "off_tmass_day_75_nlevel_40"  "off_tmass_day_75_nlevel_50"
# [37] "off_tmass_day_75_nlevel_60"  "off_tmass_day_150_nlevel_10"
# [39] "off_tmass_day_150_nlevel_20" "off_tmass_day_150_nlevel_30"
# [41] "off_tmass_day_150_nlevel_40" "off_tmass_day_150_nlevel_50"
# [43] "off_tmass_day_150_nlevel_60" "off_tmass_day_300_nlevel_10"
# [45] "off_tmass_day_300_nlevel_20" "off_tmass_day_300_nlevel_30"
# [47] "off_tmass_day_300_nlevel_40" "off_tmass_day_300_nlevel_50"
# [49] "off_tmass_day_300_nlevel_60"

# --------------------------------------------------------------------------
# step 2: summarize the number, mean mass and total mass of clonal propagules
# produced by parental plants
# --------------------------------------------------------------------------
parent_gen <- off_gen %>% group_by(seed, parent_N, survival_rate) %>% dplyr::summarise(seed = unique(seed),
                                                                                       parent_i = unique(parent_i),
                                                                                       parent_M0 = unique(parent_M0),
                                                                                       parent_r = unique(parent_r),
                                                                                       parent_N = unique(parent_N),
                                                                                       parent_K = unique(parent_K),
                                                                                       parent_time = unique(parent_time),
                                                                                       propagule_tmass = sum(propagule.indiv.mass),
                                                                                       propagule_num = n(),
                                                                                       propagule_mean.mass = propagule_tmass/propagule_num
                                                                                        )

# a overview of data
View(parent_gen)
glimpse(parent_gen)

# --------------------------------------------------------------------------
# step 3: table 2 - analysis of clonal propagules
# --------------------------------------------------------------------------
# anova
propagule_tmass.aov <- aov(propagule_tmass ~ parent_N * survival_rate, data = parent_gen)
summary(propagule_tmass.aov)

propagule_num.aov <- aov(propagule_num ~ parent_N * survival_rate, data = parent_gen)
summary(propagule_num.aov)

propagule_mean.mass.aov <- aov(propagule_mean.mass ~ parent_N * survival_rate, data = parent_gen)
summary(propagule_mean.mass.aov)

# > propagule_tmass.aov <- aov(propagule_tmass ~ parent_N * survival_rate, data = parent_gen)
# > summary(propagule_tmass.aov)
#                         Df    Sum Sq   Mean Sq F value Pr(>F)
# parent_N                 1 342180599 342180599    7098 <2e-16 ***
# survival_rate            1 542900528 542900528   11261 <2e-16 ***
# parent_N:survival_rate   1  61585872  61585872    1277 <2e-16 ***
# Residuals              116   5592279     48209
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# >
# >
# > propagule_num.aov <- aov(propagule_num ~ parent_N * survival_rate, data = parent_gen)
# > summary(propagule_num.aov)
#                         Df Sum Sq Mean Sq F value Pr(>F)
# parent_N                 1  35754   35754    8549 <2e-16 ***
# survival_rate            1 126005  126005   30129 <2e-16 ***
# parent_N:survival_rate   1   6871    6871    1643 <2e-16 ***
# Residuals              116    485       4
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# >
# >
# > propagule_mean.mass.aov <- aov(propagule_mean.mass ~ parent_N * survival_rate, data = parent_gen)
# > summary(propagule_mean.mass.aov)
#                         Df Sum Sq Mean Sq F value Pr(>F)
# parent_N                 1   9143    9143 347.960 <2e-16 ***
# survival_rate            1     82      82   3.136 0.0792 .
# parent_N:survival_rate   1     29      29   1.117 0.2928
# Residuals              116   3048      26
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# summarizes data
propagule_tmass_plot.data <- summarySE(parent_gen, measurevar = "propagule_tmass", groupvars = c("parent_N","survival_rate"))
propagule_tmass_plot.data

propagule_num_plot.data <- summarySE(parent_gen, measurevar = "propagule_num", groupvars = c("parent_N","survival_rate"))
propagule_num_plot.data

propagule_mean.mass_plot.data <- summarySE(parent_gen, measurevar = "propagule_mean.mass", groupvars = c("parent_N","survival_rate"))
propagule_mean.mass_plot.data

# --------------------------------------------------------------------------
# step 4.1: fig. 2 - plots - number, mean mass, and total mass of survived clonal propagules
# --------------------------------------------------------------------------
pd <- position_dodge(0.1) # move them .05 to the left and right

plot_propagule_tmass <- ggplot(propagule_tmass_plot.data, aes(x = factor(parent_N), y = propagule_tmass, colour = factor(survival_rate))) +
                        geom_errorbar(aes(ymin = propagule_tmass - se, ymax = propagule_tmass + se), width = .1, position = pd) +
                        geom_line(position = pd, size = 4) +
                        geom_point(position = pd, size = 4) +
                        theme_classic() +
                        theme(legend.title = element_text(family = "serif", colour = "black", size = 14),
                              legend.text = element_text(family = "serif", colour = "black", size = 14),
                              legend.background = element_rect(fill = NA),
                              legend.position = c(0.8, 0.88),
                              axis.line = element_line(size = 1, linetype = "solid"),
                              axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                              axis.text = element_text(family = "serif", colour = "black", size = 14),
                              axis.title = element_text(family = "serif", colour = "black", size = 14)) +
                              scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, by = 5000)) +
                              scale_x_discrete(breaks = seq(10, 60, 10)) +
                              labs(y = "Summed mass of survived propagules (mg)",
                                   x = "N levels in parental generation (mg/L)",
                                   color = "Survival rate (%)")

(plot_propagule_tmass)

# --------------------------------------------------------------------------
# step 4.2: plot number of clonal propagules
# --------------------------------------------------------------------------

plot_propagule_num <- ggplot(propagule_num_plot.data, aes(x = factor(parent_N), y = propagule_num, colour = factor(survival_rate))) +
                            geom_errorbar(aes(ymin = propagule_num - se, ymax = propagule_num + se), width = .1, position = pd) +
                            geom_line(position = pd, size = 4) +
                            geom_point(position = pd, size = 4) +
                            theme_classic() +
                            theme(
                                  # legend.title = element_text(family = "serif", colour = "black", size = 14),
                                  # legend.text = element_text(family = "serif", colour = "black", size = 14),
                                  # legend.position = c(0.4, 0.8),
                                  legend.position = "none",
                                  axis.line = element_line(size = 1, linetype = "solid"),
                                  axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                  axis.text = element_text(family = "serif", colour = "black", size = 14),
                                  axis.title = element_text(family = "serif", colour = "black", size = 14)) +
                                  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 50)) +
                                  scale_x_discrete(breaks = seq(10, 60, 10)) +
                                  labs(y = "Number of survived propagules",
                                       x = "N levels in parental generation (mg/L)",
                                       color = "Survival rate (%)")

(plot_propagule_num)

# --------------------------------------------------------------------------
# step 4.3: plot mean mass of clonal propagules
# --------------------------------------------------------------------------
plot_propagule_mean.mass <- ggplot(propagule_mean.mass_plot.data, aes(x = factor(parent_N), y = propagule_mean.mass, colour = factor(survival_rate))) +
                             geom_errorbar(aes(ymin = propagule_mean.mass - se, ymax = propagule_mean.mass + se), width = .1, position = pd) +
                             geom_line(position = pd, size = 4) +
                             geom_point(position = pd, size = 4) +
                             theme_classic() +
                             theme(
                                   # legend.title = element_text(family = "serif", colour = "black", size = 14),
                                   # legend.text = element_text(family = "serif", colour = "black", size = 14),
                                   # legend.position = c(0.8, 0.2),
                                  legend.position = "none",
                                  axis.line = element_line(size = 1, linetype = "solid"),
                                  axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                  axis.text = element_text(family = "serif", colour = "black", size = 14),
                                  axis.title = element_text(family = "serif", colour = "black", size = 14)) +
                                  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
                                  scale_x_discrete(breaks = seq(10, 60, 10)) +
                                  labs(y = "Mean mass of survived propagules (mg)",
                                       x = "N levels in parental generation (mg/L)",
                                       color = "Survival rate (%)")

(plot_propagule_mean.mass)

# --------------------------------------------------------------------------
# step 5: ouput plots of clonal propagules
# --------------------------------------------------------------------------
# combined figs
fig_pa_gen <- plot_propagule_tmass +
              plot_propagule_num +
              plot_propagule_mean.mass +
              plot_layout(ncol = 3) +
              plot_annotation(tag_levels = 'A')

# export the combined plot
ggexport(fig_pa_gen, filename = "./fig.2_propagules.png",
         width = 4500,
         height = 1500,
         pointsize = 12,
         res = 300)

# ==========================================================================
#
#
# analyses of performance of clonal offspring
#
#
# ==========================================================================
# --------------------------------------------------------------------------
# step 1: preparation
# --------------------------------------------------------------------------
# change wide data to long data
off_gen_01 <- off_gen %>% gather(key = 'off_treatment', value = 'off_final_mass', -c(seed:survival_rate))
off_gen_01

# N levels and growth time in offspring generation
off_Nlevels <- seq(10, 60, 10)
off_growth_periods <- c(30, 45, 60, 75, 150, 300)

# create two columns named off_N and off_time
off_gen_02 <- off_gen_01 %>% mutate(off_N = NA, off_time = NA)

# add the info of off_N and off_time in the two newly added columns
for (n in 1:length(off_Nlevels)){
                                 nlevel <- off_Nlevels[n]
                                 search_grid <- paste("nlevel_", nlevel, sep = "")
                                 off_gen_02$off_N[str_detect(off_gen_01$off_treatment, search_grid)] <- nlevel
                                 }

for (t in 1:length(off_growth_periods)){
                                        off_tm <- off_growth_periods[t]
                                        search_grid <- paste("day_", off_tm, sep = "")
                                        off_gen_02$off_time[str_detect(off_gen_02$off_treatment, search_grid)] <- off_tm
                                       }

# a overview of data
View(off_gen_02)
glimpse(off_gen_02)

# calculate summed mass and mean mass of offspring generation
off_gen_03 <- off_gen_02 %>% group_by(seed, parent_N, off_N, survival_rate, off_time) %>% dplyr::summarise(parent_i = unique(parent_i),
                                                                                                           parent_M0 = unique(parent_M0),
                                                                                                           parent_r = unique(parent_r),
                                                                                                           parent_N = unique(parent_N),
                                                                                                           parent_K = unique(parent_K),
                                                                                                           parent_time = unique(parent_time),
                                                                                                           parent_tmass = unique(parent_tmass),
                                                                                                           total_propagule_num = unique(propagule_num),
                                                                                                           survived_propagule_num = n(),
                                                                                                           off_N = unique(off_N),
                                                                                                           mean_off_mass = mean(off_final_mass)/1000,
                                                                                                           total_off_mass = sum(off_final_mass)/1000
                                                                                                           )

View(off_gen_03)
glimpse(off_gen_03)

# --------------------------------------------------------------------------
# step 2: anova
# --------------------------------------------------------------------------
total_off_mass.aov <- aov(total_off_mass ~ parent_N * off_N * survival_rate * off_time, data = off_gen_03)
summary(total_off_mass.aov)

# > summary(total_off_mass.aov)
#                                         Df    Sum Sq   Mean Sq F value   Pr(>F)
# parent_N                                 1  64790497  64790497 1605.66  < 2e-16 ***
# off_N                                    1 111137949 111137949 2754.26  < 2e-16 ***
# survival_rate                            1 193714906 193714906 4800.72  < 2e-16 ***
# off_time                                 1 113924742 113924742 2823.33  < 2e-16 ***
# parent_N:off_N                           1   7660420   7660420  189.84  < 2e-16 ***
# parent_N:survival_rate                   1  12311052  12311052  305.10  < 2e-16 ***
# off_N:survival_rate                      1  22423824  22423824  555.72  < 2e-16 ***
# parent_N:off_time                        1   4770902   4770902  118.23  < 2e-16 ***
# off_N:off_time                           1  26030624  26030624  645.10  < 2e-16 ***
# survival_rate:off_time                   1  23287334  23287334  577.12  < 2e-16 ***
# parent_N:off_N:survival_rate             1   1452403   1452403   35.99 2.14e-09 ***
# parent_N:off_N:off_time                  1   1193869   1193869   29.59 5.64e-08 ***
# parent_N:survival_rate:off_time          1    941843    941843   23.34 1.40e-06 ***
# off_N:survival_rate:off_time             1   5310143   5310143  131.60  < 2e-16 ***
# parent_N:off_N:survival_rate:off_time    1    234045    234045    5.80   0.0161 *
# Residuals                             4304 173671817     40351
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


mean_off_mass.aov <- aov(mean_off_mass ~ parent_N * off_N * survival_rate * off_time, data = off_gen_03)
summary(mean_off_mass.aov)

# > summary(mean_off_mass.aov)
#                                         Df Sum Sq Mean Sq  F value   Pr(>F)
# parent_N                                 1    104     104   17.061 3.69e-05 ***
# off_N                                    1  21226   21226 3481.078  < 2e-16 ***
# survival_rate                            1      1       1    0.099   0.7533
# off_time                                 1  22383   22383 3670.915  < 2e-16 ***
# parent_N:off_N                           1     15      15    2.539   0.1111
# parent_N:survival_rate                   1      0       0    0.046   0.8307
# off_N:survival_rate                      1      0       0    0.015   0.9035
# parent_N:off_time                        1     38      38    6.295   0.0121 *
# off_N:off_time                           1   5090    5090  834.850  < 2e-16 ***
# survival_rate:off_time                   1      0       0    0.037   0.8474
# parent_N:off_N:survival_rate             1      0       0    0.007   0.9333
# parent_N:off_N:off_time                  1      5       5    0.765   0.3819
# parent_N:survival_rate:off_time          1      0       0    0.018   0.8937
# off_N:survival_rate:off_time             1      0       0    0.005   0.9459
# parent_N:off_N:survival_rate:off_time    1      0       0    0.002   0.9623
# Residuals                             4304  26243       6
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# --------------------------------------------------------------------------
# step 3: fig.s1 - s4 - summed mass of offspring generation
# fig. s4 == fig. 3
# --------------------------------------------------------------------------
total_off_mass_plot.data <- summarySE(off_gen_03, measurevar = "total_off_mass", groupvars = c("parent_N", "off_N", "survival_rate", "off_time"))
total_off_mass_plot.data

mean_off_mass_plot.data <- summarySE(off_gen_03, measurevar = "mean_off_mass", groupvars = c("parent_N", "off_N", "survival_rate", "off_time"))
mean_off_mass_plot.data

off_growth_periods <- c(30, 45, 60, 75, 150, 300)
survival_rate <- c(25, 50, 75, 100)

# create null list
plot.list <- vector('list', length(off_growth_periods) * length(survival_rate))
plot.list

# loop for summed mass of offspring generation
for (s in 1:length(survival_rate)){
  print(s)

  off_srate <- survival_rate[s]
  total_off_mass_data_srate <- total_off_mass_plot.data %>% filter(survival_rate == off_srate)

for (t in 1:length(off_growth_periods)){

    off_tm <- off_growth_periods[t]
    total_off_mass_data_srate_tm <- total_off_mass_data_srate %>% filter(off_time == off_tm)

# main plot
    i <- (s - 1) * length(off_growth_periods) + t
    plot.list[[i]] <-  ggplot(total_off_mass_data_srate_tm, aes(x = factor(off_N), y = factor(parent_N), fill = total_off_mass)) +
                             # geom_point(stat = "identity", size = 4) +
                             # geom_line(size = 1) +
                             geom_tile() +
                             scale_fill_viridis_c() +
                             theme_tufte() +
                             theme(legend.title = element_text(family = "serif", colour = "black", size = 14),
                                   legend.text = element_text(family = "serif", colour = "black", size = 14),
                                   # legend.position = c(60, 60),
                                   axis.line = element_line(size = 1, linetype = "solid"),
                                   axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                   axis.text = element_text(family = "serif", colour = "black", size = 14),
                                   axis.title = element_text(family = "serif", colour = "black", size = 14)) +
                             labs(y = "Parental N levels (mg/L)",
                                  x = "Offspring N levels (mg/L)",
                                  fill = paste("Summed mass", "\n", "of offspring (g)", sep = "")) +
                             scale_x_discrete(breaks = seq(10, 60, 10)) +
                             scale_y_discrete(breaks = seq(10, 60, 10)) +
                             # scale_fill_gradientn(colors = c("blue","yellow"), limits = c(0,2500))+
                             ggtitle(paste("Survival rate: ", off_srate, " % ; Time: ", off_tm, " days", sep = ""))

}
}

# survival rate = 25%
fig_off_gen_tmass_srate25 <- plot.list[[1]] +
                             plot.list[[2]] +
                             plot.list[[3]] +
                             plot.list[[4]] +
                             plot.list[[5]] +
                             plot.list[[6]] +
                             plot_layout(ncol = 2) +
                             plot_annotation(tag_levels = 'A')

fig_off_gen_tmass_srate25

# survival rate = 50%
fig_off_gen_tmass_srate50 <- plot.list[[7]] +
                             plot.list[[8]] +
                             plot.list[[9]] +
                             plot.list[[10]] +
                             plot.list[[11]] +
                             plot.list[[12]] +
                             plot_layout(ncol = 2) +
                             plot_annotation(tag_levels = 'A')

fig_off_gen_tmass_srate50

# survival rate = 75%
fig_off_gen_tmass_srate75 <- plot.list[[13]] +
                             plot.list[[14]] +
                             plot.list[[15]] +
                             plot.list[[16]] +
                             plot.list[[17]] +
                             plot.list[[18]] +
                             plot_layout(ncol = 2) +
                             plot_annotation(tag_levels = 'A')

fig_off_gen_tmass_srate75

# survival rate = 100%
fig_off_gen_tmass_srate100 <- plot.list[[19]] +
                              plot.list[[20]] +
                              plot.list[[21]] +
                              plot.list[[22]] +
                              plot.list[[23]] +
                              plot.list[[24]] +
                              plot_layout(ncol = 2) +
                              plot_annotation(tag_levels = 'A')

fig_off_gen_tmass_srate100

# save the combined plot of summed mass of offspring generation
# fig.s1-s4
ggexport(fig_off_gen_tmass_srate25, filename = "./fig.s1_off_gen_tmass_srate25.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_tmass_srate50, filename = "./fig.s2_off_gen_tmass_srate50.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_tmass_srate75, filename = "./fig.s3_off_gen_tmass_srate75.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_tmass_srate100, filename = "./fig.s4_off_gen_tmass_srate100.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

# --------------------------------------------------------------------------
# step 4: fig.s5 - s8 - mean mass of offspring generation
# fig.s8 == fig.4
# --------------------------------------------------------------------------
# create null list
plot.list02 <- vector('list', length(off_growth_periods) * length(survival_rate))
plot.list02

# loop for mean mass of offspring generation
for (s in 1:length(survival_rate)){

  off_srate <- survival_rate[s]
  mean_off_mass_data_srate <- mean_off_mass_plot.data %>% filter(survival_rate == off_srate)

for (t in 1:length(off_growth_periods)){

    off_tm <- off_growth_periods[t]
    mean_off_mass_data_srate_tm <- mean_off_mass_data_srate %>% filter(off_time == off_tm)

# main plot
    i <- (s - 1) * length(off_growth_periods) + t

    plot.list02[[i]] <- ggplot(mean_off_mass_data_srate_tm, aes(x = factor(off_N), y = factor(parent_N), fill = mean_off_mass)) +
                             # geom_point(stat = "identity", size = 4) +
                             # geom_line(size = 1) +
                             geom_tile() +
                             scale_fill_viridis_c() +
                             theme_tufte() +
                             theme(legend.title = element_text(family = "serif", colour = "black", size = 14),
                                   legend.text = element_text(family = "serif", colour = "black", size = 14),
                                   # legend.position = c(60, 60),
                                   axis.line = element_line(size = 1, linetype = "solid"),
                                   axis.ticks = element_line(colour = "black", linetype = "solid", size = 1),
                                   axis.text = element_text(family = "serif", colour = "black", size = 14),
                                   axis.title = element_text(family = "serif", colour = "black", size = 14)) +
                             labs(y = "Parental N levels (mg/L)",
                                  x = "Offspring N levels (mg/L)",
                                  fill = paste("Mean mass", "\n", "of offspring (g)", sep = "")) +
                             scale_x_discrete(breaks = seq(10, 60, 10)) +
                             scale_y_discrete(breaks = seq(10, 60, 10)) +
                             # scale_fill_gradientn(colors = c("blue","yellow"), limits = c(0,2500))+
                             ggtitle(paste("Survival rate: ", off_srate, " % ; Time: ", off_tm, " days", sep = ""))

}
}

# survival rate = 25%
fig_off_gen_mean.mass_srate25 <- plot.list02[[1]] +
                                 plot.list02[[2]] +
                                 plot.list02[[3]] +
                                 plot.list02[[4]] +
                                 plot.list02[[5]] +
                                 plot.list02[[6]] +
                                 plot_layout(ncol = 2) +
                                 plot_annotation(tag_levels = 'A')

fig_off_gen_mean.mass_srate25

# survival rate = 50%
fig_off_gen_mean.mass_srate50 <- plot.list02[[7]] +
                                 plot.list02[[8]] +
                                 plot.list02[[9]] +
                                 plot.list02[[10]] +
                                 plot.list02[[11]] +
                                 plot.list02[[12]] +
                                 plot_layout(ncol = 2) +
                                 plot_annotation(tag_levels = 'A')

fig_off_gen_mean.mass_srate50

# survival rate = 75%
fig_off_gen_mean.mass_srate75 <- plot.list02[[13]] +
                                 plot.list02[[14]] +
                                 plot.list02[[15]] +
                                 plot.list02[[16]] +
                                 plot.list02[[17]] +
                                 plot.list02[[18]] +
                                 plot_layout(ncol = 2) +
                                 plot_annotation(tag_levels = 'A')

fig_off_gen_mean.mass_srate75

# survival rate = 100%
fig_off_gen_mean.mass_srate100 <- plot.list02[[19]] +
                                  plot.list02[[20]] +
                                  plot.list02[[21]] +
                                  plot.list02[[22]] +
                                  plot.list02[[23]] +
                                  plot.list02[[24]] +
                                  plot_layout(ncol = 2) +
                                  plot_annotation(tag_levels = 'A')

fig_off_gen_mean.mass_srate100

# save the combined plot of mean mass of offspring generation
# fig.s5-s8
ggexport(fig_off_gen_mean.mass_srate25, filename = "./fig.s5_off_gen_mean.mass_srate25.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_mean.mass_srate50, filename = "./fig.s6_off_gen_mean.mass_srate50.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_mean.mass_srate75, filename = "./fig.s7_off_gen_mean.mass_srate75.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)

ggexport(fig_off_gen_mean.mass_srate100, filename = "./fig.s8_off_gen_mean.mass_srate100.png",
         width = 3000,
         height = 3300,
         pointsize = 12,
         res = 300)