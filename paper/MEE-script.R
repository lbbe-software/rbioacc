# `rbioacc`: an R-package to analyse toxicokinetic data
# A. Ratier, V. Baudrot, M. Kaag, A. Siberchicot, C. Lopes, S. Charles#
# Submitted to Methods in Ecology and Evolution
# Updated on 2021-12-22 with R version 4.1.1
# Contact: sandrine.charles@univ-lyon1.fr

############################### Reset ##################################
rm(list = ls())

############################# Packages #################################
### install package `rbioacc`, if needed
if(is.element('rbioacc', installed.packages()[, 1]) == FALSE){
  install.packages('rbioacc', priority = "NA")}
### load package `rbioacc`
library(rbioacc)

############################# Brief worked examples #################################
# Example 1 (Figure 1) ####
### load the data set
data("Oncorhynchus_two")
### create an `rbioacc` object for data analysis
### Choose which exposure concentration to consider
data1 <- Oncorhynchus_two[Oncorhynchus_two$expw == 0.00440,]
### build the TK model (automatically chosen) according to the input data
# specify the time at which the accumulation phase ends (here 49 days)
modeldata1 <- modelData(data1, time_accumulation = 49)
### fit the TK model
m1 <- fitTK(modeldata1)
### plot the fitting results
plot(m1)
### get the bioaccumulation metrics (here BCF because exposure via water)
BCFk_all <- bioacc_metric(m1, "k") # "k" for kinetic BCF
### get a summary of the kinetic BCF
for(i in 1:ncol(BCFk_all)){
  BCFk <- quantile(BCFk_all[, i], c(0.5, 0.025, 0.975))
}
# Here the data has not reached the steady-state yet at the end
# of the accumulation phase: the steady-state BCF is not recommended.
### plot the BCF posterior probability distribution
plot(BCFk_all)
### get the TK parameter estimates
quantile_table(m1)
### get the time at which 95% of the chemical is eliminated
t95(m1)
### get the equations of the TK model built accordingly to the data
equations(m1, Oncorhynchus_two)

# Example 2 (Figures 2, 3) ####
### load the data set
data("Chironomus_benzoapyrene")
### create an `rbioacc` object for data analysis
data2 <- Chironomus_benzoapyrene
### build the TK model (automatically chosen) according to data
# specify the time at which the accumulation phase ends (here 3 days)
modeldata2 <- modelData(data2, time_accumulation = 3)
### fit the TK model
m2 <- fitTK(modeldata2, iter = 10000)
### plot the fitting result
plot(m2)
### get the bioaccumulation metrics (here BSAF because exposure via sediment)
BSAFk_all <- bioacc_metric(m2, "k") # for the kinetic BSAF
# Here it is reasonable to ask for the steady-state BSAF
BSAFss_all <- bioacc_metric(m2, "ss")
### get summaries of both BSAF
for(i in 1:ncol(BSAFk_all)){
  BSAFk <- quantile(BSAFk_all[, i], c(0.5, 0.025, 0.975))
}
for(i in 1:ncol(BSAFss_all)){
  BSAFss <- quantile(BSAFss_all[, i], c(0.5, 0.025, 0.975))
}
### display BSAF estimates
(BSAF <- t(cbind(BSAFk, BSAFss)))
### plot BSAF posterior probability distributions
plot(BSAFk_all)
plot(BSAFss_all)
### get the TK parameter estimates
quantile_table(m2)
### get the time at which 95% of the chemical is eliminated
t95(m2)
### get the equations of the TK model built accordingly to the data
equations(m2, Chironomus_benzoapyrene)
### get the PPC
ppc(m2)
ppc(m2)$labels$subtitle # to get the percentage
### get the plot of prior and posterior distributions for each parameter
plot_PriorPost(m2)
### get the correlations between parameters
corrMatrix(m2)
corrPlot(m2)
### get the PSRF value for each parameter
psrf(m2)
### get the traces of MCMC chains
mcmcTraces(m2)
### get the WAIC
waic(m2)

# Example 3 : calibration and validation (Figure 4) ####
# Calibration
### create the dataframe with data collected at 0.025 µg/mL
dfcalib <- structure(list(
  time = c(0, 0.083, 1, 2, 6, 7, 8, 12, 0, 0.083, 1, 2, 6, 7, 8, 12),
  conc = c(0, 1.2, 2.94, 3.57, 4.65, 1.52, 1.73, 1.89, 0, 4.7, 11.04,
           6.65, 10.91, 3.37, 3.13, 2.88),
  replicate = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2),
  expw = c(0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
           0.025, 0.025, 0.025, 0.025, 0.025, 0.025)),
  row.names = c(NA, -16L), class = c("tbl_df", "tbl", "data.frame"))
data3 <- dfcalib
### build the TK model according to the data
### specify the time at the end of exposure
modeldata3 <- modelData(data3, time_accumulation = 6)
### fit the TK model
m3 <- fitTK(modeldata3)
# Validation
### create the data frame with data collected at 0.1 µg/mL
dfvalid <- structure(list(time = c(0, 0.083, 1, 2, 6, 7, 8, 12, 0, 0.083, 1, 2, 6, 7, 8, 12), conc = c(0, 2.95, 10.02, 15.45, 12.82, 5.75, 5.08, 3.15, 0, 14.66, 29.54, 30.51, 28.65, 15.46, 11.77, 9.07), replicate = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2), expw = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1, 0.1, 0.1)), row.names = c(NA, -16L), class = c("tbl_df", "tbl", "data.frame"))
### load the dplyr package to create prediction data frame
library(dplyr)
### format data for predictions
data4pred <- data.frame(unique(dfvalid %>% select(time, expw)))
### perform simulations
predictm3 <- predict(m3, data4pred)
### plot simulations
p <- plot(predictm3)
p
### format data for validation
data4valid <- data.frame(unique(dfvalid %>% select(time, conc)))
### plot data against simulations
validation <- p + geom_point(data = data4valid, aes(x = time, y = conc))
validation

# Example 4 : prediction step with function predict() (Figure 5) ####
### prepare data for prediction at exposure 0.05 µg/mL
### The total duration of the experiment is 15 days
data4pred <- data.frame(time = 1:15, expw = 0.05)
plot(predict(m3, data4pred))

# Example 5 : prediction step with function predict_manual() (Figure 6) ####
### Extract posterior table for each parameter
posteriors <- rstan::extract(m3[["stanfit"]])
mcmcm3 <- data.frame(kee = posteriors$ke[, 1], kuw = posteriors$ku[, 1], sigmaConc = posteriors$sigmaCGpred[, 1])
### prepare data for prediction at exposure 0.01 µg/mL
# the total duration of the experiment is 75 days
data4pred <- data.frame(time = 1:75, expw = 0.01)
### perform predictions from median or mean parameter values (without uncertainty)
predictm1_manual <- predict_manual(param = data.frame(kee = 0.03834562, kuw = 10.56466351), data4pred, time_accumulation = 49)
### plot predictions performed manually
plot(predictm1_manual)

# Example 6: fit under time-variable exposure ####
## upload input files
data("Exposure_Sialis_lutaria") # exposure concentration profile
data("Internal_Sialis_lutaria") # internal concentrations
### fit the appropriate TK model
Exposure_Sialis_lutaria$value <- Exposure_Sialis_lutaria$Cwater
Internal_Sialis_lutaria$value <- Internal_Sialis_lutaria$Cinternal
modeldata_SL <- modelData_ode(Exposure_Sialis_lutaria, Internal_Sialis_lutaria, time_accumulation = 2.170)
fit_SL <- fitTK(modeldata_SL, iter = 100)
### get a summary of parameter estimates
quantile_table(fit_SL)