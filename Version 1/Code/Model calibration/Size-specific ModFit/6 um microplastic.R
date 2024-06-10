## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(gridExtra)   # Arrange plots in to one figure
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data


## Build mrgsolve-based PBTK Model
mod <- mcode ("PBTK.code", PBTK.code) 

############################################################################################################################## 
# Model calibration and plot for 7 data sets                                
# Keinänen et al. (2021)   : Single dose at 0.1 mg, 6 um, matrix: Spleen, Kidney, Urine, Blood, Liver, Lung, GI (Unit: mg/kg) 
##############################################################################################################################
## Input data set for model calibration/ oral
Obs.KS6 <-read.csv(file="KS6.csv")     # KS6 dataset: Spleen 
names(Obs.KS6) = c("Time", "CSpleen")

Obs.KK6 <-read.csv(file="KK6.csv")     # KK6 dataset: Kidney 
names(Obs.KK6) = c("Time", "CKidney")

Obs.KU6 <-read.csv(file="KU6.csv")     # KU6 dataset: Urine 
names(Obs.KU6) = c("Time", "AUrine")

Obs.KB6 <-read.csv(file="KB6.csv")     # KB6 dataset: Blood 
names(Obs.KB6) = c("Time", "CBlood")

Obs.KL6 <-read.csv(file="KL6.csv")     # KL6 dataset: Liver  
names(Obs.KL6) = c("Time", "CLiver")

Obs.KLu6 <-read.csv(file="KLu6.csv")   # KLu6 dataset: Lung  
names(Obs.KLu6) = c("Time", "CLung")

Obs.KGI6 <-read.csv(file="KGI6.csv")   # KGI6 dataset: GI 
names(Obs.KGI6) = c("Time", "CGI")

# Define the prediction function
## Define the prediction function (for least squares fit using levenberg-marquart algorithm)
pred <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24     # h
  End_time      = 4      # day
  ## Exposure scenario for Single dose at 0.1 mg
  TDOSE.K       = 1      # Dosing frequency during exposure time
  DOSEoral.K    = 0.1    # mg
  ## Oral exposure route
  ex.K <- ev(ID = 1, amt = DOSEoral.K, ii = tinterval, addl = TDOSE.K-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.K  = tgrid(0, tinterval*(TDOSE.K-1) + tinterval*End_time, 0.2)
  
  ## Get a prediction
  out.K <- 
    mod %>%                                               # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(Spleen, Kidney, Urine, Blood, Lung, Liver, GI, Brain, Rest, Feces)%>%                             # select model output
    update(atol = 1e-70, maxstep = 50000)%>%
    mrgsim_d(data = ex.K, tgrid = tsamp.K)
  out.K <-cbind.data.frame(Time     = out.K$time/24,     # day
                           CSpleen  = out.K$Spleen,
                           CKidney  = out.K$Kidney,
                           AUrine   = out.K$Urine,
                           CBlood   = out.K$Blood,
                           CLung    = out.K$Lung,
                           CLiver   = out.K$Liver,
                           CGI      = out.K$GI,
                           CBrain   = out.K$Brain,
                           CRest    = out.K$Rest,
                           AFeces   = out.K$Feces)
  
  return(list("out.K"=out.K))
}

## Cost function (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <- pred (pars)
  cost<- modCost  (model=out$out.K, obs= Obs.KS6, x="Time")
  cost<- modCost  (model=out$out.K, obs= Obs.KK6, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KU6, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KL6, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KLu6, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KGI6, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KB6, x="Time", cost=cost)
  return(cost)
}

## Local sensitivity analysis
## Choose the sensitive parameters in the model
## Initial parameters
theta.int <- log(c(
  PLu    = 0.15,       # Partition coefficient     
  PL     = 0.08,
  PK     = 0.15,
  PS     = 0.15,
  PGI    = 0.15,
  PBR    = 0.15,
  PR     = 0.15,
  PALuC  = 0.001,      #Membrane-limited permeability coefficient 
  PAGIC  = 0.001,
  PALC   = 0.001,
  PAKC   = 0.001,
  PASC   = 0.001,
  PABRC  = 0.000001,
  PARC   = 0.000001,
  KSRESrelease  = 0.003,      #Release rate constant of phagocytic cells
  KLuRESrelease = 0.005,
  KKRESrelease  = 0.01,
  KLRESrelease  = 0.0075,
  KSRESmax      = 10,        #Maximum uptake rate constant of phagocytic cells
  KLuRESmax     = 0.1, 
  KKRESmax      = 0.1, 
  KLRESmax      = 4, 
  ASREScap      = 200,       #Uptake capacity per tissue weight 
  ALuREScap     = 15,
  AKREScap      = 15,
  ALREScap      = 100,
  KD            = 206.15,   
  N             = 1,
  KGIb          = 6e-5,   #Absorption rate of GI tract (1/h)
  KbileC        = 0.0012,   #Biliary clearance (L/h/kg^0.75)
  Kfeces        = 0.2,     #Fecal clearance (L/h)
  KurineC       = 0.00012    #Urinary clearance (L/h/kg^0.75)
))

## PBPK model fitting 
Fit<- modFit(f=MCcost, p=theta.int, method ="Port",
             control = nls.lm.control(nprint=1))
## Model fitting based on initial settings
Sim.fitK = pred(theta.int)    

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)

Sen=summary(SnsPlasma)

plot(summary(SnsPlasma))

png("Sensitivity_plot.6.png", width = 6600, height = 6000, units = "px", res = 600)
plot(summary(SnsPlasma))
dev.off()

#-----------------------------------------------------
## Set up sensitivity or necessary parameter as model input
## Selected sensitive parameters;
theta <- theta.int[abs(Sen$Mean) > 1.5*mean(abs(Sen$Mean))]
theta
## fitted parameters
theta <- log(c(
  PLu     = 0.15,
  PS      = 0.15,
  PGI     = 0.15,
  PBR     = 0.15,
  PALC    = 0.001,
  PABRC   = 0.000001,
  KKRESrelease  = 0.01, 
  KLuRESrelease = 0.005,
  KLRESmax      = 4,
  KGIb    = 6e-5,
  Kfeces  = 0.2
))

## set up necessary parameter as model input
theta <- log(c(PASC = 0.004, PARC = 0.015,
               KSRESmax      = 30,  KKRESmax = 1.9, 
               KLuRESmax     = 0.45, KLRESmax = 0.004,
               KD            = 4500,
               KGIb          = 4e-5,    
               Kfeces        = 0.5,
               KurineC       = 0.00043
               ))

## PBPK model fitting 
# Port is faster than Nelder-Mead
Fit<- modFit(f = MCcost, p = theta, method = "Port",
             control = nls.lm.control(nprint=1))
summary(Fit)                            ## Summary of fit 
exp(Fit$par)                            ## Get the arithmetic value out of the log domain
res=MCcost(Fit$par)$residuals$res       ## Check the residual for each time points
sum(res^2)                              ## Total residuals                       

## Model calibration 
Sim.fitK = pred(Fit$par)          ## Simulation 

#Fitting
theta.final <- log(c(
  PLu    = 0.15,        # Partition coefficient     
  PL     = 0.0001,      # Adjusted
  PK     = 0.15,
  PS     = 0.15,
  PGI    = 0.15,
  PBR    = 0.15,
  PR     = 0.15,
  PALuC  = 0.001,       # Membrane-limited permeability coefficient 
  PAGIC  = 0.001,
  PALC   = 0.0001,      # Adjusted
  PAKC   = 0.001,     
  PASC   = 0.004,       # Fitted 
  PABRC  = 0.000001,
  PARC   = 0.015,           # Fitted
  KSRESrelease  = 0.001,    # Adjusted  # Release rate constant of phagocytic cells
  KLuRESrelease = 0.001,    # Adjusted
  KKRESrelease  = 0.01,
  KLRESrelease  = 0.0075,
  KSRESmax      = 29.91,    # Fitted    # Maximum uptake rate constant of phagocytic cells
  KKRESmax      = 1.87,     # Fitted
  KLuRESmax     = 0.445,    # Fitted
  KLRESmax      = 0.00402,  # Fitted
  KD            = 4468.93,  # Fitted
  N             = 1,
  KGIb          = 4.01e-5,  # Fitted  # Adjusted (decrease) #Absorption rate of GI tract (1/h)
  Kfeces        = 0.508,    # Fitted  # Fecal clearance (L/h)
  KurineC       = 0.000429  # Fitted  # Urinary clearance (L/h/kg^0.75)
))

Sim.fitK = pred(theta.final)$out.K    ## Simulation of exposure scenario K

## Model calibration results
df.simKS  = cbind.data.frame (Time=Sim.fitK$Time, CSpleen=Sim.fitK$CSpleen)
df.simKK  = cbind.data.frame (Time=Sim.fitK$Time, CKidney=Sim.fitK$CKidney)
df.simKU  = cbind.data.frame (Time=Sim.fitK$Time, AUrine=Sim.fitK$AUrine)
df.simKL  = cbind.data.frame (Time=Sim.fitK$Time, CLiver=Sim.fitK$CLiver)
df.simKB  = cbind.data.frame (Time=Sim.fitK$Time, CBlood=Sim.fitK$CBlood)
df.simKLu = cbind.data.frame (Time=Sim.fitK$Time, CLung=Sim.fitK$CLung)
df.simKGI = cbind.data.frame (Time=Sim.fitK$Time, CGI=Sim.fitK$CGI)
df.simK   = cbind.data.frame (Time=Sim.fitK$Time, 
                              CBrain=Sim.fitK$CBrain, 
                              CRest =Sim.fitK$CRest,
                              AFeces=Sim.fitK$AFeces)
# Output
write.csv(df.simKS, file = "Fitting Spleen6_Kein.csv", row.names = FALSE)
write.csv(df.simKK, file = "Fitting Kidney6_Kein.csv", row.names = FALSE)
write.csv(df.simKU, file = "Fitting Urine6_Kein.csv", row.names = FALSE)
write.csv(df.simKLu, file = "Fitting Lung6_Kein.csv", row.names = FALSE)
write.csv(df.simKL, file = "Fitting Liver6_Kein.csv", row.names = FALSE)
write.csv(df.simKB, file = "Fitting Blood6_Kein.csv", row.names = FALSE)
write.csv(df.simKGI, file = "Fitting GI6_Kein.csv", row.names = FALSE)
write.csv(df.simK, file = "Fitting Rest6_Kein.csv", row.names = FALSE)

## Model calibration plot using ggplot2 
plot.KS6=
  ggplot() +
  geom_line (data = df.simKS, aes(Time, CSpleen), col="firebrick", lwd=1)+
  geom_point(data = Obs.KS6, aes(Time, CSpleen),size=2.5) + ylab("MPs concentration in spleen (mg/kg)") 

plot.KK6=
  ggplot() +
  geom_line (data = df.simKK, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.KK6, aes(Time, CKidney),size=2.5) + ylab("MPs concentration in kidney (mg/kg)")

plot.KU6=
  ggplot() +
  geom_line (data = df.simKU, aes(Time, AUrine), col="firebrick", lwd=1)+
  geom_point(data = Obs.KU6, aes(Time, AUrine),size=2.5) + ylab("MPs amount in urine (mg)")

plot.KB6=
  ggplot() +
  geom_line (data = df.simKB, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.KB6, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 

plot.KL6=
  ggplot() +
  geom_line (data = df.simKL, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.KL6, aes(Time, CLiver),size=2.5) + ylab("MPs concentration in liver (mg/kg)") 

plot.KLu6=
  ggplot() +
  geom_line (data = df.simKLu, aes(Time, CLung), col="firebrick", lwd=1)+
  geom_point(data = Obs.KLu6, aes(Time, CLung),size=2.5) + ylab("MPs concentration in lung (mg/kg)") 

plot.KGI6=
  ggplot() +
  geom_line (data = df.simKGI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.KGI6, aes(Time, CGI),size=2.5) + ylab("MPs concentration in GI (mg/kg)") 

png("Fitting.6_Kein.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.KS6, plot.KK6, plot.KU6, plot.KGI6, plot.KB6, plot.KL6, plot.KLu6, nrow = 3)
dev.off()
#----------------------------------------------------------------------------
## Goodness of fit
# Specify the x values you are interested in
specific_Kx_values <- c(0.25, 0.5, 1, 2)

mean.obsKGI <- c(5.543, 1.329, 0.03, 0.129)
mean.obsKS  <- c(0.0028, 0.0023, 0.0106, 0.0335)
mean.obsKK  <- c(0.0007, 0.0103, 0.003, 0.0059)
mean.obsKU  <- c(7.3e-07, 2.8e-07, 6.3e-07, 1.4e-06)
mean.obsKL  <- c(0.0004, 0.0015, 0.0006, 0.0009)
mean.obsKLu <- c(0.0009, 0.0007, 0.0041, 0.0057)
mean.obsKB  <- c(0.0005, 0.0006, 0.0017, 0.0025)
# Extract y values for specific x values
selected_KSrows <- df.simKS[df.simKS$Time %in% specific_Kx_values, ]
selected_KKrows <- df.simKK[df.simKK$Time %in% specific_Kx_values, ]
selected_KUrows <- df.simKU[df.simKU$Time %in% specific_Kx_values, ]
selected_KLrows <- df.simKL[df.simKL$Time %in% specific_Kx_values, ]
selected_KLurows <- df.simKLu[df.simKLu$Time %in% specific_Kx_values, ]
selected_KBrows  <- df.simKB[df.simKB$Time %in% specific_Kx_values, ]
selected_KGIrows <- df.simKGI[df.simKGI$Time %in% specific_Kx_values, ]
# Extracted y values
predicted_KSvalues <- selected_KSrows$CSpleen
predicted_KKvalues <- selected_KKrows$CKidney
predicted_KUvalues <- selected_KUrows$AUrine
predicted_KLvalues <- selected_KLrows$CLiver
predicted_KLuvalues <- selected_KLurows$CLung
predicted_KBvalues  <- selected_KBrows$CBlood
predicted_KGIvalues <- selected_KGIrows$CGI

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultKS <- calculate_mape(mean.obsKS, predicted_KSvalues)
mape_resultKK <- calculate_mape(mean.obsKK, predicted_KKvalues)
mape_resultKU <- calculate_mape(mean.obsKU, predicted_KUvalues)
mape_resultKL <- calculate_mape(mean.obsKL, predicted_KLvalues)
mape_resultKLu <- calculate_mape(mean.obsKLu, predicted_KLuvalues)
mape_resultKB  <- calculate_mape(mean.obsKB, predicted_KBvalues)
mape_resultKGI <- calculate_mape(mean.obsKGI, predicted_KGIvalues)

# Print the result
cat("MAPE:", mape_resultKGI, "%\n")
cat("MAPE:", mape_resultKS, "%\n")
cat("MAPE:", mape_resultKK, "%\n")
cat("MAPE:", mape_resultKU, "%\n")
cat("MAPE:", mape_resultKL, "%\n")
cat("MAPE:", mape_resultKLu, "%\n")
cat("MAPE:", mape_resultKB, "%\n")

# Calculate R-squared
correlationKS <- cor(mean.obsKS, predicted_KSvalues)
correlationKK <- cor(mean.obsKK, predicted_KKvalues)
correlationKU <- cor(mean.obsKU, predicted_KUvalues)
correlationKL <- cor(mean.obsKL, predicted_KLvalues)
correlationKLu <- cor(mean.obsKLu, predicted_KLuvalues)
correlationKB  <- cor(mean.obsKB, predicted_KBvalues)
correlationKGI <- cor(mean.obsKGI, predicted_KGIvalues)
# Print the R-squared value
print(paste("R-squared:", round(correlationKS^2, 3)))
print(paste("R-squared:", round(correlationKK^2, 3)))
print(paste("R-squared:", round(correlationKU^2, 3)))
print(paste("R-squared:", round(correlationKL^2, 7)))
print(paste("R-squared:", round(correlationKLu^2, 3)))
print(paste("R-squared:", round(correlationKB^2, 3)))
print(paste("R-squared:", round(correlationKGI^2, 3)))

# Overall goodness of fit
observed <- c(mean.obsKS, mean.obsKK, mean.obsKU, mean.obsKL, mean.obsKLu, 
              mean.obsKB, mean.obsKGI)
predicted <- c(predicted_KSvalues, predicted_KKvalues, predicted_KUvalues, 
               predicted_KLvalues, predicted_KLuvalues, predicted_KBvalues, 
               predicted_KGIvalues)

# Fit a linear regression model
lm_model <- lm(observed ~ predicted)

# Calculate R-squared
r_squared <- summary(lm_model)$r.squared

# Calculate adjusted R-squared
n <- length(observed)
p <- length(lm_model$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse <- sqrt(mean((observed - predicted)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse, "\n")