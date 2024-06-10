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
################################################################################################################################## 
# Model calibration and plot for 7 data sets                                
# K.    Keinänen et al. (2021)   : Single dose at 0.1 mg, 20 nm, matrix: Spleen, Kidney, Urine, Blood, Liver, Lung, GI (Unit: mg/kg) 
##################################################################################################################################

## Input data set for model calibration/ oral
Obs.KB20 <-read.csv(file="KB20.csv")      # KB20 dataset: matrix: Blood (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KB20) = c("Time", "CBlood")

Obs.KLu20 <-read.csv(file="KLu20.csv")    # KLu20 dataset: matrix: Lung (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KLu20) = c("Time", "CLung")

Obs.KL20 <-read.csv(file="KL20.csv")      # KL20 dataset: matrix: Liver (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KL20) = c("Time", "CLiver")

Obs.KS20 <-read.csv(file="KS20.csv")      # KS20 dataset: matrix: Spleen (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KS20) = c("Time", "CSpleen")

Obs.KK20 <-read.csv(file="KK20.csv")      # KK20 dataset: matrix: Kidney (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KK20) = c("Time", "CKidney")

Obs.KU20 <-read.csv(file="KU20.csv")      # KU20 dataset: matrix: Urine (Unit: mg)/Keinänen et al. (2021) 
names(Obs.KU20) = c("Time", "AUrine")

Obs.KGI20 <-read.csv(file="KGI20.csv")    # KGI20 dataset: matrix: GI tract (Unit: mg)/Keinänen et al. (2021) 
names(Obs.KGI20) = c("Time", "CGI")

# Define the prediction function
## Define the prediction function (for least squares fit using levenberg-marquart algorithm)
pred <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24     # h
  End_time      = 4     # day
  
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
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.K, tgrid = tsamp.K)
  out.K <-cbind.data.frame(Time     = out.K$time/24,     # day
                           CBlood   = out.K$Blood,
                           CSpleen  = out.K$Spleen,
                           CKidney  = out.K$Kidney,
                           CLung    = out.K$Lung,
                           CLiver   = out.K$Liver,
                           CGI      = out.K$GI,
                           CBrain   = out.K$Brain,
                           CRest    = out.K$Rest,
                           AUrine   = out.K$Urine,
                           AFeces   = out.K$Feces)
  
  return(list("out.K"=out.K))
}

## Cost function (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  out <- pred (pars)
  cost<- modCost  (model=out$out.K, obs= Obs.KS20, x="Time")
  cost<- modCost  (model=out$out.K, obs= Obs.KK20, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KU20, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KB20, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KL20, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KLu20, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KGI20, x="Time", cost=cost)
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
  KGIb          = 6e-5,     #Absorption rate of GI tract (1/h)
  KbileC        = 0.0012,   #Biliary clearance (L/h/kg^0.75)
  Kfeces        = 0.2,      #Fecal clearance (L/h)
  KurineC       = 0.00012   #Urinary clearance (L/h/kg^0.75)
))


## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)

Sen=summary(SnsPlasma)

plot(summary(SnsPlasma))

png("Sensitivity_plot.20.png", width = 6600, height = 6000, units = "px", res = 600)
plot(summary(SnsPlasma))
dev.off()
#-----------------------------------------------------
## set up sensitivity or necessary parameter as model input
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

# Set boundary for each parameter
lb <- log(c(
  PLu     = 1e-6,
  PS      = 1e-6,
  PGI     = 1e-6,
  PBR     = 1e-6,
  PALC    = 1e-6,
  PABRC   = 1e-7,
  KKRESrelease  = 1e-6, 
  KLuRESrelease = 1e-6,
  KLRESmax      = 1e-6,
  KGIb    = 1e-6,
  Kfeces  = 1e-6
))

ub <- log(c(
  PLu     = 6,
  PS      = 6,
  PGI     = 6,
  PBR     = 6,
  PALC    = 6,
  PABRC   = 6,
  KKRESrelease  = 6, 
  KLuRESrelease = 6,
  KLRESmax      = 50,
  KGIb    = 6,
  Kfeces  = 6
))

## PBPK model fitting 
## least squares fit using levenberg-marquart (method "Marq") algorithm
Fit<- modFit(f=MCcost, p=theta, method ="Marq", 
             control = nls.lm.control(nprint=1))

## set up necessary parameter as model input
theta <- log(c(
  PK            = 0.37, 
  PARC          = 0.01,
  KLuRESrelease = 0.5, 
  KSRESmax      = 1,
  KD            = 500,
  KGIb          = 4.4e-5,   
  Kfeces        = 0.5
))
## PBPK model fitting 
Fit<- modFit(f=MCcost, p=theta, method ="Port",
             control = nls.lm.control(nprint=1))
summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain
res=MCcost(Fit$par)$residuals$res            ## Check the residual for each time points
sum(res^2)                                   ## Total residuals                       

## Model calibration 
Sim.fitK = pred(Fit$par)$out.K

#Fitting
theta.final <- log(c(
  PLu    = 0.15,       # Partition coefficient     
  PL     = 0.0001,     # Adjusted
  PK     = 0.374,      # Fitted  
  PS     = 0.15,
  PGI    = 0.15,
  PBR    = 0.15,
  PR     = 0.15,
  PALuC  = 0.001,      # Membrane-limited permeability coefficient 
  PAGIC  = 0.001,
  PALC   = 0.0001,     # Adjusted
  PAKC   = 0.001,       
  PASC   = 0.001,      
  PABRC  = 0.000001,
  PARC   = 0.00695,         # Fitted
  KSRESrelease  = 0.003,    # Release rate constant of phagocytic cells
  KLuRESrelease = 0.53,     # Fitted
  KKRESrelease  = 0.01,       
  KLRESrelease  = 0.0075,    
  KSRESmax      = 0.905,    # Fitted # Maximum uptake rate constant of phagocytic cells
  KLuRESmax     = 0.1,       
  KKRESmax      = 0.1,       
  KLRESmax      = 0.001,    # Adjusted  
  KD            = 558.04,   # Fitted
  KGIb          = 4.23e-5,  # Fitted  # Absorption rate of GI tract (1/h)
  Kfeces        = 0.548,    # Fitted  # Fecal clearance (L/h)
  KurineC       = 0.00012   # Urinary clearance (L/h/kg^0.75)
))

Sim.fitK = pred(theta.final)$out.K

## Model calibration plot using ggplot2 
## Calibration results of exposure scenario K
df.simKS   = cbind.data.frame (Time=Sim.fitK$Time, CSpleen=Sim.fitK$CSpleen)
df.simKK   = cbind.data.frame (Time=Sim.fitK$Time, CKidney=Sim.fitK$CKidney)
df.simKU   = cbind.data.frame (Time=Sim.fitK$Time, AUrine=Sim.fitK$AUrine)
df.simKL   = cbind.data.frame (Time=Sim.fitK$Time, CLiver=Sim.fitK$CLiver)
df.simKLu  = cbind.data.frame (Time=Sim.fitK$Time, CLung =Sim.fitK$CLung)
df.simKB   = cbind.data.frame (Time=Sim.fitK$Time, CBlood=Sim.fitK$CBlood)
df.simKGI  = cbind.data.frame (Time=Sim.fitK$Time, CGI=Sim.fitK$CGI)
df.simK    = cbind.data.frame (Time=Sim.fitK$Time, 
                               CBrain=Sim.fitK$CBrain, 
                               CRest =Sim.fitK$CRest,
                               AFeces=Sim.fitK$AFeces)
write.csv(df.simKS, file = "Fitting Spleen20_Kein.csv", row.names = FALSE)
write.csv(df.simKK, file = "Fitting Kidney20_Kein.csv", row.names = FALSE)
write.csv(df.simKU, file = "Fitting Urine20_Kein.csv", row.names = FALSE)
write.csv(df.simKLu, file = "Fitting Lung20_Kein.csv", row.names = FALSE)
write.csv(df.simKL, file = "Fitting Liver20_Kein.csv", row.names = FALSE)
write.csv(df.simKB, file = "Fitting Blood20_Kein.csv", row.names = FALSE)
write.csv(df.simKGI, file = "Fitting GI20_Kein.csv", row.names = FALSE)
write.csv(df.simK, file = "Fitting All20_Kein.csv", row.names = FALSE)

## Plot
plot.KS20=
  ggplot() +
  geom_line (data = df.simKS, aes(Time, CSpleen), col="firebrick", lwd=1)+
  geom_point(data = Obs.KS20, aes(Time, CSpleen),size=2.5) + ylab("MPs concentration in spleen (mg/kg)") 

plot.KK20=
  ggplot() +
  geom_line (data = df.simKK, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.KK20, aes(Time, CKidney),size=2.5) + ylab("MPs concentration in kidney (mg/kg)")

plot.KU20=
  ggplot() +
  geom_line (data = df.simKU, aes(Time, AUrine), col="firebrick", lwd=1)+
  geom_point(data = Obs.KU20, aes(Time, AUrine),size=2.5) + ylab("MPs amount in urine (mg)")

plot.KB20=
  ggplot() +
  geom_line (data = df.simKB, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.KB20, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 

plot.KL20=
  ggplot() +
  geom_line (data = df.simKL, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.KL20, aes(Time, CLiver),size=2.5) + ylab("MPs concentration in liver (mg/kg)") 

plot.KLu20=
  ggplot() +
  geom_line (data = df.simKLu, aes(Time, CLung), col="firebrick", lwd=1)+
  geom_point(data = Obs.KLu20, aes(Time, CLung),size=2.5) + ylab("MPs concentration in lung (mg/kg)") 

plot.KGI20=
  ggplot() +
  geom_line (data = df.simKGI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.KGI20, aes(Time, CGI),size=2.5) + ylab("MPs concentration in GI tract (mg/kg)") 

png("Fitting.20_Kein.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.KS20, plot.KK20, plot.KU20, plot.KGI20, plot.KB20, plot.KL20, plot.KLu20, nrow = 3)
dev.off()

## Statistics----------------------------------------------------------------------------
# Specify the x values you are interested in
specific_Kx_values <- c(0.25, 0.5, 1, 2)

mean.obsKS <- c(0.00205, 0.00093, 0.00152, 0.00969)
mean.obsKK <- c(0.00171, 0.0018, 0.00073, 0.00278)
mean.obsKU <- c(1.1e-06, 6.1e-07, 7e-08, 4.7e-07)
mean.obsKL <- c(0.00145, 0.00068, 0.00021, 0.00082)
mean.obsKLu <- c(0.00328, 0.00077, 0.00057, 0.00289)
mean.obsKB <- c(0.00343, 0.00242, 0.00036, 0.00239)
mean.obsKGI <- c(4.405, 0.801, 0.062, 0.027)
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
mape_resultKB <- calculate_mape(mean.obsKB, predicted_KBvalues)
mape_resultKGI <- calculate_mape(mean.obsKGI, predicted_KGIvalues)
# Print the result
cat("MAPE:", mape_resultKS, "%\n")
cat("MAPE:", mape_resultKK, "%\n")
cat("MAPE:", mape_resultKU, "%\n")
cat("MAPE:", mape_resultKL, "%\n")
cat("MAPE:", mape_resultKLu, "%\n")
cat("MAPE:", mape_resultKB, "%\n")
cat("MAPE:", mape_resultKGI, "%\n")

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
print(paste("R-squared:", round(correlationKL^2, 3)))
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