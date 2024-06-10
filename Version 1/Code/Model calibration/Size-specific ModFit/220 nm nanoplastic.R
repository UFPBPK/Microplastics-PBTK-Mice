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
# K.    Keinänen et al. (2021)   : Single dose at 0.1 mg, 220 nm, matrix: Spleen, Kidney, Urine, Blood, Liver, Lung, GI (Unit: mg/kg) 
##################################################################################################################################

## Input data set for model calibration/ oral
Obs.KB220 <-read.csv(file="KB220.csv")      # KB220 dataset: matrix: Blood (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KB220) = c("Time", "CBlood")

Obs.KLu220 <-read.csv(file="KLu220.csv")    # KLu220 dataset: matrix: Lung (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KLu220) = c("Time", "CLung")

Obs.KL220 <-read.csv(file="KL220.csv")      # KL220 dataset: matrix: Liver (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KL220) = c("Time", "CLiver")

Obs.KS220 <-read.csv(file="KS220.csv")      # KS220 dataset: matrix: Spleen (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KS220) = c("Time", "CSpleen")

Obs.KK220 <-read.csv(file="KK220.csv")      # KK220 dataset: matrix: Kidney (Unit: mg/kg)/Keinänen et al. (2021) 
names(Obs.KK220) = c("Time", "CKidney")

Obs.KU220 <-read.csv(file="KU220.csv")      # KU220 dataset: matrix: Urine (Unit: mg)/Keinänen et al. (2021) 
names(Obs.KU220) = c("Time", "AUrine")

Obs.KGI220 <-read.csv(file="KGI220.csv")    # KGI220 dataset: matrix: GI tract (Unit: mg)/Keinänen et al. (2021) 
names(Obs.KGI220) = c("Time", "CGI")

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
  cost<- modCost  (model=out$out.K, obs= Obs.KS220, x="Time")
  cost<- modCost  (model=out$out.K, obs= Obs.KK220, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KU220, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KB220, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KL220, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KLu220, x="Time", cost=cost)
  cost<- modCost  (model=out$out.K, obs= Obs.KGI220, x="Time", cost=cost)
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

## Model fitting based on initial settings
Sim.fitK = pred(theta.int)$out.K         

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)

Sen=summary(SnsPlasma)

plot(summary(SnsPlasma))

png("Sensitivity_plot.220.png", width = 6600, height = 6000, units = "px", res = 600)
plot(summary(SnsPlasma))
dev.off()
#-----------------------------------------------------
## set up sensitivity or necessary parameter as model input
## Selected sensitive parameters;
theta <- theta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
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
Fit<- modFit(f=MCcost, p=theta, method ="Marq", lower = lb, upper = ub,
             control = nls.lm.control(nprint=1))

## set up necessary parameter as model input
theta <- log(c(
  PK     = 0.25,
  PASC   = 0.0024,
  PARC   = 0.004,
  KLRESmax      = 0.007,
  KKRESmax      = 0.29,
  KLuRESmax     = 0.63,
  KD            = 3200,
  KGIb          = 4.3e-5,    
  Kfeces        = 0.55,     
  KurineC       = 0.00029
))

## PBPK model fitting 
Fit<- modFit(f=MCcost, p=theta, method ="Marq", control = nls.lm.control(nprint=1))

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
  PK     = 0.251,      # Fitted  
  PS     = 0.15,
  PGI    = 0.15,
  PBR    = 0.15,
  PR     = 0.15,
  PALuC  = 0.001,      # Membrane-limited permeability coefficient 
  PAGIC  = 0.001,
  PALC   = 0.0001,     # Adjusted
  PAKC   = 0.001,       
  PASC   = 0.00249,    # Fitted   
  PABRC  = 0.000001,
  PARC   = 0.00224,    # Fitted
  KSRESrelease  = 0.003,      # Release rate constant of phagocytic cells
  KLuRESrelease = 0.005,
  KKRESrelease  = 0.01,
  KLRESrelease  = 0.0075,
  KSRESmax      = 10,        # Maximum uptake rate constant of phagocytic cells
  KLuRESmax     = 0.643,     # Fitted
  KKRESmax      = 0.291,     # Fitted   # Adjusted (increase)
  KLRESmax      = 0.00427,   # Fitted  
  KD            = 2464.7,    # Fitted  
  KGIb          = 4.22e-5,   # Fitted  # Absorption rate of GI tract (1/h)
  Kfeces        = 0.55,      # Fitted  # Fecal clearance (L/h)
  KurineC       = 0.0003     # Fitted  # Urinary clearance (L/h/kg^0.75)
))
Sim.fitK = pred(theta.final)$out.K


## Calibration results
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
write.csv(df.simKS, file = "Fitting Spleen220_Kein.csv", row.names = FALSE)
write.csv(df.simKK, file = "Fitting Kidney220_Kein.csv", row.names = FALSE)
write.csv(df.simKU, file = "Fitting Urine220_Kein.csv", row.names = FALSE)
write.csv(df.simKLu, file = "Fitting Lung220_Kein.csv", row.names = FALSE)
write.csv(df.simKL, file = "Fitting Liver220_Kein.csv", row.names = FALSE)
write.csv(df.simKB, file = "Fitting Blood220_Kein.csv", row.names = FALSE)
write.csv(df.simKGI, file = "Fitting GI220_Kein.csv", row.names = FALSE)
write.csv(df.simK, file = "Fitting All220_Kein.csv", row.names = FALSE)

## Model calibration plot using ggplot2 
plot.KS220=
  ggplot() +
  geom_line (data = df.simKS, aes(Time, CSpleen), col="firebrick", lwd=1)+
  geom_point(data = Obs.KS220, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.KK220=
  ggplot() +
  geom_line (data = df.simKK, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.KK220, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)")

plot.KU220=
  ggplot() +
  geom_line (data = df.simKU, aes(Time, AUrine), col="firebrick", lwd=1)+
  geom_point(data = Obs.KU220, aes(Time, AUrine),size=2.5) + ylab("NPs amount in urine (mg)")

plot.KB220=
  ggplot() +
  geom_line (data = df.simKB, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.KB220, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.KL220=
  ggplot() +
  geom_line (data = df.simKL, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.KL220, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.KLu220=
  ggplot() +
  geom_line (data = df.simKLu, aes(Time, CLung), col="firebrick", lwd=1)+
  geom_point(data = Obs.KLu220, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

plot.KGI220=
  ggplot() +
  geom_line (data = df.simKGI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.KGI220, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI tract (mg/kg)") 

png("Fitting.220_Kein.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.KS220, plot.KK220, plot.KU220, plot.KGI220, plot.KB220, 
             plot.KL220, plot.KLu220, nrow = 3)
dev.off()

## Statistics----------------------------------------------------------------------------
# Specify the x values you are interested in
specific_Kx_values <- c(0.25, 0.5, 1, 2)

mean.obsKS <- c(0.002, 0.0015, 0.009, 0.0103)
mean.obsKK <- c(0.0012, 0.0012, 0.0038, 0.0022)
mean.obsKU <- c(4.8e-07, 2.3e-07, 7.8e-07, 6.3e-07)
mean.obsKL <- c(0.0009, 0.0003, 0.0012, 0.0007)
mean.obsKLu <- c(0.0017, 0.0012, 0.0036, 0.0025)
mean.obsKB  <- c(0.0011, 0.0007, 0.0043, 0.0012)
mean.obsKGI <- c(4.124, 0.938, 0.055, 0.09)
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
