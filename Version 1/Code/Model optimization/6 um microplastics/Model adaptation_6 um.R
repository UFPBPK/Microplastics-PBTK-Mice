## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(pheatmap)    # Clustered heatmap
library(viridis)     # Color palette
library(RColorBrewer) # Rcolorbrewer palette

## Build mrgsolve-based PBTK Model
mod <- mcode ("MPsPBTK.code", PBTK.code.eval.6) 

## RE-calibrate the model
# Input evaluation data
# Liang et al. (2021)  : Single dose at 250 mg/kg, 5 um, matrix: Blood, GI (Unit: mg/kg) 
# Han et al. (2023) : Repeated daily dose at 50 mg/kg, 5 um, matrix: Serum, brain, liver, kidney, intestine (Unit: mg/kg)
# Zhang et al. (2024) : Repeated daily dose at 0.5 mg, 5 um, matrix: Liver, Kidney, Spleen, Lung (Unit: mg/kg)
##############################################################################################################################
## Input data set for model calibration/ oral
Obs.L5 <-read.csv(file="Liang5.csv")     # Liang dataset 
Obs.LB5 = Obs.L5 %>% select(-CGI)
Obs.LGI5 = Obs.L5 %>% select(-CBlood)
names(Obs.LB5) = c("Time", "CBlood")
names(Obs.LGI5) = c("Time", "CGI")

Obs.H5 <-read.csv(file="Han5.csv")       # Han dataset 
Obs.HB5  = Obs.H5 %>% select(Time, CBlood)   # Blood
Obs.HGI5 = Obs.H5 %>% select(Time, CGI)      # GI
Obs.HK5  = Obs.H5 %>% select(Time, CKidney)  # Kidney
Obs.HL5  = Obs.H5 %>% select(Time, CLiver)   # Liver
Obs.HBr5 = Obs.H5 %>% select(Time, CBrain)   # Brain

Obs.Z5 <-read.csv(file="Zhang5.csv")     # Zhang dataset 
Obs.ZK5  = Obs.Z5 %>% select(Time, CKidney)  # Kidney
Obs.ZL5  = Obs.Z5 %>% select(Time, CLiver)   # Liver
Obs.ZLu5 = Obs.Z5 %>% select(Time, CLung)    # Lung
Obs.ZS5  = Obs.Z5 %>% select(Time, CSpleen)  # Spleen

# Prediction function for Liang's data
pred.eval <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 1     # day
  ## Exposure scenario for Single dose at 0.1 mg
  TDOSE.L       = 1     # Dosing frequency during exposure time
  DOSEoral.L    = 5     # mg, Single dose at 250 mg/kg bw, 6-week-old C57BL/6 J mice (bw 18–20 g)
  ## Oral exposure route
  ex.L <- ev(ID = 1, amt = DOSEoral.L, ii = tinterval, addl = TDOSE.L-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.L  = tgrid(0, tinterval*(TDOSE.L-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  out.L <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    Req(Spleen, Kidney, Urine, Blood, Lung, Liver, GI, Brain, Rest, Feces)%>% 
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.L, tgrid = tsamp.L)
  out.L <-cbind.data.frame(Time    = out.L$time/24,     # day
                           CSpleen  = out.L$Spleen,
                           CKidney  = out.L$Kidney,
                           AUrine   = out.L$Urine,
                           CBlood   = out.L$Blood,
                           CLung    = out.L$Lung,
                           CLiver   = out.L$Liver,
                           CGI      = out.L$GI,
                           CBrain   = out.L$Brain,
                           CRest    = out.L$Rest,
                           AFeces   = out.L$Feces
  )
  return(list("out.L"=out.L))
}

# Prediction function for Han's data
pred.eval.H <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 90    # day
  BW            = 0.015  # kg
  ## Exposure scenario H for Single dose at 0.1 mg
  TDOSE.H       = 90     # Dosing frequency during exposure time
  DOSEoral.H    = 0.75   # mg, Daily dose at 50 mg/kg bw, 3-week-old specific pathogen-free C57BL/6N male mice (bw 15 g)
  ## Oral exposure route
  ex.H <- ev(ID = 1, amt = DOSEoral.H, ii = tinterval, addl = TDOSE.H-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.H  = tgrid(0, tinterval*(TDOSE.H-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  # For Scenario H
  out.H <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.H, tgrid = tsamp.H)
  out.H <-cbind.data.frame(Time     = out.H$time/24,    # day
                           CKidney  = out.H$Kidney,
                           CBlood   = out.H$Blood,
                           CLiver   = out.H$Liver,
                           CGI      = out.H$GI,
                           CBrain   = out.H$Brain
  )
  return(list("out.H"= out.H))
}


# Prediction function for Zhang's data
pred.eval.Z <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 56    # day
  BW            = 0.024  # kg
  ## Exposure scenario Z for Single dose at 0.1 mg
  TDOSE.Z       = 56     # Dosing frequency during exposure time
  DOSEoral.Z    = 0.5    # mg, Daily dose at 0.5 mg, 6-8-week-old C57BL/6 male mice (bw 22-26 g)
  ## Oral exposure route
  ex.Z <- ev(ID = 1, amt = DOSEoral.Z, ii = tinterval, addl = TDOSE.Z-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.Z  = tgrid(0, tinterval*(TDOSE.Z-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  # For Scenario Z
  out.Z <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.Z, tgrid = tsamp.Z)
  out.Z <-cbind.data.frame(Time     = out.Z$time/24,    # day
                           CLiver   = out.Z$Liver,
                           CSpleen  = out.Z$Spleen,
                           CKidney  = out.Z$Kidney,
                           CLung    = out.Z$Lung,
                           CGI      = out.Z$GI
  )
  return(list("out.Z"= out.Z))
}


MCcost.L<-function (pars){
  out <- pred.eval (pars)
  cost<- modCost  (model=out$out.L, obs= Obs.LB5, x="Time")
  cost<- modCost  (model=out$out.L, obs= Obs.LGI5, x="Time", cost=cost)
  return(cost)
}

MCcost.H <-function (pars){
  out <- pred.eval.H (pars)
  cost<- modCost  (model=out$out.H, obs= Obs.HB5, x="Time")
  cost<- modCost  (model=out$out.H, obs= Obs.HL5, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HK5, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HGI5, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HBr5, x="Time", cost=cost)
  return(cost)
}

MCcost.Z <-function (pars){
  out <- pred.eval.Z (pars)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZL5, x="Time")
  cost<- modCost  (model=out$out.Z, obs= Obs.ZLu5, x="Time", cost=cost)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZS5, x="Time", cost=cost)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZK5, x="Time", cost=cost)
  return(cost)
}

## Sensitivity function (FME) 
## Check the sensitive parameters in the model based on the optimal parameters derived from model calibration
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

# Local sensitivity analysis for Liang's data
SnsPlasma.L <- sensFun(func = MCcost.L, parms = theta.final, varscale = 1)
Sen.L=summary(SnsPlasma.L)
plot(summary(SnsPlasma.L))
theta.L <- theta.final[abs(Sen.L$Mean) > 1.2*mean(abs(Sen.L$Mean))]  # Selected sensitive parameters
theta.L

# Local sensitivity analysis for Han's data
SnsPlasma.H <- sensFun(func = MCcost.H, parms = theta.final, varscale = 1)
Sen.H=summary(SnsPlasma.H)
plot(summary(SnsPlasma.H))
theta.H <- theta.final[abs(Sen.H$Mean) > 1.2*mean(abs(Sen.H$Mean))]  # Selected sensitive parameters
theta.H

# Local sensitivity analysis for Zhang's data
SnsPlasma.Z <- sensFun(func = MCcost.Z, parms = theta.final, varscale = 1)
Sen.Z=summary(SnsPlasma.Z)
plot(summary(SnsPlasma.Z))
theta.Z <- theta.final[abs(Sen.Z$Mean) > 1.2*mean(abs(Sen.Z$Mean))]  # Selected sensitive parameters
theta.Z


## Re-fitted parameters based on Liang's data
# set up necessary parameter as model input
theta <- log(c(Kfeces = 0.508))  

## PBPK model fitting 
Fit<- modFit(f=MCcost.L, p=theta, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit)                            # Summary of fit 
exp(Fit$par)                            # Get the arithmetic value out of the log domain
res=MCcost(Fit$par)$residuals$res       # Check the residual for each time points
sum(res^2)                              # Total residuals                       

#Fitting Liang's data 
theta.final <- log(c(KGIb   = 0.0026,   # Adjust
                     Kfeces = 0.1019))  # Fitted   
Sim.fitL = pred.eval(theta.final)$out.L # Simulation of exposure scenario L


## Re-fitted parameters based on Han's data
# Set up sensitivity or necessary parameter as model input
theta.H <- log(c(PK = 0.15, PASC  = 0.004, KKRESrelease  = 0.01,
                 KGIb = 4.01e-5, KurineC  = 0.000328))
theta.H <- log(c(KKRESmax = 0.07, KLRESmax = 0.06, KGIb = 1e-4, Kfeces = 0.5))
## PBPK model fitting 
Fit.H<- modFit(f=MCcost.H, p=theta.H, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit.H)                            ## Summary of fit 
exp(Fit.H$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.H(Fit.H$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

#Fitting Han's data
theta.final.H <- log(c(
  BW    = 0.015,     # Body weight changed based on Choi et al. (2021)
  PABRC = 0.0001, PBR = 0.2,    # Adjust 
  KKRESmax = 0.189, KLRESmax = 0.0623, Kfeces = 0.3293,   # Fitted
  KGIb  = 1.07e-4))  # Adjust (decrease)
    
Sim.fitH = pred.eval.H(theta.final.H)$out.H   ## Simulation for scenario H


## Re-fitted parameters based on Zhang's data
# Set up sensitivity or necessary parameter as model input
theta.Z <- log(c(PASC = 0.002712, KGIb = 4.052e-5, KurineC = 0.000309)) 
theta.Z <- log(c(KSRESmax = 0.015, KGIb = 0.008, KurineC = 0.00005)) 
## PBPK model fitting 
Fit.Z<- modFit(f=MCcost.Z, p=theta.Z, method ="Marq", control = nls.lm.control(nprint=1))

summary(Fit.Z)                            ## Summary of fit 
exp(Fit.Z$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.Z(Fit.Z$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

#Fitting Zhang's data
theta.final.Z <- log(c(
  BW = 0.024,    # Body weight changed based on Zhang et al. (2024)
  PK = 2.13, PLu = 0.0001,      # Adjust
  KSRESmax = 0.0167, KGIb = 0.0077655, KurineC = 0.00005439)) # Fitted


Sim.fitZ = pred.eval.Z(theta.final.Z)$out.Z   ## Simulation for scenario Z


## Model evaluation output
df.simLB  = cbind.data.frame (Time=Sim.fitL$Time, CBlood=Sim.fitL$CBlood)
df.simLGI = cbind.data.frame (Time=Sim.fitL$Time, CGI=Sim.fitL$CGI)
df.simL   = cbind.data.frame (Time=Sim.fitL$Time, CBlood=Sim.fitL$CBlood, CGI=Sim.fitL$CGI) 

df.simHK  = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney)
df.simHL  = cbind.data.frame (Time=Sim.fitH$Time, CLiver=Sim.fitH$CLiver)
df.simHB  = cbind.data.frame (Time=Sim.fitH$Time, CBlood=Sim.fitH$CBlood)
df.simHGI = cbind.data.frame (Time=Sim.fitH$Time, CGI=Sim.fitH$CGI)
df.simHBr = cbind.data.frame (Time=Sim.fitH$Time, CBrain=Sim.fitH$CBrain) 
df.simH   = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney,
                              CLiver=Sim.fitH$CLiver, CBlood=Sim.fitH$CBlood,
                              CGI=Sim.fitH$CGI, CBrain=Sim.fitH$CBrain)

df.simZK  = cbind.data.frame (Time=Sim.fitZ$Time, CKidney=Sim.fitZ$CKidney)
df.simZL  = cbind.data.frame (Time=Sim.fitZ$Time, CLiver=Sim.fitZ$CLiver)
df.simZS  = cbind.data.frame (Time=Sim.fitZ$Time, CSpleen=Sim.fitZ$CSpleen)
df.simZLu = cbind.data.frame (Time=Sim.fitZ$Time, CLung=Sim.fitZ$CLung)
df.simZ   = cbind.data.frame (Time=Sim.fitZ$Time, CKidney=Sim.fitZ$CKidney,
                              CLiver=Sim.fitZ$CLiver, CSpleen=Sim.fitZ$CSpleen,
                              CLung=Sim.fitZ$CLung)
# Output
write.csv(df.simL, file = "Fitting eval.5_Liang.csv", row.names = FALSE)
write.csv(df.simH, file = "Fitting eval.5_Han.csv", row.names = FALSE)
write.csv(df.simZ, file = "Fitting eval.5_Zhang.csv", row.names = FALSE)

## Model evaluation plot using ggplot2 
#Liang's results
plot.LB5=
  ggplot() +
  geom_line (data = df.simLB, aes(Time, CBlood), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LB5, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 

plot.LGI5=
  ggplot() +
  geom_line (data = df.simLGI, aes(Time, CGI), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LGI5, aes(Time, CGI),size=2.5) + ylab("MPs concentration in GI (mg/kg)") 

# Arrange plots into one figure and save as png file
png("6Fitting.5_Liang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.LB5, plot.LGI5, nrow = 1)
dev.off()


#Han's results
plot.HB5=
  ggplot() +
  geom_line (data = df.simHB, aes(Time, CBlood), col="blue", lwd=1)+
  geom_point(data = Obs.HB5, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 

plot.HL5=
  ggplot() +
  geom_line (data = df.simHL, aes(Time, CLiver), col="blue", lwd=1)+
  geom_point(data = Obs.HL5, aes(Time, CLiver),size=2.5) + ylab("MPs concentration in liver (mg/kg)") 

plot.HK5=
  ggplot() +
  geom_line (data = df.simHK, aes(Time, CKidney), col="blue", lwd=1)+
  geom_point(data = Obs.HK5, aes(Time, CKidney),size=2.5) + ylab("MPs concentration in kidney (mg/kg)") 

plot.HGI5=
  ggplot() +
  geom_line (data = df.simHGI, aes(Time, CGI), col="blue", lwd=1)+
  geom_point(data = Obs.HGI5, aes(Time, CGI),size=2.5) + ylab("MPs concentration in GI (mg/kg)") 

plot.HBr5=
  ggplot() +
  geom_line (data = df.simHBr, aes(Time, CBrain), col="blue", lwd=1)+
  geom_point(data = Obs.HBr5, aes(Time, CBrain),size=2.5) + ylab("MPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("6Fitting.5_Han.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.HB5, plot.HK5, plot.HL5, plot.HGI5, plot.HBr5, nrow = 3)
dev.off()


#Zhang's results
plot.ZS5=
  ggplot() +
  geom_line (data = df.simZS, aes(Time, CSpleen), col="purple", lwd=1)+
  geom_point(data = Obs.ZS5, aes(Time, CSpleen),size=2.5) + ylab("MPs concentration in spleen (mg/kg)") 

plot.ZL5=
  ggplot() +
  geom_line (data = df.simZL, aes(Time, CLiver), col="purple", lwd=1)+
  geom_point(data = Obs.ZL5, aes(Time, CLiver),size=2.5) + ylab("MPs concentration in liver (mg/kg)") 

plot.ZK5=
  ggplot() +
  geom_line (data = df.simZK, aes(Time, CKidney), col="purple", lwd=1)+
  geom_point(data = Obs.ZK5, aes(Time, CKidney),size=2.5) + ylab("MPs concentration in kidney (mg/kg)") 

plot.ZLu5=
  ggplot() +
  geom_line (data = df.simZLu, aes(Time, CLung), col="purple", lwd=1)+
  geom_point(data = Obs.ZLu5, aes(Time, CLung),size=2.5) + ylab("MPs concentration in lung (mg/kg)") 


# Arrange plots into one figure and save as png file
png("6Fitting.5_Zhang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.ZS5, plot.ZK5, plot.ZL5, plot.ZLu5, nrow = 2)
dev.off()

#----------------------------------------------------------------------------
## Model performance of Liang's results
specific_Lx_values <- c(1)       # Specify the observed x values
# Extract y values for specific x values
selected_LBrows  <- df.simLB[df.simLB$Time %in% specific_Lx_values, ]
selected_LGIrows <- df.simLGI[df.simLGI$Time %in% specific_Lx_values, ]
# Extracted y values
predicted_LBvalues  <- selected_LBrows$CBlood
predicted_LGIvalues <- selected_LGIrows$CGI

mean.obsLB  <- c(18.0587)
mean.obsLGI <- c(515.778)

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultLB  <- calculate_mape(mean.obsLB, predicted_LBvalues)
mape_resultLGI <- calculate_mape(mean.obsLGI, predicted_LGIvalues)

# Print the result
cat("MAPE:", mape_resultLGI, "%\n")
cat("MAPE:", mape_resultLB, "%\n")

# Overall goodness of fit
observed <- c(mean.obsLB, mean.obsLGI)
predicted <- c(predicted_LBvalues, predicted_LGIvalues)

lm_model <- lm(observed ~ predicted)      # Fit a linear regression model
r_squared <- summary(lm_model)$r.squared  # Calculate R-squared

# Calculate adjusted R-squared
n <- length(observed)
p <- length(lm_model$coefficients) - 1    # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

rmse <- sqrt(mean((observed - predicted)^2))   # Calculate Root Mean Squared Error (RMSE)

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse, "\n")


## Model performance of Han's results
# Specify the observed x values
specific_Hx_values <- c(90)

# Extract predicted y values for specific x values
selected_HBrows  <- df.simHB[df.simHB$Time %in% specific_Hx_values, ]
selected_HLrows  <- df.simHL[df.simHL$Time %in% specific_Hx_values, ]
selected_HKrows  <- df.simHK[df.simHK$Time %in% specific_Hx_values, ]
selected_HBRrows <- df.simHBr[df.simHBr$Time %in% specific_Hx_values, ]
selected_HGIrows <- df.simHGI[df.simHGI$Time %in% specific_Hx_values, ]
# Extracted y values
predicted_HBvalues  <- selected_HBrows$CBlood
predicted_HLvalues  <- selected_HLrows$CLiver
predicted_HKvalues  <- selected_HKrows$CKidney
predicted_HBRvalues <- selected_HBRrows$CBrain
predicted_HGIvalues <- selected_HGIrows$CGI
# Observed y values
mean.obsHB  <- c(0.27238)
mean.obsHL  <- c(1.588)
mean.obsHK  <- c(0.6833)
mean.obsHBR <- c(0.0613)
mean.obsHGI <- c(0.3767)

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultHL   <- calculate_mape(mean.obsHL, predicted_HLvalues)
mape_resultHK   <- calculate_mape(mean.obsHK, predicted_HKvalues)
mape_resultHB   <- calculate_mape(mean.obsHB, predicted_HBvalues)
mape_resultHBR  <- calculate_mape(mean.obsHBR, predicted_HBRvalues)
mape_resultHGI  <- calculate_mape(mean.obsHGI, predicted_HGIvalues)

# Print the result
cat("MAPE:", mape_resultHL, "%\n")
cat("MAPE:", mape_resultHK, "%\n")
cat("MAPE:", mape_resultHB, "%\n")
cat("MAPE:", mape_resultHBR, "%\n")
cat("MAPE:", mape_resultHGI, "%\n")

# Overall goodness of fit
observed.H  <- c(mean.obsHL, mean.obsHK, mean.obsHB, mean.obsHBR, mean.obsHGI)
predicted.H <- c(predicted_HLvalues, predicted_HKvalues,predicted_HBvalues, 
                 predicted_HBRvalues, predicted_HGIvalues)

# Fit a linear regression model
lm_model.H <- lm(observed.H ~ predicted.H)

# Calculate R-squared
r_squared <- summary(lm_model.H)$r.squared

# Calculate adjusted R-squared
n <- length(observed.H)
p <- length(lm_model.H$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.H <- sqrt(mean((observed.H - predicted.H)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.H, "\n")


## Model performance of Zhang's results
# Specify the observed x values
specific_Zx_values <- c(56)

# Extract predicted y values for specific x values
selected_ZLrows  <- df.simZL[df.simZL$Time %in% specific_Zx_values, ]
selected_ZSrows  <- df.simZS[df.simZS$Time %in% specific_Zx_values, ]
selected_ZKrows  <- df.simZK[df.simZK$Time %in% specific_Zx_values, ]
selected_ZLurows <- df.simZLu[df.simZLu$Time %in% specific_Zx_values, ]

# Extracted y values
predicted_ZLvalues  <- selected_ZLrows$CLiver
predicted_ZSvalues  <- selected_ZSrows$CSpleen
predicted_ZKvalues  <- selected_ZKrows$CKidney
predicted_ZLuvalues <- selected_ZLurows$CLung

# Observed y values
mean.obsZL  <- c(39.9113)
mean.obsZS  <- c(49.6674)
mean.obsZK  <- c(95.7447)
mean.obsZLu <- c(22.5745)

# Calculate MAPE 
mape_resultZL   <- calculate_mape(mean.obsZL, predicted_ZLvalues)
mape_resultZS   <- calculate_mape(mean.obsZS, predicted_ZSvalues)
mape_resultZK   <- calculate_mape(mean.obsZK, predicted_ZKvalues)
mape_resultZLu  <- calculate_mape(mean.obsZLu, predicted_ZLuvalues)

# Print the result
cat("MAPE:", mape_resultZL, "%\n")
cat("MAPE:", mape_resultZS, "%\n")
cat("MAPE:", mape_resultZK, "%\n")
cat("MAPE:", mape_resultZLu, "%\n")

# Overall goodness of fit
observed.Z  <- c(mean.obsZL, mean.obsZS, mean.obsZK, mean.obsZLu)
predicted.Z <- c(predicted_ZLvalues, predicted_ZSvalues, predicted_ZKvalues, predicted_ZLuvalues)

# Fit a linear regression model
lm_model.Z <- lm(observed.Z ~ predicted.Z)

# Calculate R-squared
r_squared <- summary(lm_model.Z)$r.squared

# Calculate adjusted R-squared
n <- length(observed.Z)
p <- length(lm_model.Z$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.Z <- sqrt(mean((observed.Z - predicted.Z)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.Z, "\n")