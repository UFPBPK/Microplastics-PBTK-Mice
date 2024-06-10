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
mod <- mcode ("MPsPBTK.code", PBTK.code.eval.1) 

## RE-calibrate the model
# Input evaluation data
# A. Chen et al. (2022): Repeated daily dosing of 0.008 mg, 1 um, matrix: Blood (Unit: mg/L)
# Liu et al. (2022) : Single dose at 30 mg/kg bw, 0.8 um, matrix: Spleen, Kidney, Blood, Liver, Lung, Brain (Unit: mg/kg) 
##############################################################################################################################
## Input data set for model calibration/ oral
Obs.L.8 <-read.csv(file="Liu.8.csv")     # Liu dataset 
Obs.LB.8  = Obs.L.8 %>% select(Time, CBlood)
Obs.LGI.8 = Obs.L.8 %>% select(Time, CGI)      # GI
Obs.LLu.8 = Obs.L.8 %>% select(Time, CLung)    # Lung
Obs.LS.8  = Obs.L.8 %>% select(Time, CSpleen)  # Spleen
Obs.LK.8  = Obs.L.8 %>% select(Time, CKidney)  # Kidney
Obs.LL.8  = Obs.L.8 %>% select(Time, CLiver)   # Liver
Obs.LBr.8 = Obs.L.8 %>% select(Time, CBrain)   # Brain

Obs.A1 <-read.csv(file="A1.csv")         # Chen dataset 


# Prediction function for Liu's data
pred.eval.L <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 35    # day
  
  ## Exposure scenario L for single dose at 0.81 mg
  TDOSE.L       = 35       # Dosing frequency during exposure time
  DOSEoral.L    = 0.81    # mg, Single dose at 30 mg/kg bw, 5-week-old SPF CD-1® ICR female mice (bw 26-28 g)
  ## Oral exposure route
  ex.L <- ev(ID = 1, amt = DOSEoral.L, ii = tinterval, addl = TDOSE.L-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.L  = tgrid(0, tinterval*(TDOSE.L-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  # For Scenario L
  out.L <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
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
                           CBrain   = out.L$Brain)
  return(list("out.L"= out.L))
}


# Prediction function for Chen's data
pred.eval.A <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 16    # day
  ## Exposure scenario A for daily dose at 0.008 mg
  TDOSE.A       = 16       # Dosing frequency during exposure time
  DOSEoral.A    = 0.008    # mg, daily dose at 0.008 mg
  ## Oral exposure route
  ex.A <- ev(ID = 1, amt = DOSEoral.A, ii = tinterval, addl = TDOSE.A-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.A  = tgrid(0, tinterval*(TDOSE.A-1) + tinterval*End_time, 0.1)
  
  
  ## Get a prediction
  # For Scenario A
  out.A <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.A, tgrid = tsamp.A)
  out.A <-cbind.data.frame(Time    = out.A$time/24,     # day
                           CSpleen  = out.A$Spleen,
                           CKidney  = out.A$Kidney,
                           AUrine   = out.A$Urine,
                           CBlood   = out.A$Blood,
                           CLung    = out.A$Lung,
                           CLiver   = out.A$Liver,
                           CGI      = out.A$GI,
                           CBrain   = out.A$Brain)
  return(list("out.A" = out.A))
}


MCcost.L <-function (pars){
  out <- pred.eval.L (pars)
  cost<- modCost (model=out$out.L, obs= Obs.LB.8, x="Time")
  cost<- modCost (model=out$out.L, obs= Obs.LL.8, x="Time", cost=cost)
  cost<- modCost (model=out$out.L, obs= Obs.LLu.8, x="Time", cost=cost)
  cost<- modCost (model=out$out.L, obs= Obs.LS.8, x="Time", cost=cost)
  cost<- modCost (model=out$out.L, obs= Obs.LK.8, x="Time", cost=cost)
  cost<- modCost (model=out$out.L, obs= Obs.LGI.8, x="Time", cost=cost)
  cost<- modCost (model=out$out.L, obs= Obs.LBr.8, x="Time", cost=cost)
  return(cost)
}

MCcost.A <-function (pars){
  out <- pred.eval.A (pars)
  cost<- modCost (model=out$out.A, obs= Obs.A1, x="Time")
  return(cost)
}

## Sensitivity function (FME) 
## Check the sensitive parameters in the model based on the optimal parameters derived from model calibration
theta.final <- log(c(
  PLu    = 0.15,       # Partition coefficient     
  PL     = 0.0001,     # Adjusted
  PK     = 0.15,
  PS     = 0.15,
  PGI    = 0.15,
  PBR    = 0.15,
  PR     = 0.15,
  PALuC  = 0.001,         # Membrane-limited permeability coefficient 
  PAGIC  = 0.001,
  PALC   = 0.0001,        # Adjusted
  PAKC   = 0.001,
  PASC   = 0.00319,       # Fitted
  PABRC  = 0.000001,
  PARC   = 0.00349,         # Fitted
  KSRESrelease  = 0.003,    # Release rate constant of phagocytic cells
  KLuRESrelease = 0.005,
  KKRESrelease  = 0.01,
  KLRESrelease  = 0.0075,
  KSRESmax      = 31.03,     # Fitted  # Maximum uptake rate constant of phagocytic cells
  KLuRESmax     = 0.592,     # Fitted
  KKRESmax      = 0.964,     # Fitted
  KLRESmax      = 0.04,      # Adjusted
  KD            = 4092.26,   # Fitted 
  KGIb          = 3.42e-5,   # Fitted   # Absorption rate of GI tract (1/h)
  Kfeces        = 0.5046,    # Fitted   # Fecal clearance (L/h)
  KurineC       = 0.000564   # Fitted   # Urinary clearance (L/h/kg^0.75)
))

# Liu's
SnsPlasma.L <- sensFun(func = MCcost.L, parms = theta.final, varscale = 1)
Sen.L=summary(SnsPlasma.L)
plot(summary(SnsPlasma.L))
theta.L <- theta.final[abs(Sen.L$Mean) > 1.2*mean(abs(Sen.L$Mean))]  # Selected sensitive parameters
theta.L 
# Chen's
SnsPlasma.A <- sensFun(func = MCcost.A, parms = theta.final, varscale = 1)
Sen.A=summary(SnsPlasma.A)
plot(summary(SnsPlasma.A))
theta.A <- theta.final[abs(Sen.A$Mean) > 1.2*mean(abs(Sen.A$Mean))]
theta.A
#-----------------------------------------------------
## Liu's results
# set up necessary parameter as model input
theta.L <- log(c(PASC = 0.00319, KGIb = 3.42e-5, Kfeces = 0.5046, KurineC = 0.000564))  
theta.L <- log(c(PS = 0.0025, PLu = 0.6, PGI = 0.6, PABRC = 0.002, KGIb = 0.5, Kfeces = 0.5)) 

## PBPK model fitting 
Fit<- modFit(f=MCcost.L, p=theta.L, method ="Port", control = nls.lm.control(nprint=1))
summary(Fit)                            ## Summary of fit 
exp(Fit$par)                            ## Get the arithmetic value out of the log domain

res=MCcost(Fit$par)$residuals$res       ## Check the residual for each time points
sum(res^2)                              ## Total residuals                       

theta.final.L <- log(c(PS = 0.00225, PLu = 0.783, KGIb = 0.35, Kfeces = 0.1125,  # Fitted
                       PK = 0.52, PGI = 0.6, PABRC = 0.02, KLRESmax = 0.0114))  # Adjust

Sim.fitL = pred.eval.L(theta.final.L)$out.L         ## Simulation 

## Chen's results
## Set up sensitivity parameter as model input
theta.A <- log(c(KGIb = 3.42e-5, Kfeces = 0.5046, KurineC = 0.000564)) 
theta.A <- log(c(KGIb = 3.42e-4, Kfeces = 0.046))  
## PBPK model fitting 
Fit<- modFit(f=MCcost.A, p=theta.A, method ="Port", control = nls.lm.control(nprint=1))
summary(Fit)                            ## Summary of fit 
exp(Fit$par)                            ## Get the arithmetic value out of the log domain

theta.final.A <- log(c(KGIb = 0.000101093, Kfeces = 0.0011))   # Fitted
Sim.fitA = pred.eval.A(theta.final.A)$out.A         ## Simulation 

## Model evaluation output
df.simLS  = cbind.data.frame (Time=Sim.fitL$Time, CSpleen=Sim.fitL$CSpleen)
df.simLK  = cbind.data.frame (Time=Sim.fitL$Time, CKidney=Sim.fitL$CKidney)
df.simLL  = cbind.data.frame (Time=Sim.fitL$Time, CLiver=Sim.fitL$CLiver)
df.simLB  = cbind.data.frame (Time=Sim.fitL$Time, CBlood=Sim.fitL$CBlood)
df.simLLu = cbind.data.frame (Time=Sim.fitL$Time, CLung=Sim.fitL$CLung)
df.simLGI = cbind.data.frame (Time=Sim.fitL$Time, CGI=Sim.fitL$CGI)
df.simLBr = cbind.data.frame (Time=Sim.fitL$Time, CBrain=Sim.fitL$CBrain) 
df.simL   = cbind.data.frame (Time=Sim.fitL$Time, CSpleen=Sim.fitL$CSpleen, CKidney=Sim.fitL$CKidney,
                              CLiver=Sim.fitL$CLiver, CBlood=Sim.fitL$CBlood,
                              CLung=Sim.fitL$CLung, CGI=Sim.fitL$CGI, CBrain=Sim.fitL$CBrain) 

df.simA  = cbind.data.frame (Time=Sim.fitA$Time, CBlood=Sim.fitA$CBlood)

# Output
write.csv(df.simL, file = "1Fitting0.8_Liu.csv", row.names = FALSE)
write.csv(df.simA, file = "1Fitting Blood_Chen.csv", row.names = FALSE)

## Model calibration plot using ggplot2 
#Liu's results
plot.LB.8=
  ggplot() +
  geom_line (data = df.simLB, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.LB.8, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 

plot.LS.8=
  ggplot() +
  geom_line (data = df.simLS, aes(Time, CSpleen), col="firebrick", lwd=1)+
  geom_point(data = Obs.LS.8, aes(Time, CSpleen),size=2.5) + ylab("MPs concentration in spleen (mg/kg)") 

plot.LL.8=
  ggplot() +
  geom_line (data = df.simLL, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.LL.8, aes(Time, CLiver),size=2.5) + ylab("MPs concentration in liver (mg/kg)") 

plot.LK.8=
  ggplot() +
  geom_line (data = df.simLK, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.LK.8, aes(Time, CKidney),size=2.5) + ylab("MPs concentration in kidney (mg/kg)") 

plot.LGI.8=
  ggplot() +
  geom_line (data = df.simLGI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.LGI.8, aes(Time, CGI),size=2.5) + ylab("MPs concentration in GI (mg/kg)") 

plot.LLu.8=
  ggplot() +
  geom_line (data = df.simLLu, aes(Time, CLung), col="firebrick", lwd=1)+
  geom_point(data = Obs.LLu.8, aes(Time, CLung),size=2.5) + ylab("MPs concentration in lung (mg/kg)") 

plot.LBr.8=
  ggplot() +
  geom_line (data = df.simLBr, aes(Time, CBrain), col="firebrick", lwd=1)+
  geom_point(data = Obs.LBr.8, aes(Time, CBrain),size=2.5) + ylab("MPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("1Fitting.8_Liu.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.LB.8, plot.LS.8, plot.LK.8, plot.LL.8, plot.LLu.8, 
             plot.LGI.8, plot.LBr.8, nrow = 3)
dev.off()

#Chen's results
plot.A=
  ggplot() +
  geom_line (data = df.simA, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.A1, aes(Time, CBlood),size=2.5) + ylab("MPs concentration in blood (mg/L)") 
png("1Fitting_Chen.png", width = 6600, height = 6000, units = "px", res = 600)
plot.A
dev.off()

#----------------------------------------------------------------------------
## Model performance

# Extract predicted y values for specific x values
# Since this was a repeated dosing regimen, Time = 1 would make CGI experience the second dose. However, the observed TK data were measured at day-1 before the second dose applied.
selected_Lrows1 <- df.simL[df.simL$Time < 1, ]    # Include only the rows where the 'Time' column is less than 1.
selected_Lrows  <- tail(selected_Lrows1 , n =1)   # Select the last row
# Specify the observed x values
specific_Ax_values <- c(4, 8, 12, 16)
selected_Arows   <- df.simA[df.simA$Time %in% specific_Ax_values, ]
# Extracted y values
predicted_LBvalues  <- selected_Lrows$CBlood
predicted_LLvalues  <- selected_Lrows$CLiver
predicted_LSvalues  <- selected_Lrows$CSpleen
predicted_LKvalues  <- selected_Lrows$CKidney
predicted_LLuvalues <- selected_Lrows$CLung
predicted_LBRvalues <- selected_Lrows$CBrain
predicted_LGIvalues <- selected_Lrows$CGI

predicted_Avalues   <- selected_Arows$CBlood
# Observed y values
mean.obsLB  <- c(135.86)
mean.obsLL  <- c(69.86)
mean.obsLS  <- c(106.31)
mean.obsLK  <- c(81.56)
mean.obsLLu <- c(103.7)
mean.obsLBR <- c(27.78)
mean.obsLGI <- c(63.39)

mean.obsA   <- c(0.017015, 0.0477271, 0.08581641, 0.113208)

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultLL   <- calculate_mape(mean.obsLL, predicted_LLvalues)
mape_resultLS   <- calculate_mape(mean.obsLS, predicted_LSvalues)
mape_resultLK   <- calculate_mape(mean.obsLK, predicted_LKvalues)
mape_resultLLu  <- calculate_mape(mean.obsLLu, predicted_LLuvalues)
mape_resultLB   <- calculate_mape(mean.obsLB, predicted_LBvalues)
mape_resultLBR  <- calculate_mape(mean.obsLBR, predicted_LBRvalues)
mape_resultLGI  <- calculate_mape(mean.obsLGI, predicted_LGIvalues)

mape_resultA   <- calculate_mape(mean.obsA, predicted_Avalues)
# Print the result
cat("MAPE:", mape_resultLL, "%\n")
cat("MAPE:", mape_resultLS, "%\n")
cat("MAPE:", mape_resultLK, "%\n")
cat("MAPE:", mape_resultLB, "%\n")
cat("MAPE:", mape_resultLLu, "%\n")
cat("MAPE:", mape_resultLBR, "%\n")
cat("MAPE:", mape_resultLGI, "%\n")

cat("MAPE:", mape_resultA, "%\n")

## Overall goodness of fit for Liu's results
observed.L  <- c(mean.obsLL, mean.obsLS, mean.obsLK, mean.obsLLu, mean.obsLB, mean.obsLBR, mean.obsLGI)
predicted.L <- c(predicted_LLvalues, predicted_LSvalues, predicted_LKvalues, predicted_LLuvalues, 
                 predicted_LBvalues, predicted_LBRvalues, predicted_LGIvalues)

# Fit a linear regression model
lm_model.L <- lm(observed.L ~ predicted.L)

# Calculate R-squared
r_squared <- summary(lm_model.L)$r.squared

# Calculate adjusted R-squared
n <- length(observed.L)
p <- length(lm_model.L$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.L <- sqrt(mean((observed.L - predicted.L)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.L, "\n")


## Overall goodness of fit for Chen's results
# Fit a linear regression model
lm_model.A <- lm(mean.obsA ~ predicted_Avalues)

# Calculate R-squared
r_squared <- summary(lm_model.A)$r.squared

# Calculate adjusted R-squared
n <- length(mean.obsA)
p <- length(lm_model.A$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.A <- sqrt(mean((mean.obsA - predicted_Avalues)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.A, "\n")
