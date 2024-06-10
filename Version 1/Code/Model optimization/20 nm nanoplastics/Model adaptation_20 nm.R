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
mod <- mcode ("NPsPBTK.code", PBTK.code.eval.20) 

## Model evaluation   
## RE-calibrate the model
# Input evaluation data
# Liang et al. (2021) : Single dose at 250 mg/kg bw, 50 nm, matrix: Spleen, Kidney, Blood, Liver, Lung, GI (Unit: mg/kg) 
# Han et al. (2023)   : Repeated daily dose at 50 mg/kg, 50 nm, matrix: Serum, brain, liver, kidney, intestine (Unit: mg/kg)
# Wang et al. (2023)  : Single dose at 2.5, 25, 250 mg/kg, 500 nm, matrix: Blood (Unit: mg/L) 
# Tao et al. (2024)   : Repeated daily dose at 224 ug, 80 nm, matrix: Liver (Unit: mg/kg)
# Ma et al. (2024)    : Repeated daily dose at 1 mg, 25, 50, 100 nm, matrix: Blood, lung, liver, spleen, kidney, brain, GI (Unit: mg/kg)
##############################################################################################################################
## Input data set for model calibration/ oral
Obs.L50 <-read.csv(file="Liang50.csv")     # Liang dataset 
Obs.LB50  = Obs.L50 %>% select(Time, CBlood)   # Blood
Obs.LGI50 = Obs.L50 %>% select(Time, CGI)      # GI
Obs.LLu50 = Obs.L50 %>% select(Time, CLung)    # Lung
Obs.LS50  = Obs.L50 %>% select(Time, CSpleen)  # Spleen
Obs.LK50  = Obs.L50 %>% select(Time, CKidney)  # Kidney
Obs.LL50  = Obs.L50 %>% select(Time, CLiver)   # Liver
Obs.LBr50 = Obs.L50 %>% select(Time, CBrain)   # Brain

Obs.H50 <-read.csv(file="Han50.csv")     # Han dataset 
Obs.HB50  = Obs.H50 %>% select(Time, CBlood)   # Blood
Obs.HGI50 = Obs.H50 %>% select(Time, CGI)      # GI
Obs.HK50  = Obs.H50 %>% select(Time, CKidney)  # Kidney
Obs.HL50  = Obs.H50 %>% select(Time, CLiver)   # Liver
Obs.HBr50 = Obs.H50 %>% select(Time, CBrain)   # Brain

Obs.W50L <-read.csv(file="Wang50L.csv")         # Wang dataset-Low dose
Obs.W50M <-read.csv(file="Wang50M.csv")         # Wang dataset-Medium dose
Obs.W50H <-read.csv(file="Wang50H.csv")         # Wang dataset-High dose

Obs.T80 <-read.csv(file="Tao80.csv")     # Tao dataset 

Obs.M25 <-read.csv(file="Ma25.csv")      # Ma dataset: 25 nm
Obs.MK25  = Obs.M25 %>% select(Time, CKidney)  # Kidney
Obs.ML25  = Obs.M25 %>% select(Time, CLiver)   # Liver
Obs.MLu25 = Obs.M25 %>% select(Time, CLung)    # Lung
Obs.MS25  = Obs.M25 %>% select(Time, CSpleen)  # Spleen
Obs.MB25  = Obs.M25 %>% select(Time, CBlood)   # Blood
Obs.MGI25 = Obs.M25 %>% select(Time, CGI)      # GI
Obs.MBr25 = Obs.M25 %>% select(Time, CBrain)   # Brain

Obs.M50 <-read.csv(file="Ma50.csv")      # Ma dataset: 50 nm
Obs.MK50  = Obs.M50 %>% select(Time, CKidney)  # Kidney
Obs.ML50  = Obs.M50 %>% select(Time, CLiver)   # Liver
Obs.MLu50 = Obs.M50 %>% select(Time, CLung)    # Lung
Obs.MS50  = Obs.M50 %>% select(Time, CSpleen)  # Spleen
Obs.MB50  = Obs.M50 %>% select(Time, CBlood)   # Blood
Obs.MGI50 = Obs.M50 %>% select(Time, CGI)      # GI
Obs.MBr50 = Obs.M50 %>% select(Time, CBrain)   # Brain

Obs.M100 <-read.csv(file="Ma100.csv")      # Ma dataset: 100 nm
Obs.MK100  = Obs.M100 %>% select(Time, CKidney)  # Kidney
Obs.ML100  = Obs.M100 %>% select(Time, CLiver)   # Liver
Obs.MLu100 = Obs.M100 %>% select(Time, CLung)    # Lung
Obs.MS100  = Obs.M100 %>% select(Time, CSpleen)  # Spleen
Obs.MB100  = Obs.M100 %>% select(Time, CBlood)   # Blood
Obs.MGI100 = Obs.M100 %>% select(Time, CGI)      # GI
Obs.MBr100 = Obs.M100 %>% select(Time, CBrain)   # Brain

# Prediction function for Liang's data
pred.eval.L <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 1    # day
  ## Exposure scenario L for Single dose at 0.1 mg
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
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.L, tgrid = tsamp.L)
  out.L <-cbind.data.frame(Time    = out.L$time/24,     # day
                           CSpleen  = out.L$Spleen,
                           CKidney  = out.L$Kidney,
                           CBlood   = out.L$Blood,
                           CLung    = out.L$Lung,
                           CLiver   = out.L$Liver,
                           CGI      = out.L$GI,
                           CBrain   = out.L$Brain
  )
  return(list("out.L"= out.L))
}

# Prediction function for Han's data
pred.eval.H <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24     # h
  End_time      = 90     # day
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

# Prediction function for Wang's data
pred.eval.W <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 1    # day
  ## Exposure scenario W for single dose at 2.5, 25, 250 mg/kg, BW = 20 g
  TDOSE.W          = 1       # Dosing frequency during exposure time
  DOSEoral.Wlow    = 0.05    # mg
  DOSEoral.Wmed    = 0.5     # mg
  DOSEoral.Whigh   = 5       # mg
  ## Oral exposure route
  ex.W1 <- ev(ID = 1, amt = DOSEoral.Wlow, ii = tinterval, addl = TDOSE.W-1, 
              cmt = "ALumen", replicate = FALSE) 
  ex.W2 <- ev(ID = 1, amt = DOSEoral.Wmed, ii = tinterval, addl = TDOSE.W-1, 
              cmt = "ALumen", replicate = FALSE) 
  ex.W3 <- ev(ID = 1, amt = DOSEoral.Whigh, ii = tinterval, addl = TDOSE.W-1, 
              cmt = "ALumen", replicate = FALSE) 
  ## set up the exposure time
  tsamp.W  = tgrid(0, tinterval*(TDOSE.W-1) + tinterval*End_time, 0.1)
  ## Get a prediction
  # For Low Dose Scenario 
  out.W1 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.W1, tgrid = tsamp.W)
  out.W1 <-cbind.data.frame(Time    = out.W1$time/24,   # day
                            CBlood  = out.W1$Blood)
  # For Medium Dose Scenario 
  out.W2 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.W2, tgrid = tsamp.W)
  out.W2 <-cbind.data.frame(Time    = out.W2$time/24,   # day
                            CBlood  = out.W2$Blood)
  # For High Dose Scenario 
  out.W3 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.W3, tgrid = tsamp.W)
  out.W3 <-cbind.data.frame(Time    = out.W3$time/24,   # day
                            CBlood  = out.W3$Blood)
  return(list("out.W1" = out.W1,"out.W2" = out.W2, "out.W3" =out.W3))
}

# Prediction function for Tao's data
pred.eval.T <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 112    # day
  BW            = 0.023  # kg
  ## Exposure scenario T for Single dose at 0.1 mg
  TDOSE.T       = 112     # Dosing frequency during exposure time
  DOSEoral.T    = 0.224   # mg, Daily dose at 224 ug, 5-week-old C57BL/6 male mice (bw 18-26 g)
  ## Oral exposure route
  ex.T <- ev(ID = 1, amt = DOSEoral.T, ii = tinterval, addl = TDOSE.T-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.T  = tgrid(0, tinterval*(TDOSE.T-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  out.T <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.T, tgrid = tsamp.T)
  out.T <-cbind.data.frame(Time     = out.T$time/24,    # day
                           CLiver   = out.T$Liver)
  return(list("out.T"= out.T))
}

# Prediction function for Ma's data
pred.eval.M <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24     # h
  End_time      = 7      # day
  ## Exposure scenario M for Single dose at 0.1 mg
  TDOSE.M       = 7     # Dosing frequency during exposure time
  DOSEoral.M    = 1     # mg, Daily dose at 1 mg, 6-week-old BALB/c male mice (bw 19-20 g)
  ## Oral exposure route
  ex.M <- ev(ID = 1, amt = DOSEoral.M, ii = tinterval, addl = TDOSE.M-1, 
             cmt = "ALumen", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp.M  = tgrid(0, tinterval*(TDOSE.M-1) + tinterval*End_time, 0.1)
  
  ## Get a prediction
  out.M <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.M, tgrid = tsamp.M)
  out.M <-cbind.data.frame(Time     = out.M$time/24,    # day
                           CLiver   = out.M$Liver,
                           CSpleen  = out.M$Spleen,
                           CKidney  = out.M$Kidney,
                           CLung    = out.M$Lung,
                           CBlood   = out.M$Blood,
                           CBrain   = out.M$Brain,
                           CGI      = out.M$GI
  )
  return(list("out.M"= out.M))
}

MCcost.L <-function (pars){
  out <- pred.eval.L (pars)
  cost<- modCost  (model=out$out.L, obs= Obs.LB50, x="Time")
  cost<- modCost  (model=out$out.L, obs= Obs.LL50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LLu50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LS50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LK50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LGI50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LBr50, x="Time", cost=cost)
  return(cost)
}

MCcost.H <-function (pars){
  out <- pred.eval.H (pars)
  cost<- modCost  (model=out$out.H, obs= Obs.HB50, x="Time")
  cost<- modCost  (model=out$out.H, obs= Obs.HL50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HK50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HGI50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HBr50, x="Time", cost=cost)
  return(cost)
}

MCcost.W <-function (pars){
  out <- pred.eval.W (pars)
  cost<- modCost  (model=out$out.W1, obs= Obs.W50L, x="Time")
  cost<- modCost  (model=out$out.W2, obs= Obs.W50M, x="Time", cost=cost)
  cost<- modCost  (model=out$out.W3, obs= Obs.W50H, x="Time", cost=cost)
  return(cost)
}

MCcost.T <-function (pars){
  out <- pred.eval.T (pars)
  cost<- modCost  (model=out$out.T, obs= Obs.T80, x="Time")
  return(cost)
}

MCcost.M1 <-function (pars){
  out <- pred.eval.M (pars)
  cost<- modCost  (model=out$out.M, obs= Obs.ML25, x="Time")
  cost<- modCost  (model=out$out.M, obs= Obs.MGI25, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MS25, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MK25, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MB25, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MLu25, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MBr25, x="Time", cost=cost)
  return(cost)
}
MCcost.M2 <-function (pars){
  out <- pred.eval.M (pars)
  cost<- modCost  (model=out$out.M, obs= Obs.ML50, x="Time")
  cost<- modCost  (model=out$out.M, obs= Obs.MGI50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MS50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MK50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MB50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MLu50, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MBr50, x="Time", cost=cost)
  return(cost)
}
MCcost.M3 <-function (pars){
  out <- pred.eval.M (pars)
  cost<- modCost  (model=out$out.M, obs= Obs.ML100, x="Time")
  cost<- modCost  (model=out$out.M, obs= Obs.MGI100, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MS100, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MK100, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MB100, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MLu100, x="Time", cost=cost)
  cost<- modCost  (model=out$out.M, obs= Obs.MBr100, x="Time", cost=cost)
  return(cost)
}

## Sensitivity function (FME) 
## Check the sensitive parameters in the model based on the optimal parameters derived from model calibration
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

# Local sensitivity analysis for Wang's data
SnsPlasma.W <- sensFun(func = MCcost.W, parms = theta.final, varscale = 1)
Sen.W=summary(SnsPlasma.W)
plot(summary(SnsPlasma.W))
theta.W <- theta.final[abs(Sen.W$Mean) > 1.2*mean(abs(Sen.W$Mean))]  # Selected sensitive parameters
theta.W

# Local sensitivity analysis for Tao's data
SnsPlasma.T <- sensFun(func = MCcost.T, parms = theta.final, varscale = 1)
Sen.T=summary(SnsPlasma.T)
plot(summary(SnsPlasma.T))
theta.T <- theta.final[abs(Sen.T$Mean) > 1.2*mean(abs(Sen.T$Mean))]  # Selected sensitive parameters
theta.T

# Local sensitivity analysis for Ma's data
SnsPlasma.M1 <- sensFun(func = MCcost.M1, parms = theta.final, varscale = 1)
Sen.M=summary(SnsPlasma.M1)
plot(summary(SnsPlasma.M1))
theta.M <- theta.final[abs(Sen.M$Mean) > 1.2*mean(abs(Sen.M$Mean))]  # Selected sensitive parameters
theta.M


## Re-fitted parameters based on Liang's data
# Set up sensitivity or necessary parameter as model input
theta.L <- log(c(PASC = 0.001, PALuC = 0.001, KSRESrelease = 0.003, KGIb = 4.23e-5))
theta.L <- log(c(PS = 5.25, PK = 7.2, KLRESmax = 0.85, KGIb = 1.14e-2, Kfeces = 0.045))  #Derived by several trials and errors
## PBPK model fitting 
Fit.L<- modFit(f=MCcost.L, p=theta.L, method ="Marq", control = nls.lm.control(nprint=1))

theta.L <- log(c(PL = 0.8, PLu = 1, PASC = 0.1, PAKC = 0.1, KGIb = 1.14e-2, Kfeces = 0.045))  
Fit.L<- modFit(f=MCcost.L, p=theta.L, method ="Port", control = nls.lm.control(nprint=1))

summary(Fit.L)                            ## Summary of fit 
exp(Fit.L$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.L(Fit.L$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals 
# Fitting Liang's data
theta.final <- log(c(PS = 5.1645, PK = 7.22, KLRESmax = 1.403,   # Fitted
                     PL = 0.67, PLu = 0.876, PAKC = 0.109, KGIb = 0.0122 ,Kfeces = 0.0437,   # Fitted
                     PBR = 0.48, PASC = 0.012, PABRC = 0.0004, KLuRESmax = 4, KD = 2100))    # Adjust

Sim.fitL = pred.eval.L(theta.final)$out.L     # Simulation 


## Re-fitted parameters based on Han's data
# Set up sensitivity or necessary parameter as model input
theta.H <- log(c(KKRESrelease = 0.01, KGIb = 4.23e-5, KurineC  = 0.00012))
theta.H <- log(c(PK = 0.26, KLRESmax = 0.05, KGIb = 1.6e-4, Kfeces = 0.15))
## PBPK model fitting 
Fit.H<- modFit(f=MCcost.H, p=theta.H, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit.H)                            ## Summary of fit 
exp(Fit.H$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.H(Fit.H$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Han's data
theta.final.H <- log(c(
  BW = 0.015,       # Body weight changed based on Han et al. (2023)
  PK = 0.2757, KLRESmax = 0.0458, Kfeces = 0.1369,    # Fitted
  PBR = 0.56, PABRC = 0.001, KGIb = 0.000165))        # Adjust

Sim.fitH = pred.eval.H(theta.final.H)$out.H   # Simulation 


## Re-fitted parameters based on Wang's data
# Set up sensitivity or necessary parameter as model input
theta.W <- log(c(KGIb = 4.23e-5))
theta.W <- log(c(KGIb = 4.23e-4, Kfeces = 0.02))

## PBPK model fitting 
Fit.W<- modFit(f=MCcost.W, p=theta.W, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit.W)                            ## Summary of fit 
exp(Fit.W$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.W(Fit.W$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Wang's data
theta.final.W <- log(c(KGIb = 0.0001568, Kfeces = 0.01283))   # Fitted

# Simulation 
Sim.fitW1 = pred.eval.W(theta.final.W)$out.W1   
Sim.fitW2 = pred.eval.W(theta.final.W)$out.W2   
Sim.fitW3 = pred.eval.W(theta.final.W)$out.W3 


## Re-fitted parameters based on Tao's data
# Set up sensitivity or necessary parameter as model input
theta.T <- log(c(KGIb = 4.23e-5, KurineC  = 0.00012))
theta.T <- log(c(KGIb = 1.5e-12))   

## PBPK model fitting 
Fit.T<- modFit(f=MCcost.T, p=theta.T, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit.T)                            ## Summary of fit 
exp(Fit.T$par)                            ## Get the arithmetic value out of the log domain

# Fitting Tao's data
theta.final.T <- log(c(
  BW = 0.023,        # Body weight changed based on Tao et al. (2024)
  KGIb = 1.66e-12))  # Adjust

Sim.fitT = pred.eval.T(theta.final.T)$out.T   # Simulation 


## Re-fitted parameters based on Ma's data
# Set up sensitivity or necessary parameter as model input
theta.M <- log(c(PR = 0.15, KSRESrelease = 0.003, KGIb = 4.23e-5)) 
theta.M1 <- log(c(PBR = 0.3, PABRC = 0.003, PS = 0.0005, PK = 0.2, 
                  KGIb = 0.01, Kfeces = 0.11))
## PBPK model fitting 
Fit.M<- modFit(f=MCcost.M1, p=theta.M1, method ="Marq", control = nls.lm.control(nprint=1))

theta.M2 <- log(c(PK = 0.13, PS = 0.001, KGIb = 0.01, Kfeces = 0.135))    
## PBPK model fitting 
Fit.M<- modFit(f=MCcost.M2, p=theta.M2, method ="Marq", control = nls.lm.control(nprint=1))

theta.M3 <- log(c(PK = 0.15, PS = 0.001, KGIb = 0.01, Kfeces = 0.135))
## PBPK model fitting 
Fit.M<- modFit(f=MCcost.M3, p=theta.M3, method ="Marq", control = nls.lm.control(nprint=1))

summary(Fit.M)                            ## Summary of fit 
exp(Fit.M$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.M1(Fit.M$par)$residuals$res    ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Ma's data
# Simulation for 25 nm
theta.final.M1 <- log(c(
  PK = 0.278, PS = 4.33e-05, PBR = 0.288, PABRC = 0.00723, KGIb = 0.01168, Kfeces = 0.1114,  # Fitted
  PR = 0.22))  # Adjust
 
Sim.fitM1 = pred.eval.M(theta.final.M1)$out.M   

# Simulation for 50 nm
theta.final.M2 <- log(c(
  PLu = 0.01, PABRC = 0.000029,    # Adjust
  PK = 0.126, PS = 0.000928, KGIb = 0.01178, Kfeces = 0.133))  # Fitted

Sim.fitM2 = pred.eval.M(theta.final.M2)$out.M   

# Simulation for 100 nm
theta.final.M3 <- log(c(
  PABRC = 0.000024,  PLu = 0.01, KGIb  = 0.009,       # Adjust
  PK    = 0.1617,  PS = 0.0005676, Kfeces = 0.1372))  # Fitted

Sim.fitM3 = pred.eval.M(theta.final.M3)$out.M   

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

df.simHK  = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney)
df.simHL  = cbind.data.frame (Time=Sim.fitH$Time, CLiver=Sim.fitH$CLiver)
df.simHB  = cbind.data.frame (Time=Sim.fitH$Time, CBlood=Sim.fitH$CBlood)
df.simHGI = cbind.data.frame (Time=Sim.fitH$Time, CGI=Sim.fitH$CGI)
df.simHBr = cbind.data.frame (Time=Sim.fitH$Time, CBrain=Sim.fitH$CBrain) 
df.simH   = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney,
                              CLiver=Sim.fitH$CLiver, CBlood=Sim.fitH$CBlood,
                              CGI=Sim.fitH$CGI, CBrain=Sim.fitH$CBrain)

df.simW1  = cbind.data.frame(Time=Sim.fitW1$Time, CBlood=Sim.fitW1$CBlood)
df.simW2  = cbind.data.frame(Time=Sim.fitW2$Time, CBlood=Sim.fitW2$CBlood)
df.simW3  = cbind.data.frame(Time=Sim.fitW3$Time, CBlood=Sim.fitW3$CBlood)
df.simW   = cbind.data.frame(Time=Sim.fitW1$Time, CBloodL=Sim.fitW1$CBlood,
                             CBloodM=Sim.fitW2$CBlood, CBloodH=Sim.fitW3$CBlood) 

df.simT  = cbind.data.frame (Time=Sim.fitT$Time, CLiver=Sim.fitT$CLiver)

df.simM1S  = cbind.data.frame (Time=Sim.fitM1$Time, CSpleen=Sim.fitM1$CSpleen)
df.simM1K  = cbind.data.frame (Time=Sim.fitM1$Time, CKidney=Sim.fitM1$CKidney)
df.simM1L  = cbind.data.frame (Time=Sim.fitM1$Time, CLiver=Sim.fitM1$CLiver)
df.simM1B  = cbind.data.frame (Time=Sim.fitM1$Time, CBlood=Sim.fitM1$CBlood)
df.simM1Lu = cbind.data.frame (Time=Sim.fitM1$Time, CLung=Sim.fitM1$CLung)
df.simM1GI = cbind.data.frame (Time=Sim.fitM1$Time, CGI=Sim.fitM1$CGI)
df.simM1Br = cbind.data.frame (Time=Sim.fitM1$Time, CBrain=Sim.fitM1$CBrain) 
df.simM1   = cbind.data.frame (Time=Sim.fitM1$Time, CSpleen=Sim.fitM1$CSpleen, 
                               CKidney=Sim.fitM1$CKidney, CLiver=Sim.fitM1$CLiver, 
                               CBlood=Sim.fitM1$CBlood, CLung=Sim.fitM1$CLung, 
                               CGI=Sim.fitM1$CGI, CBrain=Sim.fitM1$CBrain) 

df.simM2S  = cbind.data.frame (Time=Sim.fitM2$Time, CSpleen=Sim.fitM2$CSpleen)
df.simM2K  = cbind.data.frame (Time=Sim.fitM2$Time, CKidney=Sim.fitM2$CKidney)
df.simM2L  = cbind.data.frame (Time=Sim.fitM2$Time, CLiver=Sim.fitM2$CLiver)
df.simM2B  = cbind.data.frame (Time=Sim.fitM2$Time, CBlood=Sim.fitM2$CBlood)
df.simM2Lu = cbind.data.frame (Time=Sim.fitM2$Time, CLung=Sim.fitM2$CLung)
df.simM2GI = cbind.data.frame (Time=Sim.fitM2$Time, CGI=Sim.fitM2$CGI)
df.simM2Br = cbind.data.frame (Time=Sim.fitM2$Time, CBrain=Sim.fitM2$CBrain) 
df.simM2   = cbind.data.frame (Time=Sim.fitM2$Time, CSpleen=Sim.fitM2$CSpleen, 
                               CKidney=Sim.fitM2$CKidney, CLiver=Sim.fitM2$CLiver, 
                               CBlood=Sim.fitM2$CBlood, CLung=Sim.fitM2$CLung, 
                               CGI=Sim.fitM2$CGI, CBrain=Sim.fitM2$CBrain) 

df.simM3S  = cbind.data.frame (Time=Sim.fitM3$Time, CSpleen=Sim.fitM3$CSpleen)
df.simM3K  = cbind.data.frame (Time=Sim.fitM3$Time, CKidney=Sim.fitM3$CKidney)
df.simM3L  = cbind.data.frame (Time=Sim.fitM3$Time, CLiver=Sim.fitM3$CLiver)
df.simM3B  = cbind.data.frame (Time=Sim.fitM3$Time, CBlood=Sim.fitM3$CBlood)
df.simM3Lu = cbind.data.frame (Time=Sim.fitM3$Time, CLung=Sim.fitM3$CLung)
df.simM3GI = cbind.data.frame (Time=Sim.fitM3$Time, CGI=Sim.fitM3$CGI)
df.simM3Br = cbind.data.frame (Time=Sim.fitM3$Time, CBrain=Sim.fitM3$CBrain) 
df.simM3   = cbind.data.frame (Time=Sim.fitM3$Time, CSpleen=Sim.fitM3$CSpleen, 
                               CKidney=Sim.fitM3$CKidney, CLiver=Sim.fitM3$CLiver, 
                               CBlood=Sim.fitM3$CBlood, CLung=Sim.fitM3$CLung, 
                               CGI=Sim.fitM3$CGI, CBrain=Sim.fitM3$CBrain) 

write.csv(df.simL, file = "Fitting eval.50_Liang.csv", row.names = FALSE)
write.csv(df.simH, file = "Fitting eval.50_Han.csv", row.names = FALSE)
write.csv(df.simW, file = "Fitting eval.50_Wang.csv", row.names = FALSE)
write.csv(df.simT, file = "Fitting eval.80_Tao.csv", row.names = FALSE)
write.csv(df.simM1, file = "Fitting eval.25_Ma.csv", row.names = FALSE)
write.csv(df.simM2, file = "Fitting eval.50_Ma.csv", row.names = FALSE)
write.csv(df.simM3, file = "Fitting eval.100_Ma.csv", row.names = FALSE)

## Model evaluation plot using ggplot2 
# Liang's results
plot.LB50=
  ggplot() +
  geom_line (data = df.simLB, aes(Time, CBlood), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LB50, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.LS50=
  ggplot() +
  geom_line (data = df.simLS, aes(Time, CSpleen), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LS50, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.LL50=
  ggplot() +
  geom_line (data = df.simLL, aes(Time, CLiver), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LL50, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.LK50=
  ggplot() +
  geom_line (data = df.simLK, aes(Time, CKidney), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LK50, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.LGI50=
  ggplot() +
  geom_line (data = df.simLGI, aes(Time, CGI), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LGI50, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.LLu50=
  ggplot() +
  geom_line (data = df.simLLu, aes(Time, CLung), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LLu50, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

plot.LBr50=
  ggplot() +
  geom_line (data = df.simLBr, aes(Time, CBrain), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LBr50, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("20Fitting.50_Liang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.LB50, plot.LS50, plot.LK50, plot.LL50, plot.LLu50, plot.LGI50, plot.LBr50, nrow = 3)
dev.off()


#Han's results
plot.HB50=
  ggplot() +
  geom_line (data = df.simHB, aes(Time, CBlood), col="blue", lwd=1)+
  geom_point(data = Obs.HB50, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.HL50=
  ggplot() +
  geom_line (data = df.simHL, aes(Time, CLiver), col="blue", lwd=1)+
  geom_point(data = Obs.HL50, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.HK50=
  ggplot() +
  geom_line (data = df.simHK, aes(Time, CKidney), col="blue", lwd=1)+
  geom_point(data = Obs.HK50, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.HGI50=
  ggplot() +
  geom_line (data = df.simHGI, aes(Time, CGI), col="blue", lwd=1)+
  geom_point(data = Obs.HGI50, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.HBr50=
  ggplot() +
  geom_line (data = df.simHBr, aes(Time, CBrain), col="blue", lwd=1)+
  geom_point(data = Obs.HBr50, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("20Fitting.50_Han.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.HB50, plot.HK50, plot.HL50, plot.HGI50, plot.HBr50, nrow = 3)
dev.off()

 
# Wang's results
plot.W50L=
  ggplot() +
  geom_line (data = df.simW1, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.W50L, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/kg)") 

plot.W50M=
  ggplot() +
  geom_line (data = df.simW2, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.W50M, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/kg)")

plot.W50H=
  ggplot() +
  geom_line (data = df.simW3, aes(Time, CBlood), col="firebrick", lwd=1)+
  geom_point(data = Obs.W50H, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/kg)")

# Arrange plots into one figure and save as png file
png("20Fitting.50_Wang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.W50L, plot.W50M, plot.W50H, nrow = 2)
dev.off()


#Tao's results
plot.T80=
  ggplot() +
  geom_line (data = df.simT, aes(Time, CLiver), col="green", lwd=1)+
  geom_point(data = Obs.T80, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 
# Arrange plots into one figure and save as png file
png("20Fitting.80_Tao.png", width = 6600, height = 6000, units = "px", res = 600)
plot.T80
dev.off()


# Ma's results
plot.MB25=
  ggplot() +
  geom_line (data = df.simM1B, aes(Time, CBlood), col="purple", lwd=1)+
  geom_point(data = Obs.MB25, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.MBr25=
  ggplot() +
  geom_line (data = df.simM1Br, aes(Time, CBrain), col="purple", lwd=1)+
  geom_point(data = Obs.MBr25, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

plot.MGI25=
  ggplot() +
  geom_line (data = df.simM1GI, aes(Time, CGI), col="purple", lwd=1)+
  geom_point(data = Obs.MGI25, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.MS25=
  ggplot() +
  geom_line (data = df.simM1S, aes(Time, CSpleen), col="purple", lwd=1)+
  geom_point(data = Obs.MS25, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.ML25=
  ggplot() +
  geom_line (data = df.simM1L, aes(Time, CLiver), col="purple", lwd=1)+
  geom_point(data = Obs.ML25, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.MK25=
  ggplot() +
  geom_line (data = df.simM1K, aes(Time, CKidney), col="purple", lwd=1)+
  geom_point(data = Obs.MK25, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.MLu25=
  ggplot() +
  geom_line (data = df.simM1Lu, aes(Time, CLung), col="purple", lwd=1)+
  geom_point(data = Obs.MLu25, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

# Arrange plots into one figure and save as png file
png("20Fitting.25_Ma.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.MB25, plot.MS25, plot.MGI25, plot.MK25, plot.ML25, plot.MLu25, plot.MBr25, nrow = 3)
dev.off()


plot.MB50=
  ggplot() +
  geom_line (data = df.simM2B, aes(Time, CBlood), col="orange", lwd=1)+
  geom_point(data = Obs.MB50, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.MBr50=
  ggplot() +
  geom_line (data = df.simM2Br, aes(Time, CBrain), col="orange", lwd=1)+
  geom_point(data = Obs.MBr50, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

plot.MGI50=
  ggplot() +
  geom_line (data = df.simM2GI, aes(Time, CGI), col="orange", lwd=1)+
  geom_point(data = Obs.MGI50, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.MS50=
  ggplot() +
  geom_line (data = df.simM2S, aes(Time, CSpleen), col="orange", lwd=1)+
  geom_point(data = Obs.MS50, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.ML50=
  ggplot() +
  geom_line (data = df.simM2L, aes(Time, CLiver), col="orange", lwd=1)+
  geom_point(data = Obs.ML50, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.MK50=
  ggplot() +
  geom_line (data = df.simM2K, aes(Time, CKidney), col="orange", lwd=1)+
  geom_point(data = Obs.MK50, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.MLu50=
  ggplot() +
  geom_line (data = df.simM2Lu, aes(Time, CLung), col="orange", lwd=1)+
  geom_point(data = Obs.MLu50, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

# Arrange plots into one figure and save as png file
png("20Fitting.50_Ma.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.MB50, plot.MS50, plot.MK50, plot.ML50, plot.MLu50, plot.MGI50, plot.MBr50, nrow = 3)
dev.off()


plot.MB100=
  ggplot() +
  geom_line (data = df.simM3B, aes(Time, CBlood), col="navy", lwd=1)+
  geom_point(data = Obs.MB100, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.MBr100=
  ggplot() +
  geom_line (data = df.simM3Br, aes(Time, CBrain), col="navy", lwd=1)+
  geom_point(data = Obs.MBr100, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

plot.MGI100=
  ggplot() +
  geom_line (data = df.simM3GI, aes(Time, CGI), col="navy", lwd=1)+
  geom_point(data = Obs.MGI100, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.MS100=
  ggplot() +
  geom_line (data = df.simM3S, aes(Time, CSpleen), col="navy", lwd=1)+
  geom_point(data = Obs.MS100, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.ML100=
  ggplot() +
  geom_line (data = df.simM3L, aes(Time, CLiver), col="navy", lwd=1)+
  geom_point(data = Obs.ML100, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.MK100=
  ggplot() +
  geom_line (data = df.simM3K, aes(Time, CKidney), col="navy", lwd=1)+
  geom_point(data = Obs.MK100, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.MLu100=
  ggplot() +
  geom_line (data = df.simM3Lu, aes(Time, CLung), col="navy", lwd=1)+
  geom_point(data = Obs.MLu100, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

# Arrange plots into one figure and save as png file
png("20Fitting.100_Ma.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.MB100, plot.MS100, plot.MK100, plot.ML100, plot.MLu100, plot.MGI100, plot.MBr100, nrow = 3)
dev.off()
#----------------------------------------------------------------------------
# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

## Model performance of Liang's results
# Specify the observed x values
specific_Lx_values <- c(1)

# Extract predicted y values for specific x values
selected_LBrows  <- df.simLB[df.simLB$Time %in% specific_Lx_values, ]
selected_LLrows  <- df.simLL[df.simLL$Time %in% specific_Lx_values, ]
selected_LSrows  <- df.simLS[df.simLS$Time %in% specific_Lx_values, ]
selected_LKrows  <- df.simLK[df.simLK$Time %in% specific_Lx_values, ]
selected_LLurows <- df.simLLu[df.simLLu$Time %in% specific_Lx_values, ]
selected_LBRrows <- df.simLBr[df.simLBr$Time %in% specific_Lx_values, ]
selected_LGIrows <- df.simLGI[df.simLGI$Time %in% specific_Lx_values, ]
# Extracted y values
predicted_LBvalues  <- selected_LBrows$CBlood
predicted_LLvalues  <- selected_LLrows$CLiver
predicted_LSvalues  <- selected_LSrows$CSpleen
predicted_LKvalues  <- selected_LKrows$CKidney
predicted_LLuvalues <- selected_LLurows$CLung
predicted_LBRvalues <- selected_LBRrows$CBrain
predicted_LGIvalues <- selected_LGIrows$CGI
# Observed y values
mean.obsLB  <- c(33.5892)
mean.obsLL  <- c(185.338)
mean.obsLS  <- c(540)
mean.obsLK  <- c(500.715)
mean.obsLLu <- c(100.347)
mean.obsLBR <- c(27.0112)
mean.obsLGI <- c(1563.699)

mape_resultLL   <- calculate_mape(mean.obsLL, predicted_LLvalues)
mape_resultLS   <- calculate_mape(mean.obsLS, predicted_LSvalues)
mape_resultLK   <- calculate_mape(mean.obsLK, predicted_LKvalues)
mape_resultLLu  <- calculate_mape(mean.obsLLu, predicted_LLuvalues)
mape_resultLB   <- calculate_mape(mean.obsLB, predicted_LBvalues)
mape_resultLBR  <- calculate_mape(mean.obsLBR, predicted_LBRvalues)
mape_resultLGI  <- calculate_mape(mean.obsLGI, predicted_LGIvalues)

# Print the result
cat("MAPE:", mape_resultLL, "%\n")
cat("MAPE:", mape_resultLS, "%\n")
cat("MAPE:", mape_resultLK, "%\n")
cat("MAPE:", mape_resultLB, "%\n")
cat("MAPE:", mape_resultLLu, "%\n")
cat("MAPE:", mape_resultLBR, "%\n")
cat("MAPE:", mape_resultLGI, "%\n")

# Overall goodness of fit
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
mean.obsHB  <- c(5.97751)
mean.obsHL  <- c(22.026)
mean.obsHK  <- c(8.9173)
mean.obsHBR <- c(3.0053)
mean.obsHGI <- c(35.4827)

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


## Model performance of Wang's results
# Specify the observed x values
specific_Wx_values <- c(1)

# Extract predicted y values for specific x values
selected_Wrows <- df.simW[df.simW$Time %in% specific_Wx_values, ]

# Extracted y values
predicted_W1values  <- selected_Wrows$CBloodL
predicted_W2values  <- selected_Wrows$CBloodM
predicted_W3values  <- selected_Wrows$CBloodH
# Observed y values
mean.obsL  <- c(0.1314)
mean.obsM  <- c(0.4745)
mean.obsH  <- c(3.3431)

mape_resultL   <- calculate_mape(mean.obsL, predicted_W1values)
mape_resultM   <- calculate_mape(mean.obsM, predicted_W2values)
mape_resultH   <- calculate_mape(mean.obsH, predicted_W3values)

# Print the result
cat("MAPE:", mape_resultL, "%\n")
cat("MAPE:", mape_resultM, "%\n")
cat("MAPE:", mape_resultH, "%\n")

# Overall goodness of fit
observed.W  <- c(mean.obsL, mean.obsM, mean.obsH)
predicted.W <- c(predicted_W1values, predicted_W2values, predicted_W3values)

# Fit a linear regression model
lm_model.W <- lm(observed.W ~ predicted.W)

# Calculate R-squared
r_squared.W <- summary(lm_model.W)$r.squared

# Calculate adjusted R-squared
n <- length(observed.W)
p <- length(lm_model.W$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared.W <- 1 - (1 - r_squared.W) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.W <- sqrt(mean((observed.W - predicted.W)^2))

# Print the results
cat("R-squared:", r_squared.W, "\n")
cat("Adjusted R-squared:", adjusted_r_squared.W, "\n")
cat("RMSE:", rmse.W, "\n")


## Model performance of Tao's results
specific_Tx_values <- c(112)                     # Specify the observed x values
# Extract predicted y values for specific x values
selected_TLrows  <- df.simT[df.simT$Time %in% specific_Tx_values, ]
predicted_TLvalues  <- selected_TLrows$CLiver    # Extracted y values
mean.obsTL  <- c(2.51675e-09)                    # Observed y values
mape_resultTL   <- calculate_mape(mean.obsTL, predicted_TLvalues)  # Calculate MAPE 
cat("MAPE:", mape_resultTL, "%\n")               # Print the result

rmse.T <- sqrt(mean((mean.obsTL - predicted_TLvalues)^2))      # Calculate RMSE
# Print the results
cat("RMSE:", rmse.T, "\n")


## Model performance of Ma's results: 25 nm
# Specify the observed x values
specific_Mx_values <- c(7)

# Extract predicted y values for specific x values
selected_M1Brows  <- df.simM1B[df.simM1B$Time %in% specific_Mx_values, ]
selected_M1Lrows  <- df.simM1L[df.simM1L$Time %in% specific_Mx_values, ]
selected_M1Srows  <- df.simM1S[df.simM1S$Time %in% specific_Mx_values, ]
selected_M1Krows  <- df.simM1K[df.simM1K$Time %in% specific_Mx_values, ]
selected_M1Lurows <- df.simM1Lu[df.simM1Lu$Time %in% specific_Mx_values, ]
selected_M1BRrows <- df.simM1Br[df.simM1Br$Time %in% specific_Mx_values, ]
selected_M1GIrows <- df.simM1GI[df.simM1GI$Time %in% specific_Mx_values, ]
# Extracted y values
predicted_M1Bvalues  <- selected_M1Brows$CBlood
predicted_M1Lvalues  <- selected_M1Lrows$CLiver
predicted_M1Svalues  <- selected_M1Srows$CSpleen
predicted_M1Kvalues  <- selected_M1Krows$CKidney
predicted_M1Luvalues <- selected_M1Lurows$CLung
predicted_M1BRvalues <- selected_M1BRrows$CBrain
predicted_M1GIvalues <- selected_M1GIrows$CGI
# Observed y values
mean.obsM1B  <- c(144.9977)
mean.obsM1L  <- c(70.31077)
mean.obsM1S  <- c(8.304)
mean.obsM1K  <- c(57.70215)
mean.obsM1Lu <- c(57.61471)
mean.obsM1BR <- c(29.31313)
mean.obsM1GI <- c(82.19203)

mape_resultM1L   <- calculate_mape(mean.obsM1L, predicted_M1Lvalues)
mape_resultM1S   <- calculate_mape(mean.obsM1S, predicted_M1Svalues)
mape_resultM1K   <- calculate_mape(mean.obsM1K, predicted_M1Kvalues)
mape_resultM1Lu  <- calculate_mape(mean.obsM1Lu, predicted_M1Luvalues)
mape_resultM1B   <- calculate_mape(mean.obsM1B, predicted_M1Bvalues)
mape_resultM1BR  <- calculate_mape(mean.obsM1BR, predicted_M1BRvalues)
mape_resultM1GI  <- calculate_mape(mean.obsM1GI, predicted_M1GIvalues)

# Print the result
cat("MAPE:", mape_resultM1L, "%\n")
cat("MAPE:", mape_resultM1S, "%\n")
cat("MAPE:", mape_resultM1K, "%\n")
cat("MAPE:", mape_resultM1B, "%\n")
cat("MAPE:", mape_resultM1Lu, "%\n")
cat("MAPE:", mape_resultM1BR, "%\n")
cat("MAPE:", mape_resultM1GI, "%\n")

# Overall goodness of fit
observed.M1  <- c(mean.obsM1L, mean.obsM1S, mean.obsM1K, mean.obsM1Lu, 
                  mean.obsM1B, mean.obsM1BR, mean.obsM1GI)
predicted.M1 <- c(predicted_M1Lvalues, predicted_M1Svalues, predicted_M1Kvalues, 
                  predicted_M1Luvalues, predicted_M1Bvalues, 
                  predicted_M1BRvalues, predicted_M1GIvalues)

# Fit a linear regression model
lm_model.M1 <- lm(observed.M1 ~ predicted.M1)

# Calculate R-squared
r_squared <- summary(lm_model.M1)$r.squared

# Calculate adjusted R-squared
n <- length(observed.M1)
p <- length(lm_model.M1$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.M1 <- sqrt(mean((observed.M1 - predicted.M1)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.M1, "\n")


## Model performance of Ma's results: 50 nm 
# Extract predicted y values for specific x values
selected_M2Brows  <- df.simM2B[df.simM2B$Time %in% specific_Mx_values, ]
selected_M2Lrows  <- df.simM2L[df.simM2L$Time %in% specific_Mx_values, ]
selected_M2Srows  <- df.simM2S[df.simM2S$Time %in% specific_Mx_values, ]
selected_M2Krows  <- df.simM2K[df.simM2K$Time %in% specific_Mx_values, ]
selected_M2Lurows <- df.simM2Lu[df.simM2Lu$Time %in% specific_Mx_values, ]
selected_M2BRrows <- df.simM2Br[df.simM2Br$Time %in% specific_Mx_values, ]
selected_M2GIrows <- df.simM2GI[df.simM2GI$Time %in% specific_Mx_values, ]
# Extracted y values
predicted_M2Bvalues  <- selected_M2Brows$CBlood
predicted_M2Lvalues  <- selected_M2Lrows$CLiver
predicted_M2Svalues  <- selected_M2Srows$CSpleen
predicted_M2Kvalues  <- selected_M2Krows$CKidney
predicted_M2Luvalues <- selected_M2Lurows$CLung
predicted_M2BRvalues <- selected_M2BRrows$CBrain
predicted_M2GIvalues <- selected_M2GIrows$CGI
# Observed y values
mean.obsM2B  <- c(125.4455)
mean.obsM2L  <- c(64.96882)
mean.obsM2S  <- c(23.6)
mean.obsM2K  <- c(47.82213)
mean.obsM2Lu <- c(38.66253)
mean.obsM2BR <- c(12.83518)
mean.obsM2GI <- c(54.60191)

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultM2L   <- calculate_mape(mean.obsM2L, predicted_M2Lvalues)
mape_resultM2S   <- calculate_mape(mean.obsM2S, predicted_M2Svalues)
mape_resultM2K   <- calculate_mape(mean.obsM2K, predicted_M2Kvalues)
mape_resultM2Lu  <- calculate_mape(mean.obsM2Lu, predicted_M2Luvalues)
mape_resultM2B   <- calculate_mape(mean.obsM2B, predicted_M2Bvalues)
mape_resultM2BR  <- calculate_mape(mean.obsM2BR, predicted_M2BRvalues)
mape_resultM2GI  <- calculate_mape(mean.obsM2GI, predicted_M2GIvalues)

# Print the result
cat("MAPE:", mape_resultM2L, "%\n")
cat("MAPE:", mape_resultM2S, "%\n")
cat("MAPE:", mape_resultM2K, "%\n")
cat("MAPE:", mape_resultM2B, "%\n")
cat("MAPE:", mape_resultM2Lu, "%\n")
cat("MAPE:", mape_resultM2BR, "%\n")
cat("MAPE:", mape_resultM2GI, "%\n")

# Overall goodness of fit
observed.M2  <- c(mean.obsM2L, mean.obsM2S, mean.obsM2K, mean.obsM2Lu, 
                  mean.obsM2B, mean.obsM2BR, mean.obsM2GI)
predicted.M2 <- c(predicted_M2Lvalues, predicted_M2Svalues, predicted_M2Kvalues, 
                  predicted_M2Luvalues, predicted_M2Bvalues, 
                  predicted_M2BRvalues, predicted_M2GIvalues)

# Fit a linear regression model
lm_model.M2 <- lm(observed.M2 ~ predicted.M2)

# Calculate R-squared
r_squared <- summary(lm_model.M2)$r.squared

# Calculate Root Mean Squared Error (RMSE)
rmse.M2 <- sqrt(mean((observed.M2 - predicted.M2)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.M2, "\n")


## Model performance of Ma's results: 100 nm 
# Extract predicted y values for specific x values
selected_M3Brows  <- df.simM3B[df.simM3B$Time %in% specific_Mx_values, ]
selected_M3Lrows  <- df.simM3L[df.simM3L$Time %in% specific_Mx_values, ]
selected_M3Srows  <- df.simM3S[df.simM3S$Time %in% specific_Mx_values, ]
selected_M3Krows  <- df.simM3K[df.simM3K$Time %in% specific_Mx_values, ]
selected_M3Lurows <- df.simM3Lu[df.simM3Lu$Time %in% specific_Mx_values, ]
selected_M3BRrows <- df.simM3Br[df.simM3Br$Time %in% specific_Mx_values, ]
selected_M3GIrows <- df.simM3GI[df.simM3GI$Time %in% specific_Mx_values, ]
# Extracted y values
predicted_M3Bvalues  <- selected_M3Brows$CBlood
predicted_M3Lvalues  <- selected_M3Lrows$CLiver
predicted_M3Svalues  <- selected_M3Srows$CSpleen
predicted_M3Kvalues  <- selected_M3Krows$CKidney
predicted_M3Luvalues <- selected_M3Lurows$CLung
predicted_M3BRvalues <- selected_M3BRrows$CBrain
predicted_M3GIvalues <- selected_M3GIrows$CGI
# Observed y values
mean.obsM3B  <- c(122.808)
mean.obsM3L  <- c(52.61458)
mean.obsM3S  <- c(13.69966)
mean.obsM3K  <- c(43.30394)
mean.obsM3Lu <- c(19.29897)
mean.obsM3BR <- c(9.142715)
mean.obsM3GI <- c(49.27785)

mape_resultM3L   <- calculate_mape(mean.obsM3L, predicted_M3Lvalues)
mape_resultM3S   <- calculate_mape(mean.obsM3S, predicted_M3Svalues)
mape_resultM3K   <- calculate_mape(mean.obsM3K, predicted_M3Kvalues)
mape_resultM3Lu  <- calculate_mape(mean.obsM3Lu, predicted_M3Luvalues)
mape_resultM3B   <- calculate_mape(mean.obsM3B, predicted_M3Bvalues)
mape_resultM3BR  <- calculate_mape(mean.obsM3BR, predicted_M3BRvalues)
mape_resultM3GI  <- calculate_mape(mean.obsM3GI, predicted_M3GIvalues)

# Print the result
cat("MAPE:", mape_resultM3L, "%\n")
cat("MAPE:", mape_resultM3S, "%\n")
cat("MAPE:", mape_resultM3K, "%\n")
cat("MAPE:", mape_resultM3B, "%\n")
cat("MAPE:", mape_resultM3Lu, "%\n")
cat("MAPE:", mape_resultM3BR, "%\n")
cat("MAPE:", mape_resultM3GI, "%\n")

# Overall goodness of fit
observed.M3  <- c(mean.obsM3L, mean.obsM3S, mean.obsM3K, mean.obsM3Lu, 
                  mean.obsM3B, mean.obsM3BR, mean.obsM3GI)
predicted.M3 <- c(predicted_M3Lvalues, predicted_M3Svalues, predicted_M3Kvalues, 
                  predicted_M3Luvalues, predicted_M3Bvalues, 
                  predicted_M3BRvalues, predicted_M3GIvalues)

# Fit a linear regression model
lm_model.M3 <- lm(observed.M3 ~ predicted.M3)

# Calculate R-squared
r_squared <- summary(lm_model.M3)$r.squared

# Calculate adjusted R-squared
n <- length(observed.M3)
p <- length(lm_model.M3$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared <- 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.M3 <- sqrt(mean((observed.M3 - predicted.M3)^2))

# Print the results
cat("R-squared:", r_squared, "\n")
cat("Adjusted R-squared:", adjusted_r_squared, "\n")
cat("RMSE:", rmse.M3, "\n")