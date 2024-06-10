## Load libraries
library(mrgsolve)     # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)     # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)      # Needed for plot
library(gridExtra)    # Arrange plots in to one figure
library(FME)          # Package for MCMC simulation and model fitting
library(minpack.lm)   # Package for model fitting
library(reshape)      # Package for melt function to reshape the table
library(tidyr)        # R-package for tidy messy data
library(tidyverse)    # R-package for tidy messy data

## Build mrgsolve-based PBTK Model
mod <- mcode ("NPsPBTK.code", PBTK.code.eval.220) 

## RE-calibrate the model
# Input evaluation data
# Liang et al. (2021) : Single dose at 250 mg/kg, 500 nm, matrix: Blood, GI (Unit: mg/kg) 
# Choi et al. (2021) : Repeated daily dose at 0.005, 0.025, 0.05 mg, 500 nm, matrix: Kidney, Liver, GI (Unit: mg/kg) 
# Han et al. (2023) : Repeated daily dose at 50 mg/kg, 500 nm, matrix: Serum, brain, liver, kidney, intestine (Unit: mg/kg)
# Tao et al. (2024) : Repeated daily dose at 224 ug, 500 nm, matrix: Liver (Unit: mg/kg)
# Zhang et al. (2024) : Repeated daily dose at 0.5 mg, 500 nm, matrix: Liver, Kidney, Spleen, Lung (Unit: mg/kg)
#---------------------------------------------------------------------------#
## Input data set
Obs.L500 <-read.csv(file="Liang500.csv")     # Liang dataset 
Obs.LB500  = Obs.L500 %>% select(Time, CBlood)   # Blood
Obs.LGI500 = Obs.L500 %>% select(Time, CGI)      # GI
Obs.LLu500 = Obs.L500 %>% select(Time, CLung)    # Lung
Obs.LS500  = Obs.L500 %>% select(Time, CSpleen)  # Spleen
Obs.LK500  = Obs.L500 %>% select(Time, CKidney)  # Kidney
Obs.LL500  = Obs.L500 %>% select(Time, CLiver)   # Liver
Obs.LBr500 = Obs.L500 %>% select(Time, CBrain)   # Brain

Obs.C500L <-read.csv(file="Choi500L.csv")         # Choi dataset-Low dose
Obs.CGI500L = Obs.C500L %>% select(Time, CGI)     # GI
Obs.CK500L  = Obs.C500L %>% select(Time, CKidney) # Kidney
Obs.CL500L  = Obs.C500L %>% select(Time, CLiver)  # Liver

Obs.C500M <-read.csv(file="Choi500M.csv")         # Choi dataset-Medium dose
Obs.CGI500M = Obs.C500M %>% select(Time, CGI)     # GI
Obs.CK500M  = Obs.C500M %>% select(Time, CKidney) # Kidney
Obs.CL500M  = Obs.C500M %>% select(Time, CLiver)  # Liver

Obs.C500H <-read.csv(file="Choi500H.csv")         # Choi dataset-High dose
Obs.CGI500H = Obs.C500H %>% select(Time, CGI)     # GI
Obs.CK500H  = Obs.C500H %>% select(Time, CKidney) # Kidney
Obs.CL500H  = Obs.C500H %>% select(Time, CLiver)  # Liver

Obs.H500 <-read.csv(file="Han500.csv")           # Han dataset 
Obs.HB500  = Obs.H500 %>% select(Time, CBlood)   # Blood
Obs.HGI500 = Obs.H500 %>% select(Time, CGI)      # GI
Obs.HK500  = Obs.H500 %>% select(Time, CKidney)  # Kidney
Obs.HL500  = Obs.H500 %>% select(Time, CLiver)   # Liver
Obs.HBr500 = Obs.H500 %>% select(Time, CBrain)   # Brain

Obs.T500 <-read.csv(file="Tao500.csv")           # Tao dataset 
Obs.TL500  = Obs.T500 %>% select(Time, CLiver)   # Liver

Obs.Z500 <-read.csv(file="Zhang500.csv")         # Zhang dataset 
Obs.ZK500  = Obs.Z500 %>% select(Time, CKidney)  # Kidney
Obs.ZL500  = Obs.Z500 %>% select(Time, CLiver)   # Liver
Obs.ZLu500 = Obs.Z500 %>% select(Time, CLung)    # Lung
Obs.ZS500  = Obs.Z500 %>% select(Time, CSpleen)  # Spleen

# Prediction function for Liang's data
pred.eval.L <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 1     # day
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

# Prediction function for Choi's data
pred.eval.C <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24    # h
  End_time      = 14    # day
  BW            = 0.03  # kg
  ## Exposure scenario C for daily dose at 0.005, 0.025, and 0.05 mg, BW =30 g
  TDOSE.C          = 14       # Dosing frequency during exposure time
  DOSEoral.Clow    = 0.005    # mg, daily dose at 0.005 mg
  DOSEoral.Cmed    = 0.025    # mg, daily dose at 0.025 mg
  DOSEoral.Chigh   = 0.05     # mg, daily dose at 0.05 mg
  ## Oral exposure route
  ex.C1 <- ev(ID = 1, amt = DOSEoral.Clow, ii = tinterval, addl = TDOSE.C-1, 
              cmt = "ALumen", replicate = FALSE) 
  ex.C2 <- ev(ID = 1, amt = DOSEoral.Cmed, ii = tinterval, addl = TDOSE.C-1, 
              cmt = "ALumen", replicate = FALSE) 
  ex.C3 <- ev(ID = 1, amt = DOSEoral.Chigh, ii = tinterval, addl = TDOSE.C-1, 
              cmt = "ALumen", replicate = FALSE) 
  ## set up the exposure time
  tsamp.C  = tgrid(0, tinterval*(TDOSE.C-1) + tinterval*End_time, 0.1)
  ## Get a prediction
  # For Scenario C
  out.C1 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.C1, tgrid = tsamp.C)
  out.C1 <-cbind.data.frame(Time    = out.C1$time/24,   # day
                            CKidney  = out.C1$Kidney,
                            CLiver   = out.C1$Liver,
                            CGI      = out.C1$GI
  )
  out.C2 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.C2, tgrid = tsamp.C)
  out.C2 <-cbind.data.frame(Time    = out.C2$time/24,   # day
                            CKidney  = out.C2$Kidney,
                            CLiver   = out.C2$Liver,
                            CGI      = out.C2$GI
  )
  out.C3 <- 
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.C3, tgrid = tsamp.C)
  out.C3 <-cbind.data.frame(Time    = out.C3$time/24,   # day
                            CKidney  = out.C3$Kidney,
                            CLiver   = out.C3$Liver,
                            CGI      = out.C3$GI
  )
  return(list("out.C1" = out.C1,"out.C2" = out.C2, "out.C3" =out.C3))
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

# Prediction function for Tao's data
pred.eval.T <- function(pars) {
  ## Get out of log domain
  pars %<>% lapply(exp)           ## return a list of exp (parameters) from log domain
  
  ## Define the event data sets (exposure scenario) [mrgsolve CH3.2]
  tinterval     = 24     # h
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
                           CLiver   = out.T$Liver
  )
  return(list("out.T"= out.T))
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

MCcost.L <-function (pars){
  out <- pred.eval.L (pars)
  cost<- modCost  (model=out$out.L, obs= Obs.LB500, x="Time")
  cost<- modCost  (model=out$out.L, obs= Obs.LL500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LLu500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LS500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LK500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LGI500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.L, obs= Obs.LBr500, x="Time", cost=cost)
  return(cost)
}

MCcost.C <-function (pars){
  out <- pred.eval.C (pars)
  cost<- modCost  (model=out$out.C1, obs= Obs.CL500L, x="Time")
  cost<- modCost  (model=out$out.C1, obs= Obs.CK500L, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C1, obs= Obs.CGI500L, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C2, obs= Obs.CL500M, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C2, obs= Obs.CK500M, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C2, obs= Obs.CGI500M, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C3, obs= Obs.CL500H, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C3, obs= Obs.CK500H, x="Time", cost=cost)
  cost<- modCost  (model=out$out.C3, obs= Obs.CGI500H, x="Time", cost=cost)
  return(cost)
}

MCcost.H <-function (pars){
  out <- pred.eval.H (pars)
  cost<- modCost  (model=out$out.H, obs= Obs.HB500, x="Time")
  cost<- modCost  (model=out$out.H, obs= Obs.HL500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HK500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HGI500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.H, obs= Obs.HBr500, x="Time", cost=cost)
  return(cost)
}

MCcost.T <-function (pars){
  out <- pred.eval.T (pars)
  cost<- modCost  (model=out$out.T, obs= Obs.TL500, x="Time")
  return(cost)
}

MCcost.Z <-function (pars){
  out <- pred.eval.Z (pars)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZL500, x="Time")
  cost<- modCost  (model=out$out.Z, obs= Obs.ZLu500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZS500, x="Time", cost=cost)
  cost<- modCost  (model=out$out.Z, obs= Obs.ZK500, x="Time", cost=cost)
  return(cost)
}

## Sensitivity function (FME) 
## Check the sensitive parameters in the model based on the optimal parameters derived from model calibration
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

# Local sensitivity analysis for Liang's data
SnsPlasma.L <- sensFun(func = MCcost.L, parms = theta.final, varscale = 1)
Sen.L=summary(SnsPlasma.L)
plot(summary(SnsPlasma.L))
theta.L <- theta.final[abs(Sen.L$Mean) > 1.2*mean(abs(Sen.L$Mean))]  # Selected sensitive parameters
theta.L

# Local sensitivity analysis for Choi's data
SnsPlasma.C <- sensFun(func = MCcost.C, parms = theta.final, varscale = 1)
Sen.C=summary(SnsPlasma.C)
plot(summary(SnsPlasma.C))
theta.C <- theta.final[abs(Sen.C$Mean) > 1.2*mean(abs(Sen.C$Mean))]  # Selected sensitive parameters
theta.C

# Local sensitivity analysis for Han's data
SnsPlasma.H <- sensFun(func = MCcost.H, parms = theta.final, varscale = 1)
Sen.H=summary(SnsPlasma.H)
plot(summary(SnsPlasma.H))
theta.H <- theta.final[abs(Sen.H$Mean) > 1.2*mean(abs(Sen.H$Mean))]  # Selected sensitive parameters
theta.H

# Local sensitivity analysis for Tao's data
SnsPlasma.T <- sensFun(func = MCcost.T, parms = theta.final, varscale = 1)
Sen.T=summary(SnsPlasma.T)
plot(summary(SnsPlasma.T))
theta.T <- theta.final[abs(Sen.T$Mean) > 1.2*mean(abs(Sen.T$Mean))]  # Selected sensitive parameters
theta.T

# Local sensitivity analysis for Zhang's data
SnsPlasma.Z <- sensFun(func = MCcost.Z, parms = theta.final, varscale = 1)
Sen.Z=summary(SnsPlasma.Z)
plot(summary(SnsPlasma.Z))
theta.Z <- theta.final[abs(Sen.Z$Mean) > 1.2*mean(abs(Sen.Z$Mean))]  # Selected sensitive parameters
theta.Z


## Re-fitted parameters based on Liang's data
# Set up sensitivity or necessary parameter as model input
theta.L <- log(c(PASC = 0.00249, KGIb = 4.22e-5, Kfeces = 0.55))
theta.L <- log(c(PBR= 0.8, PK = 2.5, PAKC =0.01, KLRESmax = 0.66,
                 KGIb = 2.75e-3, Kfeces = 0.05))  #Derived by several trials and errors
## PBPK model fitting 
Fit.L<- modFit(f=MCcost.L, p=theta.L, method ="Port", control = nls.lm.control(nprint=1))
summary(Fit.L)                            ## Summary of fit 
exp(Fit.L$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.L(Fit.L$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals 
# Fitting Liang's data
theta.final <- log(c(PLu = 0.0001, PABRC = 0.1,   # Adjusted
                     PBR = 0.791, PK = 2.486, PAKC = 0.01, KLRESmax = 0.6515, # Fitted
                     KGIb = 0.00274, Kfeces = 4.774e-02))                     # Fitted
Sim.fitL = pred.eval.L(theta.final)$out.L         # Simulation 


## Re-fitted parameters based on Choi's data
# Set up sensitivity or necessary parameter as model input
theta.C <- log(c(PK = 0.251, PASC = 0.00249, KKRESrelease = 0.01, KGIb = 4.22e-5, Kfeces = 0.55, KurineC = 0.0003))
theta.C <- log(c(PK = 1.1, KLRESrelease = 0.5, KGIb = 4e-3, Kfeces = 0.03))
# Set boundary for each parameter
# lb <- log(c(PAKC = 1e-2, KLRESrelease = 1e-1, KGIb = 1e-3, Kfeces = 1e-3))
ub <- log(c(PK = 1.4, KLRESrelease = 0.8, KGIb = 4.4e-3, Kfeces = 0.0325))

## PBPK model fitting 
Fit.C<- modFit(f=MCcost.C, p=theta.C, method ="Marq", upper = ub,
               control = nls.lm.control(nprint=1))
summary(Fit.C)                            ## Summary of fit 
exp(Fit.C$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.C(Fit.C$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Choi's data
theta.final.C <- log(c(
  BW = 0.03,     # Body weight changed based on Choi et al. (2021)
  PK = 1.221, KLRESrelease = 0.3711, KGIb = 0.004282, Kfeces = 0.0324  # Fitted
  ))

# Simulation 
Sim.fitC1 = pred.eval.C(theta.final.C)$out.C1   
Sim.fitC2 = pred.eval.C(theta.final.C)$out.C2   
Sim.fitC3 = pred.eval.C(theta.final.C)$out.C3   


## Re-fitted parameters based on Han's data
# Set up sensitivity or necessary parameter as model input
theta.H <- log(c(PK = 0.251, KSRESrelease  = 0.003, KKRESrelease  = 0.01,
                 KGIb = 4.22e-5, KurineC  = 0.0003))
theta.H <- log(c(PBR = 0.001, KLRESmax = 0.1, KGIb = 1e-4, Kfeces = 0.15))
## PBPK model fitting 
Fit.H <- modFit(f=MCcost.H, p=theta.H, method ="Marq", control = nls.lm.control(nprint=1))
summary(Fit.H)                            ## Summary of fit 
exp(Fit.H$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.H(Fit.H$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Han's data
theta.final.H <- log(c(
  BW        = 0.015,     # Body weight changed based on Han et al. (2023)
  PBR       = 0.02,  KGIb   = 1e-04,     # Adjusted
  KLRESmax  = 0.113, Kfeces = 0.1562     # Fitted
))

Sim.fitH = pred.eval.H(theta.final.H)$out.H   # Simulation 


## Re-fitted parameters based on Tao's data
# Set up sensitivity or necessary parameter as model input
theta.T <- log(c(KSRESrelease = 0.003, KLRESrelease = 0.0075, KLRESmax = 0.00427,
                 KGIb = 4.22e-5, KurineC = 0.0003))
theta.T <- log(c(KGIb = 2e-10))
## PBPK model fitting 
Fit.T<- modFit(f=MCcost.T, p=theta.T, method ="Nelder-Mead", control = nls.lm.control(nprint=1))
summary(Fit.T)                            ## Summary of fit 
exp(Fit.T$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.T(Fit.T$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Tao's data
theta.final.T <- log(c(
  BW   = 0.023,       # Body weight changed based on Tao et al. (2024)
  KGIb = 2.37e-10))    # Adjust   
  
Sim.fitT = pred.eval.T(theta.final.T)$out.T  # Simulation  


## Re-fitted parameters based on Zhang's data
# Set up sensitivity or necessary parameter as model input
theta.Z <- log(c(PASC = 0.00249, KSRESrelease = 0.003, 
                 KGIb = 4.22e-5, KurineC = 0.0003)) 
theta.Z <- log(c(PK = 1.9, PS = 1, KGIb = 0.33))  #1.8856980 0.8741364 0.2894973 
theta.Z <- log(c(PLu = 0.0001, PK = 1.9, KGIb = 0.33))  #2.556521e-06 1.700353e+00 3.304154e-01 
## PBPK model fitting 
Fit.Z<- modFit(f=MCcost.Z, p=theta.Z, method ="Marq", control = nls.lm.control(nprint=1))

summary(Fit.Z)                            ## Summary of fit 
exp(Fit.Z$par)                            ## Get the arithmetic value out of the log domain
res=MCcost.Z(Fit.Z$par)$residuals$res     ## Check the residual for each time points
sum(res^2)                                ## Total residuals  

# Fitting Zhang's data
theta.final.Z <- log(c(
  BW   = 0.024,    # Body weight changed based on Zhang et al. (2024)
  PS   = 0.874, PK = 1.8857, KGIb = 0.33,  # Fitted
  PLu  = 0.0001, KLRESrelease = 0.0001     # Adjusted
  ))

Sim.fitZ = pred.eval.Z(theta.final.Z)$out.Z   # Simulation 


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

df.simC1K  = cbind.data.frame(Time=Sim.fitC1$Time, CKidney=Sim.fitC1$CKidney)
df.simC1L  = cbind.data.frame(Time=Sim.fitC1$Time, CLiver=Sim.fitC1$CLiver)
df.simC1GI = cbind.data.frame(Time=Sim.fitC1$Time, CGI=Sim.fitC1$CGI)
df.simC2K  = cbind.data.frame(Time=Sim.fitC2$Time, CKidney=Sim.fitC2$CKidney)
df.simC2L  = cbind.data.frame(Time=Sim.fitC2$Time, CLiver=Sim.fitC2$CLiver)
df.simC2GI = cbind.data.frame(Time=Sim.fitC2$Time, CGI=Sim.fitC2$CGI)
df.simC3K  = cbind.data.frame(Time=Sim.fitC3$Time, CKidney=Sim.fitC3$CKidney)
df.simC3L  = cbind.data.frame(Time=Sim.fitC3$Time, CLiver=Sim.fitC3$CLiver)
df.simC3GI = cbind.data.frame(Time=Sim.fitC3$Time, CGI=Sim.fitC3$CGI)
df.simC    = cbind.data.frame (Time=Sim.fitC1$Time, CKidneyL=Sim.fitC1$CKidney, CLiverL=Sim.fitC1$CLiver,
                               CGIL=Sim.fitC1$CGI, CKidneyM=Sim.fitC2$CKidney, CLiverM=Sim.fitC2$CLiver,
                               CGIM=Sim.fitC2$CGI, CKidneyH=Sim.fitC3$CKidney, 
                               CLiverH=Sim.fitC3$CLiver, CGIH=Sim.fitC3$CGI) 

df.simHK  = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney)
df.simHL  = cbind.data.frame (Time=Sim.fitH$Time, CLiver=Sim.fitH$CLiver)
df.simHB  = cbind.data.frame (Time=Sim.fitH$Time, CBlood=Sim.fitH$CBlood)
df.simHGI = cbind.data.frame (Time=Sim.fitH$Time, CGI=Sim.fitH$CGI)
df.simHBr = cbind.data.frame (Time=Sim.fitH$Time, CBrain=Sim.fitH$CBrain) 
df.simH   = cbind.data.frame (Time=Sim.fitH$Time, CKidney=Sim.fitH$CKidney,
                              CLiver=Sim.fitH$CLiver, CBlood=Sim.fitH$CBlood,
                              CGI=Sim.fitH$CGI, CBrain=Sim.fitH$CBrain)

df.simTL  = cbind.data.frame (Time=Sim.fitT$Time, CLiver=Sim.fitT$CLiver)

df.simZK  = cbind.data.frame (Time=Sim.fitZ$Time, CKidney=Sim.fitZ$CKidney)
df.simZL  = cbind.data.frame (Time=Sim.fitZ$Time, CLiver=Sim.fitZ$CLiver)
df.simZS  = cbind.data.frame (Time=Sim.fitZ$Time, CSpleen=Sim.fitZ$CSpleen)
df.simZLu = cbind.data.frame (Time=Sim.fitZ$Time, CLung=Sim.fitZ$CLung)
df.simZ   = cbind.data.frame (Time=Sim.fitZ$Time, CKidney=Sim.fitZ$CKidney,
                              CLiver=Sim.fitZ$CLiver, CSpleen=Sim.fitZ$CSpleen,
                              CLung=Sim.fitZ$CLung)


write.csv(df.simL, file = "Fitting eval.500_Liang.csv", row.names = FALSE)
write.csv(df.simC, file = "Fitting eval.500_Choi.csv", row.names = FALSE)
write.csv(df.simH, file = "Fitting eval.500_Han.csv", row.names = FALSE)
write.csv(df.simTL, file = "Fitting eval.500_Tao.csv", row.names = FALSE)
write.csv(df.simZ, file = "Fitting eval.500_Zhang.csv", row.names = FALSE)

## Model evaluation plot using ggplot2 
# Liang's results
plot.LB500=
  ggplot() +
  geom_line (data = df.simLB, aes(Time, CBlood), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LB500, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.LS500=
  ggplot() +
  geom_line (data = df.simLS, aes(Time, CSpleen), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LS500, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.LL500=
  ggplot() +
  geom_line (data = df.simLL, aes(Time, CLiver), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LL500, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.LK500=
  ggplot() +
  geom_line (data = df.simLK, aes(Time, CKidney), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LK500, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.LGI500=
  ggplot() +
  geom_line (data = df.simLGI, aes(Time, CGI), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LGI500, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.LLu500=
  ggplot() +
  geom_line (data = df.simLLu, aes(Time, CLung), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LLu500, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

plot.LBr500=
  ggplot() +
  geom_line (data = df.simLBr, aes(Time, CBrain), col="darkgreen", lwd=1)+
  geom_point(data = Obs.LBr500, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("220Fitting.500_Liang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.LB500, plot.LS500, plot.LK500, plot.LL500, plot.LLu500, 
             plot.LGI500, plot.LBr500, nrow = 3)
dev.off()

# Choi's results
plot.CL500L=
  ggplot() +
  geom_line (data = df.simC1L, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.CL500L, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.CK500L=
  ggplot() +
  geom_line (data = df.simC1K, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.CK500L, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.CGI500L=
  ggplot() +
  geom_line (data = df.simC1GI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.CGI500L, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.CL500M=
  ggplot() +
  geom_line (data = df.simC2L, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.CL500M, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.CK500M=
  ggplot() +
  geom_line (data = df.simC2K, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.CK500M, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.CGI500M=
  ggplot() +
  geom_line (data = df.simC2GI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.CGI500M, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)")

plot.CL500H=
  ggplot() +
  geom_line (data = df.simC3L, aes(Time, CLiver), col="firebrick", lwd=1)+
  geom_point(data = Obs.CL500H, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.CK500H=
  ggplot() +
  geom_line (data = df.simC3K, aes(Time, CKidney), col="firebrick", lwd=1)+
  geom_point(data = Obs.CK500H, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.CGI500H=
  ggplot() +
  geom_line (data = df.simC3GI, aes(Time, CGI), col="firebrick", lwd=1)+
  geom_point(data = Obs.CGI500H, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)")

# Arrange plots into one figure and save as png file
png("220Fitting.500L_Choi.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.CL500L, plot.CK500L, plot.CGI500L, nrow = 2)
dev.off()

png("220Fitting.500M_Choi.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.CL500M, plot.CK500M, plot.CGI500M, nrow = 2)
dev.off()

png("220Fitting.500H_Choi.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.CL500H, plot.CK500H, plot.CGI500H, nrow = 2)
dev.off()

# Han's results
plot.HB500=
  ggplot() +
  geom_line (data = df.simHB, aes(Time, CBlood), col="blue", lwd=1)+
  geom_point(data = Obs.HB500, aes(Time, CBlood),size=2.5) + ylab("NPs concentration in blood (mg/L)") 

plot.HL500=
  ggplot() +
  geom_line (data = df.simHL, aes(Time, CLiver), col="blue", lwd=1)+
  geom_point(data = Obs.HL500, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.HK500=
  ggplot() +
  geom_line (data = df.simHK, aes(Time, CKidney), col="blue", lwd=1)+
  geom_point(data = Obs.HK500, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.HGI500=
  ggplot() +
  geom_line (data = df.simHGI, aes(Time, CGI), col="blue", lwd=1)+
  geom_point(data = Obs.HGI500, aes(Time, CGI),size=2.5) + ylab("NPs concentration in GI (mg/kg)") 

plot.HBr500=
  ggplot() +
  geom_line (data = df.simHBr, aes(Time, CBrain), col="blue", lwd=1)+
  geom_point(data = Obs.HBr500, aes(Time, CBrain),size=2.5) + ylab("NPs concentration in brain (mg/kg)") 

# Arrange plots into one figure and save as png file
png("220Fitting.500_Han.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.HB500, plot.HK500, plot.HL500, plot.HGI500, plot.HBr500, nrow = 3)
dev.off()


# Tao's results
plot.TL500=
  ggplot() +
  geom_line (data = df.simTL, aes(Time, CLiver), col="green", lwd=1)+
  geom_point(data = Obs.TL500, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 
# Arrange plots into one figure and save as png file
png("220Fitting.500_Tao.png", width = 6600, height = 6000, units = "px", res = 600)
plot.TL500
dev.off()


# Zhang's results
plot.ZS500=
  ggplot() +
  geom_line (data = df.simZS, aes(Time, CSpleen), col="purple", lwd=1)+
  geom_point(data = Obs.ZS500, aes(Time, CSpleen),size=2.5) + ylab("NPs concentration in spleen (mg/kg)") 

plot.ZL500=
  ggplot() +
  geom_line (data = df.simZL, aes(Time, CLiver), col="purple", lwd=1)+
  geom_point(data = Obs.ZL500, aes(Time, CLiver),size=2.5) + ylab("NPs concentration in liver (mg/kg)") 

plot.ZK500=
  ggplot() +
  geom_line (data = df.simZK, aes(Time, CKidney), col="purple", lwd=1)+
  geom_point(data = Obs.ZK500, aes(Time, CKidney),size=2.5) + ylab("NPs concentration in kidney (mg/kg)") 

plot.ZLu500=
  ggplot() +
  geom_line (data = df.simZLu, aes(Time, CLung), col="purple", lwd=1)+
  geom_point(data = Obs.ZLu500, aes(Time, CLung),size=2.5) + ylab("NPs concentration in lung (mg/kg)") 

# Arrange plots into one figure and save as png file
png("220Fitting.500_Zhang.png", width = 6600, height = 6000, units = "px", res = 600)
grid.arrange(plot.ZS500, plot.ZK500, plot.ZL500, plot.ZLu500, nrow = 2)
dev.off()
#----------------------------------------------------------------------------
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
mean.obsLB  <- c(17.8781)
mean.obsLL  <- c(77.22417777)
mean.obsLS  <- c(65.22446912)
mean.obsLK  <- c(55.26748181)
mean.obsLLu <- c(6.023075241)
mean.obsLBR <- c(14.37992126)
mean.obsLGI <- c(1775.361059)

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


## Model performance of Choi's results
# Specify the observed x values
specific_Cx_values <- c(14)

# Extract predicted y values for specific x values
selected_Crows <- df.simC[df.simC$Time %in% specific_Cx_values, ]
# Extracted y values
predicted_C1Lvalues  <- selected_Crows$CLiverL
predicted_C1Kvalues  <- selected_Crows$CKidneyL
predicted_C1GIvalues <- selected_Crows$CGIL
predicted_C2Lvalues  <- selected_Crows$CLiverM
predicted_C2Kvalues  <- selected_Crows$CKidneyM
predicted_C2GIvalues <- selected_Crows$CGIM
predicted_C3Lvalues  <- selected_Crows$CLiverH
predicted_C3Kvalues  <- selected_Crows$CKidneyH
predicted_C3GIvalues <- selected_Crows$CGIH
# Observed y values
mean.obsLL  <- c(0.274)
mean.obsLK  <- c(3.653)
mean.obsLGI <- c(4.064)
mean.obsML  <- c(1.05)
mean.obsMK  <- c(17.078)
mean.obsMGI <- c(13.196)
mean.obsHL  <- c(2.466)
mean.obsHK  <- c(28.95)
mean.obsHGI <- c(29.041)

# Calculate MAPE function
calculate_mape <- function(actual, predicted) {
  n <- length(actual)
  mape <- (1/n) * sum(abs((actual - predicted) / actual)) * 100
  return(mape)
}

mape_resultLL   <- calculate_mape(mean.obsLL, predicted_C1Lvalues)
mape_resultLK   <- calculate_mape(mean.obsLK, predicted_C1Kvalues)
mape_resultLGI  <- calculate_mape(mean.obsLGI, predicted_C1GIvalues)
mape_resultML   <- calculate_mape(mean.obsML, predicted_C2Lvalues)
mape_resultMK   <- calculate_mape(mean.obsMK, predicted_C2Kvalues)
mape_resultMGI  <- calculate_mape(mean.obsMGI, predicted_C2GIvalues)
mape_resultHL   <- calculate_mape(mean.obsHL, predicted_C3Lvalues)
mape_resultHK   <- calculate_mape(mean.obsHK, predicted_C3Kvalues)
mape_resultHGI  <- calculate_mape(mean.obsHGI, predicted_C3GIvalues)

# Print the result
cat("MAPE:", mape_resultLL, "%\n")
cat("MAPE:", mape_resultLK, "%\n")
cat("MAPE:", mape_resultLGI, "%\n")
cat("MAPE:", mape_resultML, "%\n")
cat("MAPE:", mape_resultMK, "%\n")
cat("MAPE:", mape_resultMGI, "%\n")
cat("MAPE:", mape_resultHL, "%\n")
cat("MAPE:", mape_resultHK, "%\n")
cat("MAPE:", mape_resultHGI, "%\n")

# Overall goodness of fit
observed.C  <- c(mean.obsLL, mean.obsLK, mean.obsLGI, mean.obsML, mean.obsMK, 
                 mean.obsMGI, mean.obsHL, mean.obsHK, mean.obsHGI)
predicted.C <- c(predicted_C1Lvalues, predicted_C1Kvalues, predicted_C1GIvalues, predicted_C2Lvalues, 
                 predicted_C2Kvalues, predicted_C2GIvalues, predicted_C3Lvalues, predicted_C3Kvalues,
                 predicted_C3GIvalues)

# Fit a linear regression model
lm_model.C <- lm(observed.C ~ predicted.C)

# Calculate R-squared
r_squared.C <- summary(lm_model.C)$r.squared

# Calculate adjusted R-squared
n <- length(observed.C)
p <- length(lm_model.C$coefficients) - 1  # Number of predictors excluding the intercept
adjusted_r_squared.C <- 1 - (1 - r_squared.C) * ((n - 1) / (n - p - 1))

# Calculate Root Mean Squared Error (RMSE)
rmse.C <- sqrt(mean((observed.C - predicted.C)^2))

# Print the results
cat("R-squared:", r_squared.C, "\n")
cat("Adjusted R-squared:", adjusted_r_squared.C, "\n")
cat("RMSE:", rmse.C, "\n")


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
mean.obsHB  <- c(1.05013)
mean.obsHL  <- c(9.814)
mean.obsHK  <- c(4.7713)
mean.obsHBR <- c(0.0533)
mean.obsHGI <- c(21.6607)

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


## Model performance of Tao's results
specific_Tx_values <- c(112)                     # Specify the observed x values
# Extract predicted y values for specific x values
selected_TLrows  <- df.simTL[df.simTL$Time %in% specific_Tx_values, ]
predicted_TLvalues  <- selected_TLrows$CLiver    # Extracted y values
mean.obsTL  <- c(1.37123e-7)                     # Observed y values
mape_resultTL   <- calculate_mape(mean.obsTL, predicted_TLvalues)  # Calculate MAPE 
cat("MAPE:", mape_resultTL, "%\n")               # Print the result

rmse.T <- sqrt(mean((mean.obsTL - predicted_TLvalues)^2))      # Calculate RMSE
# Print the results
cat("RMSE:", rmse.T, "\n")


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
mean.obsZL  <- c(482.92683)
mean.obsZS  <- c(569.40133)
mean.obsZK  <- c(715.425532)
mean.obsZLu <- c(144.717757)

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