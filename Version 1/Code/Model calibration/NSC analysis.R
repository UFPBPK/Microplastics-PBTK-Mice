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
mod <- mcode ("PBTK", PBTK.code) 
## Normalized sensitivity analysis---------------------------------------------
# Prediction model w/ AUC output
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
    mod %>%                                             # model object
    param(pars) %>%                                     # to update the parameters in the model subject
    #Req(AUC_Spleen, AUC_Kidney, AUC_Blood, AUC_Lung, AUC_Liver, AUC_GI, AUC_Brain, AUC_Rest)%>%                             # select model output
    update(atol = 1e-70, maxstep = 50000) %>%
    mrgsim_d(data = ex.K, tgrid = tsamp.K)
  out.K <-cbind.data.frame(Time    = out.K$time/24,     # day
                           AUC_Blood   = out.K$AUC_B,
                           AUC_GI      = out.K$AUC_GI,
                           AUC_Spleen  = out.K$AUC_S,
                           AUC_Kidney  = out.K$AUC_K,
                           AUC_Liver   = out.K$AUC_L,
                           AUC_Lung    = out.K$AUC_Lu,
                           AUC_Brain   = out.K$AUC_BR,
                           AUC_Rest    = out.K$AUC_Rt
  )
  return(out.K)
}

## Define the sensitivity function
# Normalized Sensitivity Coefficient (NSC)
NSC_func <- function (pars) {
  n <- length(pars)
  NSC = matrix(NA, nrow = n , ncol = 8)
  
  for (i in 1:n) {
    New.pars      <- pars %>% replace(i, log(exp((pars[i]))*1.01))  # Each parameter was increased by 1%
    Rnew          <- pred(New.pars)
    R             <- pred(pars)
    delta.pars    <- exp(pars[i])/(exp(pars[i])*0.01)
    
    ## Estimated the AUC
    AUC.new  <- Rnew %>% filter(Time == 1)%>% select(-contains("Time"))  #24 h = 1 day
    AUC.ori  <- R    %>% filter(Time == 1)%>% select(-contains("Time"))
    
    for (j in 1:dim(AUC.new)[2]) {
      delta.AUC    =  AUC.new [,j]- AUC.ori[,j]
      NSC[i, j]   <- (delta.AUC/AUC.ori[,j]) * delta.pars
    }
  }
  rownames(NSC)= c(names(pars))
  colnames(NSC)= c("AUCB_24","AUCGI_24","AUCS_24","AUCK_24","AUCL_24","AUCLu_24","AUCBR_24","AUCRt_24")
  return (NSC = NSC)
}

## Calculate NSC------------------------------------------------------------------
# 6 micron
theta.final.6 <- log(c(
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
  ASREScap      = 200,      # Uptake capacity per tissue weight 
  ALuREScap     = 15,
  AKREScap      = 15,
  ALREScap      = 100,
  KD            = 4468.93,  # Fitted
  N             = 1,
  KGIb          = 4.01e-5,  # Fitted  # Adjusted (decrease) #Absorption rate of GI tract (1/h)
  Kfeces        = 0.508,    # Fitted  # Fecal clearance (L/h)
  KbileC        = 0.0012,   #Biliary clearance (L/h/kg^0.75)
  KurineC       = 0.000429  # Fitted  # Urinary clearance (L/h/kg^0.75)
))

NSC.6<- NSC_func(theta.final.6)
write.csv(NSC.6, file = "NCS_6.csv")

# 1 micron
theta.final.1 <- log(c(
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
  ASREScap      = 200,      # Uptake capacity per tissue weight 
  ALuREScap     = 15,
  AKREScap      = 15,
  ALREScap      = 100,
  KD            = 4092.26,   # Fitted 
  N             = 1,
  KGIb          = 3.42e-5,   # Fitted   # Absorption rate of GI tract (1/h)
  Kfeces        = 0.5046,    # Fitted   # Fecal clearance (L/h)
  KbileC        = 0.0012,   #Biliary clearance (L/h/kg^0.75)
  KurineC       = 0.000564   # Fitted   # Urinary clearance (L/h/kg^0.75)
))

NSC.1<- NSC_func(theta.final.1)
write.csv(NSC.1, file = "NCS_1.csv")

# 220 nm
theta.final.220 <- log(c(
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
  ASREScap      = 200,      # Uptake capacity per tissue weight 
  ALuREScap     = 15,
  AKREScap      = 15,
  ALREScap      = 100,
  KD            = 2464.7,    # Fitted  
  N             = 1,
  KGIb          = 4.22e-5,   # Fitted  # Absorption rate of GI tract (1/h)
  Kfeces        = 0.55,      # Fitted  # Fecal clearance (L/h)
  KbileC        = 0.0012,   #Biliary clearance (L/h/kg^0.75)
  KurineC       = 0.0003     # Fitted  # Urinary clearance (L/h/kg^0.75)
))

NSC.220<- NSC_func(theta.final.220)
write.csv(NSC.220, file = "NCS_220.csv")

# 20 nm
theta.final.20 <- log(c(
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
  ASREScap      = 200,      # Uptake capacity per tissue weight 
  ALuREScap     = 15,
  AKREScap      = 15,
  ALREScap      = 100,
  KD            = 558.04,   # Fitted
  N             = 1,
  KGIb          = 4.23e-5,  # Fitted  # Absorption rate of GI tract (1/h)
  Kfeces        = 0.548,    # Fitted  # Fecal clearance (L/h)
  KbileC        = 0.0012,   # Biliary clearance (L/h/kg^0.75)
  KurineC       = 0.00012   # Urinary clearance (L/h/kg^0.75)
))

NSC.20<- NSC_func(theta.final.20)
write.csv(NSC.20, file = "NCS_20.csv")


## Plot for 6 micron
# Define the breaks 
# Compute percentiles
percentiles <- quantile(NSC.6, probs = seq(0, 1, 0.1))

# Print the percentiles
print(percentiles)
breaks <- c(-1, -0.5, -0.2, -0.05, -0.01, -0.005, -5e-4, -1e-5, -1e-8, 
            0, 1e-8, 1e-5, 5e-4, 0.005, 0.01, 0.05, 0.2, 0.5, 1)

NSC_plot <- pheatmap(as.matrix(NSC.6), 
                     cutree_rows = 4,          # hclust: Hierarchical Clustering
                     clustering_distance_rows = "correlation",   # distance measure
                     cluster_cols = F,
                     fontsize = 5.5,
                     color = magma(length(breaks) - 1),                      
                     breaks = breaks,     
                     border_color = "grey20",
                     legend = T,
                     cellwidth = 12, cellheight = 5,
                     main="NSC")

# Re-order original data to match ordering in heatmap (top-to-bottom)
rownames(NSC.6[NSC_plot$tree_row[["order"]],])


# Sort parameters into 4 groups  
row_P_grp <- sort(cutree(NSC_plot$tree_row, k=4))    
row_P <- data.frame(Cluster = factor(row_P_grp)) # Rename and create a data frame
levels(row_P$Cluster) <- paste("Cluster", 1:4, sep = "")


# Create row_col list with specified colors
row_col <- list(Cluster = c(Cluster1 = "mediumpurple1", Cluster2 = "lavender", 
                            Cluster3 = "khaki1", Cluster4 = "gold2"))

# Plot color swatches for heatmap
color_scale <- colorRampPalette(brewer.pal(5, "RdBu"))(length(breaks) - 1)
coul<- alpha(color_scale, alpha = 0.9)  # Adjust transparency of colors

par(mar = c(0, 0, 1, 0))
image(1:length(breaks), 1, as.matrix(1:length(breaks)), col = coul,
      xlab = "Predefined breaks", ylab = "", axes = FALSE)

# Add annotation for row cluster
NSC_plot2 <- pheatmap(as.matrix(NSC.6), cutree_rows = 4,
                      clustering_distance_rows = "correlation",  
                      color = coul,
                      breaks = breaks,      # Specify the logarithmic breaks
                      annotation_row = row_P,
                      annotation_colors = row_col,
                      cluster_cols = FALSE,
                      fontsize = 5.5,
                      legend = TRUE,
                      cellwidth = 15, cellheight = 5,
                      border_color = "grey30",
                      main="NSC-6 µm")

png("NSC.cluster_6.png", width = 7500, height = 7500, units = "px", res = 1800)
NSC_plot2
dev.off()

## Plot for 1 micron
NSC_plot <- pheatmap(as.matrix(NSC.1), cutree_rows = 4, 
                     clustering_distance_rows = "correlation",
                     cluster_cols = F,
                     fontsize = 5.5,
                     color = magma(length(breaks) - 1), 
                     breaks = breaks,
                     border_color = "grey20",
                     legend = T,
                     cellwidth = 12, cellheight = 5,
                     main="NSC")

##-----ignore this part!! Use the same cluster from 6 µm to observe the changes----
# Re-order original data to match ordering in heatmap (top-to-bottom)
rownames(NSC.1[NSC_plot$tree_row[["order"]],])
# Sort parameters into 4 groups
row_P_grp <- sort(cutree(NSC_plot$tree_row, k=4))
row_P <- data.frame(Cluster = factor(row_P_grp)) # Rename and create a data frame
levels(row_P$Cluster) <- paste("Cluster", 1:4, sep = "")


# Create row_col list with specified colors
row_col <- list(Cluster = c(Cluster1 = "violetred", Cluster2 = "cyan3", 
                            Cluster3 = "violet", Cluster4 = "olivedrab1"))
##-----ignore this part!! Use the same cluster from 6 µm to observe the changes----

# Add annotation for row cluster
NSC_plot2 <- pheatmap(as.matrix(NSC.1), cutree_rows = 4,
                      clustering_distance_rows = "correlation",
                      color = coul, 
                      breaks = breaks,
                      annotation_row = row_P,
                      annotation_colors = row_col,
                      cluster_cols = FALSE,
                      fontsize = 5.5,
                      legend = TRUE,
                      cellwidth = 15, cellheight = 5,
                      border_color = "grey30",
                      main="NSC-1 µm")

png("NSC.cluster_1.png", width = 7500, height = 7500, units = "px", res = 1800)
NSC_plot2
dev.off()

## Plot for 220 nm
NSC_plot <- pheatmap(as.matrix(NSC.220), cutree_rows = 4, 
                     clustering_distance_rows = "correlation",
                     cluster_cols = F,
                     fontsize = 5.5,
                     color  = magma(length(breaks) - 1), 
                     breaks = breaks,
                     border_color = "grey20",
                     legend = T,
                     cellwidth = 12, cellheight = 5,
                     main="NSC")

# Add annotation for row cluster
NSC_plot2 <- pheatmap(as.matrix(NSC.220), cutree_rows = 4, 
                      clustering_distance_rows = "correlation",
                      color = coul, 
                      breaks = breaks,
                      annotation_row = row_P,
                      annotation_colors = row_col,
                      cluster_cols = FALSE,
                      fontsize = 5.5,
                      legend = TRUE,
                      cellwidth = 15, cellheight = 5,
                      border_color = "grey30",
                      main="NSC-220 nm")

png("NSC.cluster_220.png", width = 7500, height = 7500, units = "px", res = 1800)
NSC_plot2
dev.off()

## Plot for 20 nm
NSC_plot <- pheatmap(as.matrix(NSC.20), cutree_rows = 4, 
                     clustering_distance_rows = "correlation",
                     cluster_cols = F,
                     fontsize = 5.5,
                     color = magma(length(breaks) - 1), 
                     breaks = breaks,
                     border_color = "grey20",
                     legend = T,
                     cellwidth = 12, cellheight = 5,
                     main="NSC")

# Add annotation for row cluster
NSC_plot2 <- pheatmap(as.matrix(NSC.20), cutree_rows = 4,
                      clustering_distance_rows = "correlation",
                      color = coul, 
                      breaks = breaks,
                      annotation_row = row_P,
                      annotation_colors = row_col,
                      cluster_cols = FALSE,
                      fontsize = 5.5,
                      legend = TRUE,
                      cellwidth = 15, cellheight = 5,
                      border_color = "grey30",
                      main="NSC-20 nm")

png("NSC.cluster_20.png", width = 7500, height = 7500, units = "px", res = 1800)
NSC_plot2
dev.off()
