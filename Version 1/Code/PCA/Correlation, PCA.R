## -Spleen--------------------------------------------2024.5.26
S <- read.csv('regression.S.csv')
str(S)   
head(S)

# Linear regression
Smodel <- lm(CSpleen ~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = S)
library(car) #first install.packages("car")
cor(S[, 2:6]) #correlation  #cor(x, method= c("pearson", "kendall", "spearman"))
predictors <- S[, 2:6]

require(GGally)
require(ggplot2)
lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_Spleen.png", width=20, height=18, units=c("cm"), dpi=400)

# K-Cluster method------------------------------------
library(factoextra)  #install.packages("factoextra") #factoextra 
library(NbClust)     #install.packages("NbClust")  # to estimate the number of clusters
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
             )  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
library(FactoMineR)  
require(factoextra)  
predictors.pca <- PCA(predictors, graph = FALSE) # FactoMineR package: PCA 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_Spleen.png", width=11, height=10, units=c("cm"), dpi=400)


## -Liver-----------------------------------------------------------------------
L <- read.csv('regression.L.csv')
str(L)   
head(L)
# Linear regression
Lmodel <- lm(CLiver~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = L)
cor(L[, 2:6])   # correlation  
predictors <- L[, 2:6]

lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_Liver.png", width=20, height=18, units=c("cm"), dpi=400)  

# K-Cluster method------------------------------------
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
             )  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
predictors.pca <- PCA(predictors, graph = FALSE) 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_Liver.png", width=11, height=10, units=c("cm"), dpi=400)


## -Blood-----------------------------------------------------------------------
B <- read.csv('regression.B.csv')
str(B)   
head(B)
# Linear regression
Bmodel <- lm(CBlood~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = B)
cor(B[, 2:6]) 
predictors <- B[, 2:6]

lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_Blood.png", width=20, height=18, units=c("cm"), dpi=400)  

# K-Cluster method------------------------------------
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
             )  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
predictors.pca <- PCA(predictors, graph = FALSE) 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_Blood.png", width=11, height=10, units=c("cm"), dpi=400)

## -GI-----------------------------------------------------------------------
GI <- read.csv('regression.GI.csv')
str(GI)   
head(GI)
# Linear regression
GImodel <- lm(CGI~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = GI)
cor(GI[, 2:6]) 
predictors <- GI[, 2:6]

lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_GI.png", width=20, height=18, units=c("cm"), dpi=400)  

# K-Cluster method------------------------------------
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
)  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
predictors.pca <- PCA(predictors, graph = FALSE) 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_GI.png", width=11, height=10, units=c("cm"), dpi=400)

## -Kidney-----------------------------------------------------------------------
K <- read.csv('regression.K.csv')
str(K)   
head(K)
# Linear regression
Kmodel <- lm(CKidney~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = K)
cor(K[, 2:6]) 
predictors <- K[, 2:6]

lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_Kidney.png", width=20, height=18, units=c("cm"), dpi=400)  

# K-Cluster method------------------------------------
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
)  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
predictors.pca <- PCA(predictors, graph = FALSE) 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_Kidney.png", width=11, height=10, units=c("cm"), dpi=400)


## -Lung-----------------------------------------------------------------------
Lu <- read.csv('regression.Lu.csv')
str(Lu)   
head(Lu)
# Linear regression
Lumodel <- lm(CLung~ KGIb + Kfeces + Particle.size + Exposure.time + Exposure.dose, data = Lu)
cor(Lu[, 2:6]) 
predictors <- Lu[, 2:6]

lower_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size = 1, color="black", alpha=0.6) +
    geom_smooth(method=loess, fill="red", color="red") +
    geom_smooth(method=lm, fill="blue", color="blue")
  p
}
g <- ggpairs(predictors, columns = 1:5, 
             lower = list(continuous = lower_fn), 
             diag = list(continuous = wrap("densityDiag", size = 1)),
             upper = list(continuous = wrap("cor", size = 4.5, fontface="bold", color="#FC4E07"))) +
  theme_minimal()+
  theme_bw(base_size = 15)
g+ theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text = element_text(size=10))

ggsave("scatterplot_Lung.png", width=20, height=18, units=c("cm"), dpi=400)  

# K-Cluster method------------------------------------
fviz_nbclust(predictors, FUNcluster = kmeans,    # K-Means
             method = "wss"      # total within sum of square
)  # Elbow method: 

kmeans.cluster <- kmeans(predictors, centers=2)   

# Variables factor map - Principal Component Analysis
predictors.pca <- PCA(predictors, graph = FALSE) 
summary(predictors.pca)  
predictors.pca.var <- get_pca_var(predictors.pca)  
predictors.pca.var$contrib      
fviz_pca_var(predictors.pca,    
             col.var = "contrib",       
             repel = TRUE,               
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
grp <- as.factor(kmeans.cluster$cluster)   
fviz_pca_biplot(predictors.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = grp, col.ind = "black",
                pointshape = 21, pointsize = 2.5,
                palette = "jco",
                # Variables
                col.var = "cos2",    # cos2 = the quality of the individuals on the factor map
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                legend.title = list(fill = "Cluster", color = "cos2"), repel = TRUE) 

ggsave("Variables factor map_Lung.png", width=11, height=10, units=c("cm"), dpi=400)
