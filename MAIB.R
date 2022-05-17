#MAIB - Methodology for Assessing the Uncertainty of Bathymetric Data
###########################################################################


#Data path
setwd("D:/LAURA/Manuscript/MAIB")

getwd() 

#List of packages to be used in the Script
pkg <- c("geoR","moments","scatterplot3d","tcltk2",
         "sp", "rgdal", "ggplot2" , "cluster" ,
         "bootstrap", "plyr", "robustbase", "MBESS", "rgeos", "gstat")

sapply(pkg, require, character.only=TRUE)



#Data reading (without outliers)
dados <- read.table("manuscript_data.txt", header=T, dec=",")
names(dados)
dados
length(dados$dz)

######################################################################
#Loading Functions
######################################################################

#Function to calculate uncertainty
theta <- function(x){(sqrt((sd(x)^2)+(mean(x)^2)))}

#Function to calculate RMSE (divided by n)
theta1 <- function(x){(sqrt((sum(x^2))/length(x)))}

#Function to calculate TCL Uncertainty
theta2 <- function(x,y){(sqrt(((sd(x)^2)*y)+(mean(x)^2)))}

#Function to calculate robust uncertainty
theta3 <- function(x){(sqrt((mad(x)^2)+(median(x)^2)))} 

######################################################################
######################### Amostra dependente #########################
######################################################################

#####################################################################
#Block Bootstrap: Generate IC with 95%
#####################################################################

#Obtaining the data
#Load data without outliers
dados <- read.table("xyzestdz_0106.txt", header=T, dec=",")
coordinates(dados) <- c("X", "Y")

#Calculation of statistics
ivt3 = theta(dados$dz)
rms1 = theta1(dados$dz)
ivt_robs = theta3(dados$dz)

#Generate blocks with a predefined diagonal size
#suggestion: Range obtained from Geostatistical analysis
tamanho <- 80

# data.frame
#data(dados$dz)
#coordinates(dados$dz) <- ~x+y
#gridded(dados$dz) <- TRUE
#bbox(dados$dz)

#data(dados)
#coordinates(dados) <- ~x+y
#proj4string(dados) <- CRS("+init=epsg:28992")

#Delimit the number and location of each block
Bloco <- makegrid(bbox(dados), cellsize = (tamanho*sqrt(2)), pretty = FALSE)
coordinates(Bloco) <- c("x1","x2")
gridded(Bloco) <- TRUE
Bloco <- as.SpatialPolygons.GridTopology(Bloco@grid)
plot(Bloco)     #Plotando os Blocos
points(dados) #Plotando os dados originais

#Extract the Block number where each point is overlapping
ptsInBloco <- as.numeric(gIntersects(dados, Bloco, byid=TRUE, 
                                     returnDense=FALSE, checkValidity=TRUE))  #Todos que tem intersecao para cada buffer

#Number of Bootstrap replications
n_vezes <- 1000  

tab_boot <- tab_boot_dz <- NULL 

for (i in 1:n_vezes)   
{
  #FIRST - Block draw
  Grid <- sample(unique(ptsInBloco),dim(dados)[1], replace = TRUE)
  #SECOND - drawing a point within each Block selected previously
  pontos <- (as.numeric(lapply(Grid,function(x) sample(which(ptsInBloco==x),1))))  #Os pontos repetidos sao contabilizados apenas uma vez
  tab_boot <- rbind(tab_boot,(pontos))
  tab_boot_dz <- rbind(tab_boot_dz,dados@data$dz[pontos])
}

#Convert data to data.frame
tab_boot <- as.data.frame(tab_boot)
tab_boot_dz <- as.data.frame(tab_boot_dz)


#Generate a new dataset after the block bootstrap 
dados_novos_ivt <- apply(tab_boot_dz,1,theta) 
dados_novos_rms <- apply(tab_boot_dz,1,theta1) 
dados_novos_robs <- apply(tab_boot_dz,1,theta3) 

#Block bootstrap confidence intervals
IC_ivt <- quantile(dados_novos_ivt,c(0.025, 0.975)) 
IC_rms <- quantile(dados_novos_rms,c(0.025, 0.975)) 
IC_robs <- quantile(dados_novos_robs,c(0.025, 0.975)) 

#bias
vies_ivt <- ivt3 - median(dados_novos_ivt)
vies_rms <- rms1 - median(dados_novos_rms)
vies_robs <- ivt_robs - median(dados_novos_robs)

windows(8,8,title="Gráficos para análise exploratória")
par(mfrow=c(3,2), family="serif")
hist(dados_novos_ivt, xlab="Incerteza (m)", ylab= "Frequência", main=" Histograma (bootstrap)")
qqnorm(dados_novos_ivt, xlab="Quantis Teóricos", ylab= "Quantis Amostrados", main=" Normal Q-Q Plot (bootstrap)")
qqline(dados_novos_ivt,lty=2, col='red')
hist(dados_novos_rms, xlab="RMSE (m)", ylab= "Frequência", main=" Histograma (bootstrap)")
qqnorm(dados_novos_rms, xlab="Quantis Teóricos", ylab= "Quantis Amostrados", main=" Normal Q-Q Plot (bootstrap)")
qqline(dados_novos_rms,lty=2, col='red')
hist(dados_novos_robs, xlab="Incerteza Robusta (m)", ylab= "Frequência", main=" Histograma (bootstrap)")
qqnorm(dados_novos_robs, xlab="Quantis Teóricos", ylab= "Quantis Amostrados", main=" Normal Q-Q Plot (bootstrap)")
qqline(dados_novos_robs,lty=2, col='red')
par(mfrow=c(1,1), family="serif")

#Exporting information:
sink("Resultados.txt", type="output", append=T)

cat(" Incerteza Vertical \n Amostra Dependente - Bloco Bootstrap","\n",
    "N?mero de replica??es: " ,n_vezes,"\n",
    "Tamanho do lado do Bloco (m): " ,round (tamanho*sqrt(2),3),"\n",
    "\n Intervalo de Confian?a de 95%","\n",
    "------------------------------------------------------","\n",
    "Incerteza (m): "         ,round(ivt3,3)   ,"\n",
    "IC (m): "         ,"[",round(IC_ivt[1],3),";",round(IC_ivt[2],3),"]","\n",
    "Vi?s Bootstrap (m): "         ,round(vies_ivt,3)   ,"\n",
    
    "\n RMSE (m): "         ,round(rms1,3)   ,"\n",
    "IC (m): "         ,"[",round(IC_rms[1],3),";",round(IC_rms[2],3),"]","\n",
    "Vi?s Bootstrap (m): "         ,round(vies_rms,3)   ,"\n",
    
    "\n Incerteza Robusta (m): "         ,round(ivt_robs,3)   ,"\n",
    "IC (m): "         ,"[",round(IC_robs[1],3),";",round(IC_robs[2],3),"]","\n",
    "Vi?s Bootstrap (m): "         ,round(vies_robs,3)   ,"\n",
    "------------------------------------------------------","\n",
    fill=F)
sink()
shell.exec("Resultados.txt")