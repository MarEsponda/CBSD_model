#Code 1
#Cassava disease dispersion model
#Matrix version loaded by blocks
#Maria del Mar Esponda 2024

library(sf)           # Spatial data manipulation
library(terra)        # Handling and analysis of raster data
library(raster)       # Library for raster data
library(geosphere)    # Geographic calculations
library(dplyr)        # Data manipulation
library(units)        # Handling units
library(bigmemory)    # Handling large datasets in memory

# Creation of the distance matrix -----------------------------------------------------------------------------------------------------------
setwd("your_file_location")                                                     # Set directory
raster <- raster("Study_zone1.tif")                                             # Base Raster (provinces that produce cassava)
raster_wgs84 <- projectRaster(raster,                                          
                              crs = CRS("+proj=longlat +datum=WGS84 +no_defs")) # Project raster to geographic coordinates
coords <- rasterToPoints(raster_wgs84)                                          # Generate one point per pixel
coords1 <- as.data.frame(coords)                                                # Convert to dataframe
coords1 <- coords1[,1:2]                                                        # Select only the first two columns
pix <- 140810                                                                   # Number of pixels in the raster
time_creation <- system.time({                                                  # Create an empty matrix
  distances <- big.matrix(nrow=pix, ncol=pix, type = "double", 
                          backingfile = "back.bin", descriptorfile = "descr.desc")
})
descriptor <- describe(distances)
dput(descriptor, file = "descr")
print(time_creation)                                                            # Elapsed: 7142.87 - 11528.86 s = 2 - 3.2 hours
object.size(distances)
# user: The CPU time used by the user to execute the code in seconds.
# system: The CPU time used by the system (operating system) to execute the code in seconds.
# elapsed: The total elapsed time from the start to the end of the code execution in seconds.
# Measure execution time with system.time()
dist_matrix <- attach.big.matrix("descr.desc")
time_filling <- system.time({                                                   
  for (i in 1:pix){
    dist_matrix[i,] <- distm(coords1[i,], coords1, fun = distHaversine)
    print(paste("Currently processing cycle", i))
  }
})
print(time_filling)                                                             # Elapsed: 7106.92 s - 14062.83 - 46474.05 = 1.97 hours - 3.9 h - 12.9 h
object.size(dist_matrix)

dist_matrix <- attach.big.matrix("descr.desc")

# Matriz en km--------------------------------------------------------------------------------------------------- 
backingfile <- "distanciaskm.bin"
descriptorfile <- "distanciaskm.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
filebacked_block[,] <- dist[,]
distkm <- attach.big.matrix("distanciaskm.desc")
for(i in 1:pix){
  distkm[i,] <- round((dist_matrix[i,]/1000),3)                                             #Dividimos por 1000
  print(paste0("lapso",i))
}

# Matriz de entrada para la ley de potencia inversa k=2---------------------------------------------------------------
backingfile <- "distpi2.bin"
descriptorfile <- "distpi2.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
for(i in 1:pix){
  distpi2[i,] <- distkm[i,]
  print(paste0("lapso",i))
}
distpi2 <- attach.big.matrix("distpi2.desc")
for (i in 1:pix) {                                                           #LLenar mi diagonal de 1 (evitar la division por 0)
  distpi2[i, i] <- 1
  print(paste0("fila",i))
}
k<-2
for(i in 1:pix){
  distpi2[i,] <- round((1/(distpi2[i,]^k)),4)                                   #Matriz de probabilidad de distancias por potencia inversa, luego creamos una matriz con los pesos con el k optimizado 
  print(paste0("lapso",i))
}
#Esta matriz se guarda automaticamente
#La matriz de entrada para la exponencial negativa no tiene restricciones, es entonces la misma de "distkm.desc"
#Una vez se hayan guardado estas matrices, es necesario almacenarlas en otra carpeta para tener su Respaldo

# Matriz de entrada para la ley de potencia inversa k=2.5---------------------------------------------------------------
distkm <- attach.big.matrix("distanciaskm.desc")
backingfile <- "distpi2.5.bin"
descriptorfile <- "distpi2.5.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
distpi2.5 <- attach.big.matrix("distpi2.5.desc")
for(i in 1:pix){
  distpi2.5[i,] <- distkm[i,]
  print(paste0("lapso",i))
}
for (i in 1:pix) {                                                           #LLenar mi diagonal de 1 (evitar la division por 0)
  distpi2.5[i, i] <- 1
  print(paste0("fila",i))
}
k<-2.5
for(i in 1:pix){
  distpi2.5[i,] <- round((1/(distpi2.5[i,]^k)),4)                               #Matriz de probabilidad de distancias por potencia inversa, luego creamos una matriz con los pesos con el k optimizado 
  print(paste0("lapso",i))
}

for(i in 1:pix){
  condicion <- distkm[i,] > 350                                                 #Todos los que tengan la distancia mayor a 350 km
  indices <- which(condicion, arr.ind = TRUE)                                   #Listamos las posiciones donde se cumple la condicion
  distpi2.5[i,indices] <- 0                                                     #Asignar 0 a las distancias mayores a 350 km
  print(paste0("lapso",i))
}

#Esta matriz se guarda automaticamente
#La matriz de entrada para la exponencial negativa no tiene restricciones, es entonces la misma de "distkm.desc"
#Una vez se hayan guardado estas matrices, es necesario almacenarlas en otra carpeta para tener su Respaldo

# Matriz de entrada para la ley de potencia inversa k=3---------------------------------------------------------------
backingfile <- "distpi3.bin"
descriptorfile <- "distpi3.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
for(i in 1:pix){
  distpi3[i,] <- distkm[i,]
  distpi3[i, i] <- 1
  print(paste0("lapso",i))
}
distpi3 <- attach.big.matrix("distpi3.desc")
for (i in 1:pix) {                                                           #LLenar mi diagonal de 1 (evitar la division por 0)
  distpi3[i, i] <- 1
  print(paste0("fila",i))
}
k<-3
for(i in 1:pix){
  distpi3[i,] <- round((1/(distpi3[i,]^k)),4)                                         #Matriz de probabilidad de distancias por potencia inversa, luego creamos una matriz con los pesos con el k optimizado 
  print(paste0("lapso",i))
}
#Esta matriz se guarda automaticamente
#La matriz de entrada para la exponencial negativa no tiene restricciones, es entonces la misma de "distkm.desc"
#Una vez se hayan guardado estas matrices, es necesario almacenarlas en otra carpeta para tener su Respaldo


# Matriz de entrada para la ley de potencia inversa k=3.5 (mas cercano)---------------------------------------------------------------
backingfile <- "distpi3.5.bin"
descriptorfile <- "distpi3.5.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
for(i in 1:pix){
  distpi3.5[i,] <- distkm[i,]
  distpi3.5[i, i] <- 1
  print(paste0("lapso",i))
}
distpi3.5 <- attach.big.matrix("distpi3.5.desc")
for (i in 1:pix) {                                                           #LLenar mi diagonal de 1 (evitar la division por 0)
  distpi3.5[i, i] <- 1
  print(paste0("fila",i))
}
k<-3.5
for(i in 1:pix){
  distpi3.5[i,] <- round((1/(distpi3.5[i,]^k)),4)                                         #Matriz de probabilidad de distancias por potencia inversa, luego creamos una matriz con los pesos con el k optimizado 
  print(paste0("lapso",i))
}
#Esta matriz se guarda automaticamente
#La matriz de entrada para la exponencial negativa no tiene restricciones, es entonces la misma de "distkm.desc"
#Una vez se hayan guardado estas matrices, es necesario almacenarlas en otra carpeta para tener su Respaldo


# Matriz de entrada para la ley de potencia inversa k=4---------------------------------------------------------------
backingfile <- "distpi4.bin"
descriptorfile <- "distpi4.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
for(i in 1:pix){
  distpi4[i,] <- distkm[i,]
  distpi4[i, i] <- 1
  print(paste0("lapso",i))
}
distpi4 <- attach.big.matrix("distpi4.desc")
k<-4
for(i in 1:pix){
  distpi4[i,] <- round((1/(distpi4[i,]^k)),4)                                         #Matriz de probabilidad de distancias por potencia inversa, luego creamos una matriz con los pesos con el k optimizado 
  print(paste0("lapso",i))
}
#Esta matriz se guarda automaticamente
#La matriz de entrada para la exponencial negativa no tiene restricciones, es entonces la misma de "distkm.desc"
#Una vez se hayan guardado estas matrices, es necesario almacenarlas en otra carpeta para tener su Respaldo



# Matriz de entrada para propagacion por mosca blanca---------------------------
backingfile <- "distmosca.bin"
descriptorfile <- "distmosca.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = pix, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
for(i in 1:pix){
  filebacked_block[i,] <- distkm[i,]
  print(paste0("lapso",i))
}

distkmosca <- attach.big.matrix("distmosca.desc")                               #Distancias menores a 5 km
for(i in 1:pix){
  condicion <- distkmosca[i,] > 5                                               #Todos los que tengan la distancia mayor a 5 km (considerando pixeles diagonales)
  indices <- which(condicion, arr.ind = TRUE)                                   #Listamos las posiciones donde se cumple la condicion
  distkmosca[i,indices] <- 0                                                    #Asignar 0 a las distancias mayores a 5km
  print(paste0("lapso",i))
}

for(i in 1:pix){
  condicion <- distkmosca[i,] == 0                                               #Todos los que tengan la distancia mayor a 5 km (considerando pixeles diagonales)
  indices <- which(condicion, arr.ind = TRUE)                                   #Listamos las posiciones donde se cumple la condicion
  distkmosca[i,indices] <- NA                                                    #Asignar 0 a las distancias mayores a 5km
  print(paste0("lapso",i))
}

# Matriz de within_distance-------------------------------------------------
backingfile <- "within_distance.bin"
descriptorfile <- "within_distance.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = 1, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
within_distance <- attach.big.matrix("within_distance.desc")

# Matriz de prob_infection--------------------------------------------------
backingfile <- "prob_infection.bin"
descriptorfile <- "prob_infection.desc"
pix<-140810
filebacked_block <- big.matrix(nrow = pix, ncol = 1, type = "double", backingfile = backingfile, descriptorfile = descriptorfile)
prob_infection <- attach.big.matrix("prob_infection.desc")