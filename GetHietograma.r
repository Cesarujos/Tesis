
library(chirps)
library(moments)
library(dplyr)

library(sf)
library(sp)
library(rgdal)
library(raster)
library(gstat)

library(openxlsx)
###Insercción de ruta de trabajo
setwd("C:\\Users\\uriar\\OneDrive\\UNMSM\\9no Ciclo\\Elaboración de tesis\\Tesis")

######DATOS DE ENTRADA############Subcuenca exportada de HEC-GeoHMS
subcuencas <- st_read("C:\\Users\\uriar\\OneDrive\\UNMSM\\Tesis\\Files\\Subcuencas.shp")


################PROCESO#####################
#Centroides

centroides<- SpatialPoints(coords=coordinates(as_Spatial(subcuencas)),
                            proj4string=CRS(proj4string(as_Spatial(subcuencas))))
centroides$DrainID <- subcuencas$DrainID

centroides_wgs84 <- spTransform(centroides, CRS("+init=epsg:4326"))
coord <- coordinates(centroides_wgs84)
coord <- as.data.frame(coord)
coord <- cbind(subcuenca = paste("W",centroides$DrainID,"0", sep = ""), coord)
#Años
anios <- list(inicio = 1982, final = 2021) #Año de inicio y final

########PROCESO##########
#PRUEBA DE DATOS DUDOSOS - outlier
dudoso <- function(ppmax){
  Kn <- 2.682 ##Para 40 datos
  ppmax_log <- c()
  for(item in ppmax){
    ppmax_log <- c(ppmax_log, log10(item))
  }
  promedio<- mean(ppmax_log)
  desviacion<- sd(ppmax_log)
  xhs <- promedio + Kn*desviacion
  xhi <- promedio - Kn*desviacion
  top_pp <- 10^xhs
  down_pp <- 10^xhi
  for(item in ppmax){
    if (top_pp < item){
      item <- -99
    }
    if(item < down_pp){
      item <- -99
    }
  }
  return(ppmax)
}

##CALCULO DE ASIMETRÍA
asimetria <- function(cadena){
  cantidad <- as.numeric(length(cadena))
  promtemp <- mean(cadena)
  desv_temp <- sd(cadena)
  temp <- 0
  for (item in cadena) {
    x <- (item - promtemp)^3
    temp <- temp + x
  }
  asimetria <- cantidad*temp/((cantidad-1)*(cantidad-2)*desv_temp^3)
  return(asimetria)
}

#ELABORACIÓN DE TABLA
get_p24 <- function(){
  lonlat <- data.frame(lon = coord$coords.x1,
                       lat = coord$coords.x2) #Dataframe de latitud y longitud
  dates <- c(paste(anios$inicio,"01-01",sep = "-"), paste(anios$final,"12-31",sep = "-"))
  dat <- get_chirps(lonlat, dates, server = "ClimateSERV")
  result_dataframe <- data.frame(anio = c(anios$inicio:anios$final))
  for( i in c(1:nrow(lonlat))){
    punto <- dat %>%
      filter(id == i)
    ppmax <- numeric()
    for(a in c(anios$inicio:anios$final)){
      datos <- punto %>%
        filter(as.numeric(format(date,'%Y')) == a) %>%
        summarise(max(chirps))
      datos <- datos$`max(chirps)`
      ppmax <- c(ppmax, datos)
    }
    ppmax <- dudoso(ppmax)
    result_dataframe = cbind(result_dataframe, i = ppmax)
  }
  colnames(result_dataframe) <- c('anio', coord$subcuenca )
  return(result_dataframe)
}

p24 <- get_p24()

#write.csv(p24, file="Prueba.csv")
tr = c(2, 5, 10, 20, 50, 100, 200, 500, 1000)

###Distribución NORMAL
normal <- function(cadena){
  prom <- mean(cadena)
  desv <- sd(cadena)
  prep_normal<- c()
  for (tr_i in tr) {
    a <- qnorm((1-(1/(tr_i))), prom, desv)
    prep_normal <- c(prep_normal,a)
  }
  cadena <- sort(cadena, decreasing = FALSE)
  tr_cadena <- c()
  for (i in c(1:length(cadena))) {
    n <- i/(length(cadena)+1)
    t <- 1/(1-n)
    tr_cadena <- c(tr_cadena, t)
  }
  tr_cadena <- qnorm((1-(1/(tr_cadena))), prom, desv)
  err <- (sum((cadena - tr_cadena)^2)/length(cadena))^0.5
  return(list(data = prep_normal, error = err))
}

###Distribución LogNormal###
log_normal <- function(cadena){
  cadena_log <- log(cadena)
  prom <- mean(cadena_log)
  desv <- sd(cadena_log)
  
  cadena <- sort(cadena, decreasing = FALSE)
  prep_logNormal <-c()
  for (tr_i in tr) {
    a <- qlnorm((1-(1/tr_i)),prom,desv, lower.tail = TRUE)
    prep_logNormal <- c(prep_logNormal,a)
  }
  
  tr_cadena <- c()
  for (i in c(1:length(cadena))) {
    n <- i/(length(cadena)+1)
    t <- 1/(1-n)
    tr_cadena <- c(tr_cadena, t)
  }
  tr_cadena <- qlnorm((1-(1/(tr_cadena))), prom, desv)
  err <- (sum((cadena - tr_cadena)^2)/length(cadena))^0.5
  
  return(list(data = prep_logNormal, error = err))
}


###Pearson tipo III ########
pearsonIII <- function(cadena){
  g <- asimetria(cadena)
  prom <- mean(cadena)
  desv <- sd(cadena)
  
  pearson <- c()
  for (tr_i in tr) {
    kt <- 2/g*((g/6*(qnorm(1-1/tr_i)-g/6)+1)^3-1) 
    pearson <- c(pearson, (prom + desv*kt))
  }
  #Calculo Error
  tr_cadena <- c()
  for (i in c(1:length(cadena))) {
    n <- i/(length(cadena)+1)
    t <- 1/(1-n)
    tr_cadena <- c(tr_cadena, t)
  }
  cadena <- sort(cadena, decreasing = FALSE)
  tr_cadena <- (prom + desv*(2/g*((g/6*(qnorm(1-1/tr_cadena)-g/6)+1)^3-1)))
  err <- (sum((cadena - tr_cadena)^2)/length(cadena))^0.5
  return(list(data = pearson, error = err))
}
####Log pearson#####
log_personIII <- function(cadena){
  cadena_log <- log(cadena)
  prom <- mean(cadena_log)
  desv <- sd(cadena_log)
  g <- asimetria(cadena_log)
  
  ####Extrapolación
  log_pearson <- c()
  for (tr_i in tr) {
    kt <- 2/g*((g/6*(qnorm(1-1/tr_i)-g/6)+1)^3-1) 
    log_pearson <- c(log_pearson, exp(prom + desv*kt))
  }
  
  #Error
  tr_cadena <- c()
  for (i in c(1:length(cadena))) {
    n <- i/(length(cadena)+1)
    t <- 1/(1-n)
    tr_cadena <- c(tr_cadena, t)
  }
  cadena <- sort(cadena, decreasing = FALSE)
  tr_cadena <- exp(prom + desv*(2/g*((g/6*(qnorm(1-1/tr_cadena)-g/6)+1)^3-1)))
  err <- (sum((cadena - tr_cadena)^2)/length(cadena))^0.5
  
  return(list(data = log_pearson, error = err))
}

#####Gumbell####

gumbell<- function(cadena){
  #Parametros necesarios
  prom <- mean(cadena)
  desv <- sd(cadena)
  g <- asimetria(cadena)
  alfa <- 0.7797*desv
  beta <- prom - 0.45*desv
  
  #Extrapolation
  gumbell_tr <- c()
  for (tr_i in tr) {
    y <- beta-alfa*log(-log(1-1/tr_i))
    gumbell_tr <- c(gumbell_tr, y)
  }
  
  #Error
  cadena <- sort(cadena, decreasing = FALSE)
  tr_cadena <- c()
  for (i in c(1:length(cadena))) {
    n <- i/(length(cadena)+1)
    t <- 1/(1-n)
    tr_cadena <- c(tr_cadena, t)
  }
  tr_cadena <- beta-alfa*log(-log(1-1/tr_cadena))
  err <- (sum((cadena - tr_cadena)^2)/length(cadena))^0.5
  
  return(list(data = gumbell_tr, error = err))
}

##Obtener extrapolacion de datos

get_extrapolation <- function(dataframe){
  extrapolation <- data.frame( retorno = tr)
  err_punt <- c()
  for(item in c(2:ncol(dataframe))){
    temp <- dataframe[item]
    colnames(temp) <- "A"
    temp<-temp$A
    
    n <- normal(temp)
    ln <- log_normal(temp)
    p <- pearsonIII(temp)
    lp <- log_personIII(temp)
    gum <- gumbell(temp)
    errors <- c(n$error, ln$error, p$error, lp$error, gum$error)
    min_error <- min(errors)
    if(min_error == errors[1]){
      datos <- n$data
    }
    if(min_error==errors[2]){
      datos <- ln$data
    }
    if(min_error==errors[3]){
      datos <- p$data
    }
    if(min_error==errors[4]){
      datos <- lp$data
    }
    if(min_error==errors[5]){
      datos <- gum$data
    }
  extrapolation <- cbind(extrapolation, i = datos)
  err_punt <- c(err_punt, min_error)
  }
  colnames(extrapolation) <- c('T', coord$subcuenca)
  return(list(data = extrapolation, list_errores = err_punt))
}

p24_extrapolation <- get_extrapolation(p24)$data
errores <- get_extrapolation(p24)$list_errores

###Correccion de extrapolacion
p24_extrapolation <- p24_extrapolation*1.13


get_mean_pp_basin <- function(subbasin, centroids, extrapolation){
  #Grills
  basin <- as_Spatial(st_union(subbasin))
  ext <- extent(basin)
  x_values<- seq(from = ext[1], to = ext[2], by = 100)
  y_values<- seq(from = ext[3], to = ext[4], by = 100)
  grd <- expand.grid(x= x_values, y= y_values)
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE
  crs(grd)<- crs(basin)
  # Assignation de localization a TR
  extp_trans <- data.frame(t(extrapolation[-1]))
  colnames(extp_trans) <- paste("TR",tr, sep = "")
  extp_trans <- cbind(extp_trans, as.data.frame(coordinates(centroids)))
  coordinates(extp_trans)<- ~coords.x1 + coords.x2
  crs(extp_trans) <- crs(grd)
  
  mean_pp_basin <- data.frame( Subcuenca = paste("W",centroides$DrainID,"0", sep = ""))
  for (returnTime in c(1:9)) {
    #IDW DE POR TIEMPO DE RETORNO
    temp <- extp_trans[returnTime]
    colnames(temp@data)<- c("TR")
    idw.p <-gstat::idw(temp$TR ~ 1, temp, grd)
    idw.p <- raster(idw.p)
    idw.p <- raster::mask(idw.p, basin)
    
    subbasin <- cbind(subbasin, id = c(1:nrow(subbasin)))
    
    #PROMEDIO POR CADA SUBCUENCA
    caden <- c()
    for (item in c(1:nrow(subbasin))) {
      sdx <- as_Spatial(st_union(filter(subbasin, id == item)))
      mean_prep <- extract(idw.p, sdx, fun=mean)
      caden <- c(caden, mean_prep)
    }
    mean_pp_basin <- cbind(mean_pp_basin, TR = caden)
  }
  colnames(mean_pp_basin) <- c("Subcuenca", paste("TR",tr, sep = ""))
  mean_pp_basin <- data.frame(t(mean_pp_basin[-1]))
  colnames(mean_pp_basin) <- c(coord$subcuenca)
  return(mean_pp_basin)
}

mean_pp_basin<-get_mean_pp_basin(subcuencas, centroides, p24_extrapolation)

######Creación de hietogramas

get_hietogramas <- function(cadena){
  duracion <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                120, 180, 240, 300, 360, 420, 480, 540, 600, 660,
                720, 780, 840, 900, 960, 1020, 1080, 1140, 1200,
                1260, 1320, 1380, 1440)
  
  #Intensidad y curvas IDF --- Modelo de Dick y Peschke
  duracion_min <- duracion/60
  
  x3 <- c()
  x2 <- c()
  y <- c()
  for (pp in c(1:length(cadena))) {
    pp_dp <- c()
    for (item in c(1:length(duracion))) {
      n <- cadena[pp]*(duracion[item]/1440)^0.25/duracion_min[item]
      pp_dp <- c(pp_dp, n)
      x2 <- c(x2, tr[pp])
    }
    y <- c(y, pp_dp)
    x3 <- c(x3, log10(duracion))
  }
  x2 <- log10(x2)
  y <- log10(y)
  
  ec_int <- data.frame(x3 = x3, x2 = x2, y = y)
  
  #Regresion lineal multiple
  reg <- lm(ec_int$y ~ ec_int$x3 + ec_int$x2)
  k <- as.numeric(10^reg$coefficients[1])
  n <- as.numeric(-reg$coefficients[2])
  m <- as.numeric(reg$coefficients[3])
  
  #Elaboracion de hietogramas
  tr_hieto <- c(2,5,10,50,100)
  hietograma <- data.frame(matrix(nrow = 24))
  colnames(hietograma)<- "Duration"
  for (tr_i in tr_hieto) {
    dur_hieto <- c()
    intensidad <- c()
    for (i in 1:24) {
      dur_hieto <- c(dur_hieto, i*60)
      intensidad <- c(intensidad, k*tr_i^m/(i*60)^n)
    }
    hietograma$Duration <- dur_hieto
    prof_ac <- intensidad*dur_hieto/60
    prof_inc<-c()
    for (i in 1:24) {
      if (i ==1){
        prof_inc <- prof_ac[i]
      }
      else{
        prof_inc <- c(prof_inc, prof_ac[i] - prof_ac[i-1])
      }
    }
    
    pp_bAlterno <-c()
    for (i in sort((1:12)*2, decreasing = TRUE)) {
      pp_bAlterno <- c(pp_bAlterno, prof_inc[i])
    }
    for (i in ((1:12)*2-1)) {
      pp_bAlterno <- c(pp_bAlterno, prof_inc[i])
    }
    
    hietograma <- cbind(hietograma, i=pp_bAlterno)
  }
  colnames(hietograma)<- c("Duration", paste("TR",tr_hieto, sep = ""))
  return(hietograma)
}


get_list_hietogramas <- function(data_mean){
  data_pp <- list()
  for (item in c(1:ncol(data_mean))) {
    temp <- mean_pp_basin[item]
    colnames(temp)<- "A"
    temp <- temp$A
    
     n_pp <- get_hietogramas(temp)
    data_pp[[names(data_mean[item])]] <- n_pp
  }
  return(data_pp)
}

hietogramas_24 <- get_list_hietogramas(mean_pp_basin)


openxlsx::write.xlsx(hietogramas_24, file = "data.xlsx")
