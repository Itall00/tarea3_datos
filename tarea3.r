# cargar librerias
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(readxl)
library(vctrs)
library(Kendall)
library(raster)
library(ggrepel)
library(lwgeom)

options(scipen = 999)
#aqui poner el path de la carpeta para correr todo sin cambiar a cada rato
#path <- ""
#path = "C:/Users/alanp/Documents/5to/cs datos espaciales/tarea3" 
path <- "/Users/itallo/Documents/GitHub/tarea3_datos"


# Cargar funciones --------------------------------------------------------
# Leemos el raster y lo cortamos al area de interes
km <- read_sf(paste0(path, "/cuenca.kml"))
km <- mutate(km, Description = "Rio Aconcagua", altura = 1021)
write_sf(km, paste0(path, "/cuenca.geojson"), overwrite=TRUE)
img.folder <- paste0(path, "/landsat")
files <- list.files(img.folder, pattern = "SR_B", full.names = TRUE)
imgs <- rast(files)

img_ext <- ext(imgs)

img_ext[2] - img_ext[1] # medimos ancho y largo de la imagen
img_ext[4] - img_ext[3]
st_crs(imgs) == st_crs(km)

img.crs <- st_crs(imgs)
st_crs(imgs)
v <- st_transform(x = km, crs = img.crs)

st_crs(imgs) == st_crs(v) # verificamos que se tenga el mismo crs
imgs.c <- crop(imgs, vect(v))

imgs.m <- mask(imgs, vect(v)) # creamos mascara de la cuenca

imgs.cm <- crop(imgs.m, vect(v))

imgs.r <- project(imgs.cm, vect(km), method = "near")
cuenca <- imgs.r[[5]] # porque 5???????
plot(cuenca)

# Importamos el LandCover Zhao
lc.file <- paste0(path, "/LC_CHILE_2014_b.tif")
lc <- rast(lc.file)

# Cortamos el Lancover al area de nuestra cuenca
lc.proj <- project(lc, cuenca, method = "near") # Metodo near para datos categoricos
lc.mask <- mask(lc.proj, cuenca)
lc.crop <- crop(lc.mask, cuenca)
plot(lc.crop)

"CR2MET2.5"

# carpeta donde se encuentran los archivos de PP CR2MET
pp.dir <- paste0(path, "/PPCR2MET2.5") # directorio datos de precipitación

files <- list.files(pp.dir, full.names = TRUE, pattern = "nc$")
files
pp <- rast(files, subds = "pr")
pp # datos de CR2MET almacenados en pp

# crear vector de fechas
fechas.pp <- seq(
  ymd("2000-01-01"), ## cambio de fecha
  ymd("2022-12-31"),
  by = "days"
)
# asignamos las fechas como nombre de las capas
names(pp) <- fechas.pp
length(fechas.pp) == length(names(pp))

# extraer valores dentro de la cuenca
extr <- terra::extract(pp, km)
as_tibble(extr)

################# paso a datos diarios la precipitacion de la cuenca
extr <- extr %>% .[, -which(names(.) == "ID")]
extr <- extr %>% drop_na()
pp.day <- extr %>% summarise_all(mean)
pp.day <- pp.day %>% pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "pp")
pp.day


# PP mensual
pp.month <- pp.day %>%
  mutate(
    fecha = as_date(fecha),
    fecha = floor_date(fecha, unit = "month")
  ) %>%
  group_by(fecha) %>%
  summarise(pp = sum(pp))
pp.month

# PP anual
pp.year <- pp.month %>%
  mutate(
    fecha = as_date(fecha),
    fecha = floor_date(fecha, unit = "year")
  ) %>%
  group_by(fecha) %>%
  summarise(pp = sum(pp))
pp.year

"MODIS"
# Importamos los datos modis del area y les ponemos el formato que necesitamos
dir <- paste0(path, "/modis")
separate_eight_day_composite <- function(x, fechas) {
  # dias que faltan para terminar el mes
  days.dif <- 1 + (days_in_month(fechas) - day(fechas))

  # si faltan menos de 8 dias entonces la imagen comprende datos que pertenecen al mes siguiente
  index <- grep(TRUE, days.dif < 8)

  # no considerar la ultima imagen del año
  index <- index[-length(index)]

  # seleccionar imagenes, fechas y dias de las imagenes que abarcan mas de un mes
  x_i <- x[[index]]
  fechas_img <- fechas[index]
  days.left <- days.dif[index]

  # raster vacio para guardar resultados
  r1 <- rast()

  # guardar las imagenes que abarcan solo 1 mes dentro de la composicion de 8 días
  r2 <- x[[-index]]

  # cantidad de imagenes que abarcan mas de un mes
  n <- length(fechas_img)

  # ciclo de 1 a n
  for (i in 1:n) {
    # obtener la fecha del último dia que abarca la serie de 8 días
    fecha2 <- fechas_img[i] %m+% days(7)

    # Multiplicamos la imagen por la proporcion de días que abarca de cada mes
    # mes actual
    x1 <- x_i[[i]] * (days.left[i] / 8)
    # mes siguiente
    x2 <- x_i[[i]] * (1 - days.left[i] / 8)
    # asignar fecha correspondiente a cada imagen
    names(x1) <- fechas_img[i]
    names(x2) <- fecha2
    # guardar las nuevas imagenes en un vector
    r1 <- c(r1, x1, x2)
  }
  r <- c(r1, r2)
  return(r)
}





# Obtenemos las imagenes mensuales a partir de diarias
daily_to_monthly <- function(x, dates, fun = "mean") {
  mes <- floor_date(as_date(dates), unit = "month")
  lista.meses <- mes %>% unique()
  n <- length(lista.meses)
  x.mensual <- rast()

  for (i in 1:n) {
    posicion <- grep(lista.meses[i], mes)
    if (fun == "mean") {
      r <- mean(x[[posicion]], na.rm = TRUE)
    } else if (fun == "sum") {
      r <- sum(x[[posicion]], na.rm = TRUE)
    } else {
      message("No se reconoce la funcion utilizada, escriba: fun = sum o mean como argumento")
    }
    x.mensual <- c(x.mensual, r)
  }
  names(x.mensual) <- lista.meses

  return(x.mensual)
}

# Obtenemos las imagenes anuales a partir de mensuales
to_yearly <- function(x, dates, fun = "mean") {
  y <- floor_date(as_date(dates), unit = "year")
  lista.y <- y %>% unique()
  n <- length(lista.y)
  x.anual <- rast()

  for (i in 1:n) {
    posicion <- grep(lista.y[i], y)
    if (fun == "mean") {
      r <- mean(x[[posicion]], na.rm = TRUE)
    } else if (fun == "sum") {
      r <- sum(x[[posicion]], na.rm = TRUE)
    } else {
      message("No se reconoce la funcion utilizada, escriba: fun = sum o mean como argumento")
    }
    x.anual <- c(x.anual, r)
  }
  names(x.anual) <- lista.y

  return(x.anual)
}

files <- list.files(dir, full.names = TRUE, pattern = "_aid0001");files # datos evotranspiracion

et <- rast(files)

fechas.et <- names(et) %>%
  str_sub(start = 26, end = 32) %>%
  as.Date("%Y%j")
names(et) <- fechas.et

# separar las imagenes que abarcan mas de un mes dentro de los 8 dias del composite
et <- separate_eight_day_composite(et, fechas.et)

# actualizar vector de fechas
fechas.et <- as_date(names(et))

et.m <- daily_to_monthly(et, dates = fechas.et, fun = "sum")
et.m
et.y <- to_yearly(et.m, dates = names(et.m), fun = "sum")
et.y

hist(et.y[[1]])
summary(et.y)
et.y[et.y > 24000] <- NA
hist(et.y[[1]], breaks = 20)

hist(et.m[[1]])
summary(et.m)
et.m[et.m > 4000] <- NA
hist(et.m[[1]], breaks = 20)


# Tabla con todos los valores del raster
values.lc <- values(lc.crop, dataframe = TRUE) %>%
  drop_na()
#values.lc
# calcular numero de pixeles por valor
cat.npix <- table(values.lc) %>%
  as.vector()

# reproyectar raster
lc.r <- terra::project(lc.crop, crs(et.y), method = "near")
lc.r
plot(lc.r, main = "Land Cover reproyectado")







rcl <- c(
  100, 150, 1,
  210, 232, 2,
  241, 252, 3,
  310, 330, 4,
  410, 450, 5,
  500, 1210, 6
)
rcl
rcl.mat <- matrix(rcl, ncol = 3, byrow = TRUE)
rcl.mat
lc.cats <- tibble(
  ID = 1:6,
  nombre = c("Cultivos", "Bosque Nativo", "Plantaciones", "Praderas", "Matorrales", "Otros")
)
lc.cats
lc.reclass <- classify(lc.r, rcl.mat, right = NA)
levels(lc.reclass) <- lc.cats

plot(lc.reclass, main = "LandCover Reclasificado", col = c("yellow", "purple", "red", "blue", "green", "white"))

# cuanto mas grande son los pixeles de ET de modis con respecto a la resolucion del LandCover?
fac <- res(et.y)[1] / res(lc.reclass)[1]
fac
# cambiar resolucion a un raster
lc.agg <- aggregate(lc.reclass, fact = fac, fun = "modal")

plot(lc.reclass, main = "LandCover resolucion original")
plot(lc.agg, main = "LandCover de baja resolucion")

# resamplear para que los raster coincidan pixel a pixel
lc.r <- resample(lc.agg, et.y, method = "near")
plot(lc.r, main = "LandCover resampleado", col = c("yellow", "purple", "red", "blue", "green", "white"))

# Consumo medio por cobertura
etr_mean <- zonal(et.y, lc.r, fun = "mean", na.rm = TRUE) %>%
  pivot_longer(cols = 2:17, names_to = "fecha", values_to = "ET") %>%
  mutate(fecha = as_date(fecha))
etr_mean

# graficar serie de tiempo por cobertura
ggplot(etr_mean, aes(x = fecha, y = ET, color = nombre)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Año", y = "Etr (mm)", title = "Consumo medio por cobertura de suelo",
    color = "Cobertura"
  )

# Consumo total por cobertura
etr_total <- zonal(et.y, lc.r, fun = "sum", na.rm = TRUE) %>%
  pivot_longer(cols = 2:17, names_to = "fecha", values_to = "ET") %>%
  mutate(fecha = as_date(fecha))
etr_total

# graficar serie de tiempo por cobertura
ggplot(etr_total, aes(x = fecha, y = ET, color = nombre)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Año", y = "Etr (mm)", title = "Consumo total por cobertura de suelo",
    color = "Cobertura"
  )

# Asignar ET anual de las plantaciones a los pixeles de matorral
# raster vacio para guardar et modificada
et.mod <- rast()
# numero de capaz sobre las que iterar
n <- nlyr(et.y)
# vector con la ETr anual de las plantaciones forestales
et_pf <- etr_mean %>%
  filter(nombre == "Matorrales") %>%
  pull(ET)
et_pf

# ciclo de 1 a n
for (i in 1:n) {
  # seleccionar imagen i de ETr
  img_i <- et.y[[i]]
  # modificar ETr de las plantaciones a los pixeles que tienen categoria 4 en el LC (matorrales)a
  img_i[lc.r == 4] <- et_pf[i]
  # guardar la imagen modificada junto con las anteriores
  et.mod <- c(et.mod, img_i)
}

# plot ET original
plot(et.y[[1]], main = paste0(names(et.y)[1], " ET original"))
# plot ET modificada
plot(et.mod[[1]], main = paste0(names(et.y)[1], " ET modificada"))
plot(mask(et.mod[[1]], lc.r), main = paste0(names(et.y)[1], " ET modificada"))

# Media de todos los pixeles por categoria
etrmod_mean <- zonal(et.mod, lc.r, fun = "mean", na.rm = TRUE)
etrmod_mean

# Consumo total de todos los pixeles por cobertura
etrmod_total <- zonal(et.mod, lc.r, fun = "sum", na.rm = TRUE) %>%
  pivot_longer(cols = 2:17, names_to = "fecha", values_to = "ET") %>%
  mutate(fecha = as_date(fecha))
etrmod_total

# graficar serie de tiempo por cobertura
ggplot(etrmod_total, aes(x = fecha, y = ET, color = nombre)) +
  geom_line(linewidth = 1) +
  labs(
    x = "A?o", y = "Etr (mm)", title = "Consumo total por cobertura de suelo",
    color = "Cobertura"
  )

area_pixel <- 500 * 500 # m2

et.y.m3 <- et.y * (0.001 * area_pixel)
et.mod.m3 <- et.mod * (0.001 * area_pixel)
plot(et.mod.m3)
plot(et.y.m3)

# Consumo total de todos los pixeles por cobertura
etrmod_total <- zonal(et.mod.m3, lc.r, fun = "sum", na.rm = TRUE) %>%
  pivot_longer(cols = 2:17, names_to = "fecha", values_to = "ET") %>%
  mutate(fecha = as_date(fecha))
etrmod_total

etrmod_total <- etrmod_total %>%
  group_by(fecha) %>%
  mutate(
    prop = 100 * ET / sum(ET)
  )

# graficar serie de tiempo por cobertura
ggplot(etrmod_total, aes(x = fecha, y = ET, color = nombre)) +
  geom_line(linewidth = 1) +
  labs(
    x = "A?o", y = "Etr (m3/a?o)", title = "Consumo total por cobertura de suelo",
    color = "Cobertura"
  )

etrmod_total <- etrmod_total %>%
  group_by(fecha) %>%
  mutate(
    prop = 100 * ET / sum(ET)
  )

# grafico de torta (pie)
etrmod_total %>% ggplot(aes(x = "", y = prop, fill = nombre)) +
  geom_bar(
    stat = "identity", width = 1,
    color = "white"
  ) +
  coord_polar("y", start = 0) +
  theme_void() +
  facet_wrap(. ~ fecha, drop = TRUE) + # crea un grafico por cada valor unico de fecha
  scale_fill_brewer(palette = "Set1") + # remove background, grid, numeric labels
  labs(fill = "Cobertura", title = "Evapotranspiracion real anual por cada cobertura de suelo")

# Graficar solo un año
etrmod_total %>%
  filter(fecha == ymd("2010-01-01")) %>%
  ggplot(aes(x = "", y = prop, fill = nombre)) +
  geom_bar(
    stat = "identity", width = 1,
    color = "white"
  ) +
  coord_polar("y", start = 0) +
  theme_void() +
  facet_wrap(. ~ fecha, drop = TRUE) +
  scale_fill_brewer(palette = "Set1") + # remove background, grid, numeric labels
  labs(fill = "Cobertura", title = "Evapotranspiracion real anual por cada cobertura de suelo")

# exportar tabla con datos de etr total por cobertura
# dir.create("resultados")
# write_csv(etrmod_total, "resultados/evapotranspiracion_total_por_cobertura_modificada.csv")

# Etr anual de la cuenca modificada para el balance hidrico

extr <- terra::extract(et.mod, km)
et.year.mod <- extr %>%
  dplyr::select(-ID) %>%
  drop_na() %>%
  summarise_all(median) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "et_mod") %>%
  mutate(fecha = as_date(fecha))
et.year.mod

# Etr anual de la cuenca riginalpara el balance hidrico
extr <- terra::extract(et.y, km)
et.year <- extr %>%
  dplyr::select(-ID) %>%
  drop_na() %>%
  summarise_all(median) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "et") %>%
  mutate(fecha = as_date(fecha))
et.year


"LANDCOVER ZHAO"
# Proyectamos las categorias solicitadas en nuestro raster segun el LandCover Zhao
# Consideramos solo las categorias de cultivos(100), bosques(200), plantaciones(250), praderas(300) y matorrales(400)

# Modificamos pixeles de un raster seleccionaremos todos los pixeles diferentes a 100, 200, 300 y 400 y les asiga NA
# Tambien asignamos los pixeles especificos de cada seccion del LandCover a solo uno por categoria general
lc.crop[lc.crop >= 500] <- NA
lc.crop[lc.crop < 100] <- NA
lc.crop[lc.crop < 200 & lc.crop >= 100] <- 1
lc.crop[lc.crop < 235 & lc.crop >= 200] <- 2
lc.crop[lc.crop < 300 & lc.crop >= 235] <- 3
lc.crop[lc.crop < 400 & lc.crop >= 300] <- 4
lc.crop[lc.crop < 500 & lc.crop >= 400] <- 5

writeRaster(lc.crop, paste0(path, "/LC_reclasificado.tif"), overwrite = TRUE)

plot(lc.crop,
  col = c("yellow", "purple", "red", "blue", "green"),
  main = "Landcover por categorias"
) # Ploteamos lc cortado, cambiamos los colores a las categorias

lc.crop_vec <- as.vector(lc.crop) # Convertimos lc.crop en un vector

color_counts <- table(lc.crop_vec) # Obtenenemos la tabla de cantidades de las categor?as de color

# Creamos un gr?fico de barras de la cantidad por categoria
text(
  x = barplot(table(lc.crop_vec),
    col = c("yellow", "purple", "red", "blue", "green"),
    main = "Cantidad por categoria",
    xlab = "Categoria", ylab = "Cantidad",
    ylim = c(0, 1200000)
  ),
  y = table(lc.crop_vec),
  labels = table(lc.crop_vec),
  pos = 3
)

######## Cuenca de estudio -------------------------------------------------------

# leer shape de la cuenca
cuenca.t = read_sf(paste0(path, "/cuenca.kml")) %>%
  # reproyectar a coordenadas geográficas
  st_transform(4326)
######## Modelo de ElevaciÃ³n Digital (DEM) ---------------------------------------

# DEM
dem = rast(paste0(path, "/SRTMGL1_NC.003_SRTMGL1_DEM_doy2000042_aid0001.TIF"))

# graficar
plot(dem)
# explorar valores del raster
hist(dem)
summary(dem)
summary(values(dem))

# calcular derivadas topograficas
der_topo = terrain(dem, v = c("slope","aspect","TPI","TRI","roughness","flowdir"))
der_topo
# graficar
plot(der_topo)

# Reclasificar mapa de elevacion
# crear matriz de reclasificaciÃ³n
rcl.matrix = matrix(
  c(-1, 250,1,
    250,Inf,2),
  ncol = 3, byrow = TRUE)

# reclasificar
dem.rec = classify(dem, rcl.matrix)

# tabla con categorias (Look Up Table)
lut = tibble(ID = 1:2, ladera = c("Bajo","Alto"))

# asignar categorias a los valores del raster
levels(dem.rec) = lut
plot(mask(dem.rec, cuenca.t))

# agregar dem al stack de imagenes de derivadas topograficas
der_topo$dem = dem.rec
plot(der_topo)

######## Topographic Position Index ----------------------------------------------
hist(der_topo$TPI)
summary(values(der_topo$TPI))

# generar breaks para una escala de cuantiles
quantile_scale = function(x, quantiles  = c(0.0001, 0.25, 0.5, 0.75, 0.9999)){
  q = c()
  for (i in 1:length(quantiles)) {
    q = c(q, quantile(x, quantiles[i], na.rm =TRUE))
  }
  return(q)
}

# plotear TPI con escala de cuantiles
qscale = quantile_scale(values(der_topo$TPI))
plot(der_topo$TPI, breaks = qscale)

# Mejorar resultado cambiando el tamaño de la ventana movil
# tamaño de la ventana
n = 11
# matriz con la posicion de cada celda de la ventana movil
f <- matrix(1:(n*n), nrow=n, ncol=n)
f
# cual es el indice del pixel del centro?
centro = f[ceiling(n/2),ceiling(n/2)]
centro
# crear matriz con solo valores 1
f <- matrix(1, nrow=n, ncol=n)
f
# aplicar una funcion con una ventana movil con la funcion focal
# el valor de cada celda de la matriz representa su peso
TPI <- focal(dem, w=f, fun=function(x, ...) x[centro] - mean(x[-centro], na.rm = TRUE))
hist(TPI)
summary(values(TPI))

# plotear usando una escala de colores en cuantiles
qscale = quantile_scale(values(TPI))
plot(TPI, breaks = qscale, main = str_c("TPI con ventana movil de ", n,"x",n," pixels"))

# exportar resultado
# direccion carpeta de salida
dir = paste0(path, "/derivadas_topo")
# crear carpeta
dir.create(dir)
# nombre del archivo
fname = str_c(dir, "/TPI_",n,".tif");fname
# guardar archivo tif de la imagen
writeRaster(TPI, fname, overwrite = TRUE)

# # crear mapa interactivo simple
# remotes::install_github("rstudio/leaflet")
# library(leaflet)
# plet(TPI)

# reemplazar TPI con la ventana movil corregida
der_topo$TPI = TPI
der_topo


######## ExposiciÃ³n --------------------------------------------------------------

# Extraer la imagen de exposicion
aspect = der_topo$aspect
plot(aspect, main = "ExposiciÃ³n")
hist(aspect)
summary(values(aspect))

# reclasificar exposicion en categorias Norte, Sur, Este y Oeste
# crear matriz de reclasificacion
rcl.matrix = matrix(
  c(0, 45,1,
    45,135,2,
    135,225,3,
    225,315,4,
    315,360, 1), 
  ncol = 3, byrow = TRUE)
# reclasificar en laderas
aspect.rec = classify(aspect, rcl.matrix)
# tabla con categorias (Look Up Table)
lut = tibble(ID = 1:4, ladera = c("N","E","S","O"))
# asignar categorias a los valores del raster
levels(aspect.rec) = lut
# observar valores y etiqueta
unique(aspect.rec)
# paleta de colores categorica
colores = terrain.colors(4, alpha = 0.8)
#graficar
plot(aspect.rec, col = colores, main = "ExposiciÃ³n reclasificada")

# reemplazar exposicion por la reclasificacion creada
der_topo$aspect = aspect.rec

####### Pendiente ---------------------------------------------------------------
slope = der_topo$slope
plot(slope, main = "Pendiente")
hist(slope)

# reclasificar exposicion en categorias Norte, Sur, Este y Oeste
# crear matriz de reclasificacion
rcl.matrix = matrix(
  c(0,5,1,
    5,15,2,
    15,30,3,
    30, Inf,4), 
  ncol = 3, byrow = TRUE)

# reclasificar en laderas
slope.rec = classify(slope, rcl.matrix)
plot(slope.rec)

# tabla con categorias (Look Up Table)
lut = tibble(ID = 1:4, 
             pendiente = c("Suave","Moderada","Pronunciada","Muy Pronunciada"))
# asignar categorias a los valores del raster
levels(slope.rec) = lut
# observar valores y etiqueta
unique(slope.rec)
# paleta de colores categorica
colores = terrain.colors(4, alpha = 0.8)
#graficar
plot(slope.rec, col = colores, main = "Pendiente reclasificada")
# reemplazar pendiente por la reclasificacion creada
der_topo$slope = slope.rec
plot(der_topo)


####### Exportar derivadas topográficas -----------------------------------------
# nombre del archivo
fname = str_c(dir, "/derivadas_topograficas.tif");fname
# guardar archivo tif de la imagen
writeRaster(der_topo, fname, overwrite = TRUE, datatype = "INT2S")


######## Evapotranspiración real -------------------------------------------------

# carpeta donde se encuentran los archivos de evapotranspiracion real de MODIS
dir = paste0(path, "/modis")

# obtener la direccion de los archivos en la carpeta
files = list.files(dir, full.names = TRUE, pattern = "aid0001");files

# leer rasters
et = rast(files)
et

# crear vector de fechas
fechas.et = names(et) %>% 
  #recortamos los nombres  desde el caracter 26 al caracter 32.
  str_sub(start = 26, end = 32) %>%
  # le damos formato de fecha, %Y = YYYY, %j=ddd (numero de dia del año o dia juliano)
  as.Date("%Y%j")
fechas.et

# asignamos las fechas como nombre de las capas
names(et) = fechas.et

# separar las imagenes que abarcan mas de un mes dentro de los 8 dias del composite
et = separate_eight_day_composite(et, fechas.et);et

# actualizar vector de fechas
fechas.et = as_date(names(et))

# calcular imagenes mensuales y anuales de etr
et.m = daily_to_monthly(et, dates = fechas.et, fun = "sum");et.m
et.y = to_yearly(et.m, dates = names(et.m), fun = "sum");et.y

# histogramas
hist(et.y[[1]])

# resumen estadistico
summary(et.y)

# eliminar pixeles anormales
et.y[et.y > 1500] = NA

####### Land Cover --------------------------------------------------------------

# paleta de colores
colores = RColorBrewer::brewer.pal(7, "Set1");colores

# Resamplear LandCover para que coincida con los pixeles de ETr
lc.res = resample(lc.agg, et.y, method = "near")
plot(lc.res, main = 'Land Cover Cauquenes 2018 con agregacion', col = colores)
et.y
# Extraer pixeles de plantaciones forestales
lc.pf = lc.res
lc.pf[lc.res != 2] = NA
plot(lc.pf, col = c("#000000"))

###### ETr de plantaciones forestales segun su elevacion --------------------------

# Mapa de elevacion de plantaciones forestales
elev = resample(der_topo$dem, lc.pf, method = "near")

elev
lc.pf



##
st_crs(elev) == st_crs(lc.pf)

st_crs(elev)
st_crs(lc.pf)

lc.pf=lc.pf_reproyectado





# Exportamos lc.pf a un archivo GeoTIFF
writeRaster(lc.pf, "lc_pf2.tif")
# Exportamos lc.pf a un archivo GeoTIFF, sobrescribiendo si es necesario
writeRaster(lc.pf, paste0(path, "/lc_pf2.tif"), overwrite=TRUE)



# Exportamos lc.pf a un archivo GeoTIFF
writeRaster(elev, "elev2.tif")
# Exportamos lc.pf a un archivo GeoTIFF, sobrescribiendo si es necesario
writeRaster(elev, paste0(path, "/elev2.tif"), overwrite=TRUE)



###
library(sf)
library(terra)

# Reproyectamos lc.pf al CRS de elev (WGS 84)
lc.pf_reproyectado <- project(lc.pf, crs(elev))

# Verificamos el nuevo CRS
crs(lc.pf_reproyectado)

lc.pf=lc.pf_reproyectado
###

st_crs(lc.pf) == st_crs(elev) # verificamos que se tenga el mismo crs

## Aqui buscamos solucionar que el [mask] extents do not match

# Cargamos la biblioteca
library(terra)

# Recortamos elev para que tenga la misma extensión que lc.pf
elev_cropped <- crop(elev, ext(lc.pf))

# Ahora deberías poder enmascarar sin problemas
elev.pf = mask(elev_cropped, lc.pf)




#### Visualizamos el problema 


elev.pf = mask(elev, lc.pf)
plot(elev.pf, col = c("#377EB8","#E41A1C"),
     main = "Plantaciones forestales según su elevación")


# evapotranspiracion de plantaciones forestales clasificadas por pendiente
et.mean.pf = zonal(et.y, elev.pf, fun = "mean", na.rm = TRUE) %>% 
  pivot_longer(cols = 2:23, names_to = "fecha", values_to = "ET") %>% 
  mutate(fecha = as_date(fecha))

# grafico de et anual
ggplot(et.mean.pf)+
  geom_line(aes(x = fecha, y = ET, color = dem), linewidth = 0.8)+
  labs(y = "ETr (mm)", color = "Altura",
       title = "Evapotranspiración real anual de plantaciones forestales según su elevación")


# Modificar LandCover
# Eliminar pixeles del LC original

# copiar LandCover original
lc.orig = lc.res

# Eliminar pixeles de plantaciones forestales
lc.orig[lc.res == 2] = NA

# LandCover sin plantaciones
plot(lc.orig, main = "LC sin plantaciones")
as.numeric(lc.orig) %>% unique

# Agregar pixeles de Plantaciones reclasificados por altura
elev.pf[elev.pf == 1] = 7

# Transformar en un raster de numeros para visualizar los nuevos valores
elev.pf = as.numeric(elev.pf)
plot(elev.pf, main = "Valores de las plantaciones reclasificadas")

# Pegar ambas imagenes con la funcion mosaic
lc.mod = mosaic(lc.orig, elev.pf, fun = "max")

# Nueva tabla de clases
lut = tibble(cat = 0:7,
             name = c('Bosque Nativo','Agricultura de secano',
                      'Plantaciones Forestales Altas','Matorrales',
                      'Pastizales','Suelo Desnudo',
                      'Agricultura de riego', 'Plantaciones Forestales Bajas'))

# Asignar clases a los valores del raster
levels(lc.mod) = lut

# Graficar Land Cover
# Original
plot(lc.res, col = RColorBrewer::brewer.pal(7, "Set2"), main = "LC original")
# Modificado
plot(lc.mod, col = RColorBrewer::brewer.pal(8, "Set2"), main = "LC modificado")

# evapotranspiracion de plantaciones forestales clasificadas por pendiente
et.mean = zonal(et.y, lc.mod, fun = "mean", na.rm = TRUE) %>% 
  pivot_longer(cols = 2:23, names_to = "fecha", values_to = "ET") %>% 
  mutate(fecha = as_date(fecha))

# grafico de et anual
ggplot(et.mean)+
  geom_line(aes(x = fecha, y = ET, color = name), linewidth = 0.8)+
  labs(y = "ETr (mm)", color = "Altura",
       title = "Evapotranspiración real anual por cobertura de suelo")













# Modificar EvapotranspiraciÃ³n real ---------------------------------------

# Asignar ET anual de las plantaciones a los pixeles de matorral
lut
# raster vacio para guardar et modificada
et.mod = rast()
# numero de capaz sobre las que iterar
n = nlyr(et.y)
# vector con la ETr anual de las plantaciones forestales
et_pf = et.mean %>% 
  filter(name == 'Plantaciones Forestales Bajas') %>% 
  pull(ET);et_pf

# ciclo de 1 a n
for (i in 1:n) {
  # seleccionar imagen i de ETr
  img_i = et.y[[i]]
  # modificar ETr de las plantaciones a los pixeles que tienen categoria 4 en el LC (matorrales)a
  img_i[lc.mod == 3] = et_pf[i]
  # guardar la imagen modificada junto con las anteriores
  et.mod = c(et.mod, img_i)
}

# plot ET original
plot(et.y[[1]], main = paste0(names(et.y)[1], ' ET original'))
# plot ET modificada
plot(et.mod[[1]], main = paste0(names(et.y)[1], ' ET modificada'))
plot(mask(et.mod[[1]], lc.mod), main = paste0(names(et.y)[1], ' ET modificada'))

# Land Cover Modificado
lc.esc.1 = lc.mod
# Cambiar coberturas de matorral por plantacion forestal.
lc.esc.1[lc.mod == 3] = 7

# Transformar valores de mm a m3/aÃ±o
# calcular el area de la cuenca
# resolucion espacial del raster de ETr es de 500m
# res(et)*111*1000 #m
area_pixel = 500*500 #m2

# transformar unidades
et.y.m3 = et.y * (0.001*area_pixel)
et.mod.m3 = et.mod * (0.001*area_pixel)
plot(et.mod.m3)
plot(et.y.m3)

# serie de tiempo de cobertura modificada
etrmod_total = zonal(et.mod.m3, lc.esc.1, fun = 'sum', na.rm = TRUE) %>% 
  pivot_longer(cols = 2:23, names_to = 'fecha', values_to = 'ET') %>% 
  mutate(fecha = as_date(fecha))
etrmod_total

# graficar serie de tiempo por cobertura
ggplot(etrmod_total, aes(x = fecha, y = ET, color = name))+
  geom_line(linewidth = 1)+
  labs(x = 'AÃ±o', y = 'Etr (m3/aÃ±o)', title = 'Consumo total por cobertura de suelo',
       color = 'Cobertura')

etrmod_total = etrmod_total %>% group_by(fecha) %>% 
  mutate(
    prop = 100*ET/sum(ET)
  )

# grafico de torta (pie)
etrmod_total %>% ggplot(aes(x = '', y = prop, fill = name))+
  geom_bar( 
    stat = 'identity', width = 1,
    color="white"
  )+
  coord_polar("y", start=0)+
  theme_void()+
  facet_wrap(.~fecha, drop = TRUE)+ #crea un grafico por cada valor unico de fecha
  scale_fill_brewer(palette="Set1")+ # remove background, grid, numeric labels
  labs(fill = 'Cobertura', title = 'EvapotranspiraciÃ³n real anual por cada cobertura de suelo')

# Graficar solo un aÃ±o
etrmod_total %>% 
  filter(fecha == ymd("2021-01-01")) %>% 
  ggplot(aes(x = '', y = prop, fill = name))+
  geom_bar( 
    stat = 'identity', width = 1,
    color="white"
  )+
  coord_polar("y", start=0)+
  theme_void()+
  facet_wrap(.~fecha, drop = TRUE)+
  scale_fill_brewer(palette="Set1")+ # remove background, grid, numeric labels
  labs(fill = 'Cobertura', title = 'EvapotranspiraciÃ³n real anual por cada cobertura de suelo')

# exportar tabla con datos de etr total por cobertura
#dir.create("resultados")
#write_csv(etrmod_total, "resultados/evapotranspiracion_total_por_cobertura_modificada.csv")

plot(et.mod)
plot(cuenca.t)
# Etr anual de la cuenca modificada para el balance hidrico
extr = terra::extract(et.mod, cuenca.t)
et.year.mod = extr %>%
  select(-ID) %>% 
  drop_na() %>% 
  summarise_all(median) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "et_mod") %>% 
  mutate(fecha = as_date(fecha))
et.year.mod

# Etr anual de la cuenca riginalpara el balance hidrico
extr = terra::extract(et.y, cuenca.t)
et.year = extr %>%
  select(-ID) %>% 
  drop_na() %>% 
  summarise_all(median) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "et") %>% 
  mutate(fecha = as_date(fecha))
et.year

# Precipitacion CR2MET -------------------------------------------------------

# carpeta donde se encuentran los archivos de evapotranspiracion real de MODIS
pp.dir = "data/raster/PPCR2MET2.5"

# obtener la direccion de los archivos en la carpeta
files = list.files(pp.dir, full.names = TRUE, pattern = "nc$");files

# leer rasters
pp = rast(files)
pp

# crear vector de fechas
fechas.pp = seq(
  ymd("1960-01-01"),
  ymd("2021-12-31"),
  by = "days"
)

# asignamos las fechas como nombre de las capas
names(pp) = fechas.pp

# extraer valores dentro de la cuenca
extr = terra::extract(pp, cuenca.t)
as_tibble(extr)

# calcular precipitacion promedio de la cuenca
pp.day = extr %>%
  select(-ID) %>% 
  drop_na() %>% 
  summarise_all(mean) %>% 
  pivot_longer(cols = 1:ncol(.), names_to = "fecha", values_to = "pp")
pp.day

# PP mensual
pp.month = pp.day %>% 
  mutate(fecha = as_date(fecha),
         fecha = floor_date(fecha, unit = "month")) %>% 
  group_by(fecha) %>% 
  summarise(pp = sum(pp))
pp.month

# PP anual
pp.year = pp.month %>% 
  mutate(fecha = as_date(fecha),
         fecha = floor_date(fecha, unit = "year")) %>% 
  group_by(fecha) %>% 
  summarise(pp = sum(pp))
pp.year

# Datos descargados de la DGA
############


q.month <- readr::read_delim(paste0(path, "/Caudal_Aconcagua.csv"), delim = ";") %>%
  pivot_longer(cols = 2:13, names_to = "mes", values_to = "caudal") 

q.month

# vector de fechas mensuales
fechas = seq(ym("2000-01"), ym("2015-12"), by = "months")

# crear columna con fechas
q.month = q.month %>% 
  mutate(fecha = fechas) %>% 
  dplyr::select(fecha, caudal)
q.month


# calcular caudal medio anual
q.year = q.month %>% 
  mutate(fecha = floor_date(fechas,unit = "year")) %>% 
  group_by(fecha) %>% 
  summarise_all(mean)
q.year



# Juntar variables del Balance Hidrico ------------------------------------
data.year = full_join(pp.year, et.year.mod, by = "fecha") %>% 
  full_join(et.year, by = "fecha") %>% 
  full_join(q.year, by = "fecha") %>% 
  mutate(
    disp = pp-et-caudal,
    disp_mod = pp-et_mod-caudal
    
  )

# Balance a hidrico anual
# serie de tiempo de datos anuales
ggplot(data.year)+
  geom_line(aes(x = fecha, y = pp, color = "PrecipitaciÃ³n"), linewidth = 0.8)+
  geom_line(aes(x = fecha, y = et, color = "ETr"), linewidth = 0.8)+
  geom_line(aes(x = fecha, y = et_mod, color = "ETr_modificada"), linewidth = 0.8)+
  geom_line(aes(x = fecha, y = caudal,color = "Caudal"), linewidth = 0.8)+
  geom_point(aes(x = fecha, y = caudal,color = "Caudal"), linewidth = 0.8)+
  geom_line(aes(x = fecha, y = disp, color = "Disponibilidad"), linewidth = 0.8)+
  geom_line(aes(x = fecha, y = disp_mod, color = "Disponibilidad_mod"), linewidth = 0.8)+
  geom_hline(yintercept = 0, linewidt = 0.8, linetype = "dashed")+
  scale_x_date(limits = c(ymd("2000-01-01"), ymd("2021-12-31")),
               date_labels = "%Y", date_breaks = "2 year")+
  labs(x = "tiempo", y = "(mm)", title = "Serie de tiempo mensual de componentes del BH",
       subtitle = "Los aÃ±os sin medicion de caudal, faltan datos en algunos meses",
       color = "")

