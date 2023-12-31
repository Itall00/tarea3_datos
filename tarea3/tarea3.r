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

path = "C:/Users/alanp/Documents/5to/cs datos espaciales/tarea3" #aqui poner el path de la carpeta para correr todo sin cambiar a cada rato
#path <- "/Users/itallo/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Tarea2_final"
# Cargar funciones --------------------------------------------------------
# Leemos el raster y lo cortamos al area de interes
km <- read_sf(paste0(path, "/polygon.kml"))
km <- mutate(km, Description = "Rio Aconcagua", altura = 1021)
write_sf(km, paste0(path, "/polygon.geojson"))
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






