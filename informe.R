# Librerías
library(tidyverse, quietly = TRUE)

# Lectura de datos y dispersión
df <- read_csv(
  file = "ortann.csv",
  col_types = cols(
    station = col_character(),
    latitude = col_double(),
    longitude = col_double(),
    elevation = col_integer(),
    tann = col_double()
))

nombres_castellano <- c(
  "station" = "estacion",
  "latitude" = "latitud",
  "longitude" = "longitud",
  "elevation" = "elevacion",
  "tann" = "temp_anual")

df <- plyr::rename(df, nombres_castellano)

df %>%
  ggplot(aes(x = elevacion, y = temp_anual)) +
  geom_point() -> fig1

ggsave("fig1_dispersion.png", fig1)
fig1

# Estimación paramétrica
lin <- lm(temp_anual ~ elevacion, df)
cuad <- lm(temp_anual ~ poly(elevacion, 2), df)

# Estimación no paramétrica

#' Va a ser más fácil trabajar con vectores que tibbles, así que extraigo las columnas relevantes
X <- df$elevacion
Y <- df$temp_anual

#' Listo primero los núcleos que mencionamos en clase.
nucleo <- list(
  "epanechnikov" = function(u) { ifelse( abs(u) <= 1, (3/4) * (1 - u^2), 0 ) },
  "triangular"   = function(u) { ifelse( abs(u) <= 1, 1 - abs(u)       , 0 ) },
  "uniforme"     = function(u) { ifelse( abs(u) <= 1, 1                , 0 ) },
  "gaussiano"    = function(u) { 1 / (sqrt(2*pi)) * exp(-u^2/2) }
)

# Grafiquémoslos para ver que tengan sentido
rango <- seq(-1.1, 1.1, by = 0.01)

tibble(
  x = rep(rango, each = length(nucleo)),
  nombre_nucleo = rep(names(nucleo), length(rango)),
  y =  map2_dbl(nombre_nucleo, x, ~nucleo[[.x]](.y))
) %>%
  ggplot(aes(x, y, color = nombre_nucleo)) +
  geom_line() -> fig2

ggsave("fig2_nucleos.png", fig2)
fig2

# Ahora vamos a necesitar un "generador de generadores" de funciones de densidad:

generador_estimador_densidad <- function(nucleo, h) {
  #' Dado un nucleo y un ancho de banda 'h', devuelve una funcion que..
  function(X) {
    #' dado un vector X, devuelve una función que...
    n <- length(X)
    function(x) {
      #' computa la densidad de X en x
      1 / (2 * n * h) * sum(nucleo( (X - x) / h ))
    }
  }
}

#' estimador_densidad_epanech es una funcion que dada una muestra X, devuelve
#' una funcion de densidad empirica usando un nucleo epanechnikov y h=200
estimador_densidad_epanech <- generador_estimador_densidad(nucleo$epanechnikov, 200)
#' usamos la funcion anterior, para generar la ecdf (empirical cumulative density function) de X
densidad_epanech_elevacion <- estimador_densidad_epanech(X)

#' Y ploteamos a ver qué pasa
data <- tibble(
  x = -100:2100,
  y = map_dbl(x, densidad_epanech_elevacion)
)

data %>% 
  ggplot(aes(x,y)) +
  geom_line() -> fig3

ggsave("fig3_densidad_epanech_elevacion.png", fig3)
fig3

#' ¡Razonable!