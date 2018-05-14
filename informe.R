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

df <- rename(df, nombres_castellano)

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

indicadora <- function(z) { ifelse(abs(z) <= 1, 1, 0) }

nucleos <- list(
    'uniforme' = function(z) { indicadora(z) * (1/2) },
    'triangular' = function(z) { indicadora(z) * (1 - abs(z)) },
    'epanechnikov' = function(z) { indicadora(z) * (3/4) * (1 - z^2) },
    'gaussiano' = function(z) { (1/sqrt(2 * pi)) * exp(-z^2/2) }
)

# Grafiquémoslos para ver que tengan sentido
rango <- seq(-1.1, 1.1, by = 0.01)

tibble(
  x = rep(rango, each = length(nucleo)),
  nombre_nucleo = rep(names(nucleo), length(rango)),
  y =  map2_dbl(nombre_nucleo, x, ~nucleo[[.x]](.y))
) %>%
  ggplot(aes(x, y, color = nombre_nucleo)) +
  geom_line() +
  coord_fixed() -> fig2

ggsave("fig2_nucleos.png", fig2, width = 6, height = 2.5)
fig2

# Exploremos cómo estiman la densidad los distintos núcleos, con
# algunos valores tentativos de ancho de banda (h):

estimar_densidad <- function(X, nombre_nucleo, h, x) {
    n <- length(X)
    nucleo <- nucleos[[nombre_nucleo]]
    sum(nucleo((X - x) / h)) / (2 * n * h)
}

grilla <- min(df$elevacion):max(df$elevacion)
valores_h <- c(50, 200, 400)

for (h in valores_h) {
    tibble(
        elevacion = rep(grilla, length(names(nucleos))),
        nucleo = rep(names(nucleos), each = length(grilla)),
        densidad = map2_dbl(elevacion, nucleo, function(x, nucleo) {
            estimar_densidad(df$elevacion, nucleo, h, x)
        })
    ) %>%
        ggplot() +
        aes(x = elevacion, y = densidad, color = nucleo) +
        geom_line() +
        labs(title = glue('h = {h}')) +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) -> fig

    ggsave(glue('fig_densidad_con_h_{h}.png'), fig)
    print(fig)
}

# Como se aprecia en las figuras, los núcleos de Epanechnikov y triangular
# ofrecen resultados similares. El valor de h = 200 parece adecuado.

# Realicemos una regresión no paramétrica con el método Nadaraya-Watson.

construir_estimador_NW <- function(X, Y, nombre_nucleo, h) {
    function(x) {
        nucleo <- nucleos[[nombre_nucleo]]
        pesos <- nucleo((X - x) / h)
        sum(pesos * Y) / sum(pesos)
    }
}
