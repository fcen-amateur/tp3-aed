## ----inicialización, include = FALSE-------------------------------------
library(tidyverse, quietly = TRUE)
library(stringr)
library(knitr)

## ----ingreso de datos----------------------------------------------------
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
  geom_point() +
  labs(title = "Fig. 1: Temperatura anual promedio y elevación.",
        y = "Temperatura anual promedio", x = "Elevación") -> fig1

ggsave("fig1_dispersion.png", fig1)
fig1

## ------------------------------------------------------------------------
df %>%
  arrange(desc(elevacion)) %>%
  mutate(
    par_consecutivo = paste(estacion, lead(estacion), sep = "-"),
    delta_elevaciones = elevacion - lead(elevacion)) %>%
  select(par_consecutivo, delta_elevaciones) %>%
  arrange(desc(delta_elevaciones)) %>%
  head() -> distancias

kable(x = distancias, caption = "Tabla 1: Distancia entre pares consecutivos de estaciones.")

## ----eliminación de atípicos---------------------------------------------
df <- filter(df, estacion !="CLK", estacion !="HAR")

## ----núcleos-------------------------------------------------------------
indicadora <- function(z) { ifelse(abs(z) <= 1, 1, 0) }

nucleos <- list(
    'uniforme' = function(z) { indicadora(z) * (1/2) },
    'triangular' = function(z) { indicadora(z) * (1 - abs(z)) },
    'epanechnikov' = function(z) { indicadora(z) * (3/4) * (1 - z^2) },
    'gaussiano' = function(z) { (1/sqrt(2 * pi)) * exp(-z^2/2) } # exactamente equivalente a dnorm(z)
)

## ----gráfico de los núcleos----------------------------------------------
rango <- seq(-1.1, 1.1, by = 0.01)
aplicar_nucleo <- function(nombre_nucleo, x) { nucleos[[nombre_nucleo]](x) }

cross_df(.l = list(
  x = rango,
  nombre_nucleo = names(nucleos))) %>%
  mutate(y = pmap_dbl(list(nombre_nucleo, x), aplicar_nucleo)) %>%
  ggplot(aes(x, y, color = nombre_nucleo)) +
  geom_line() +
  labs(title = "Fig. 2: Núcleos considerados, en el intervalo [-1, 1]",
       color = "Núcleo 'K'", y = "K(x)") +
  coord_fixed() -> fig2

ggsave("fig2_nucleos.png", fig2, width = 6, height = 3.5)
fig2

## ------------------------------------------------------------------------
estimador_densidad <- function(X, nombre_nucleo, h, x) {
    n <- length(X)
    nucleo <- nucleos[[nombre_nucleo]]   
    sum(nucleo((X - x) / h)) / (2 * n * h)
}

## ----densidad elelevacion para h seleccionados---------------------------
elevacion <- df$elevacion
grilla <- min(elevacion):max(elevacion)
h_prueba <- c(50L, 200L, 400L)

# Construimos el estimador de densidad para los datos de elevación
estimador_densidad_elevacion <- function(...) {
  estimador_densidad(X = elevacion, ...)
}

# `cross_df` devuelve el producto cartesiano de los vectores de `.l`  en forma de dataframe (por eso el sufijo _df)
cross_df(.l = list(
  nombre_nucleo = names(nucleos),
  h = h_prueba,
  x = grilla)) %>%
  mutate(
    # pmap llama la función .f para cada tupla de elementos presente en .l (nombre_nucleo, h, x),
    # y con el sufijo _dbl devuelve un vector de valores de punto flotate.
    densidad = pmap_dbl(.l = list(nombre_nucleo, h, x),
                        .f = estimador_densidad_elevacion)) -> densidades

## ----comparacion_nucleos, fig.height = 12--------------------------------
agregar_h <- function(df, h) {df %>% mutate(h = h)}

map_df(h_prueba, ~agregar_h(df, .)) %>%
  ggplot(mapping = aes(x = elevacion, alpha = 0.1)) +
  geom_histogram(bins = 30) +
  geom_line(data = densidades, mapping = aes(x = x, y = densidad*10600, color = nombre_nucleo), inherit.aes = F) +
  facet_grid(h ~ .) +
  guides(alpha = F) +
  labs(title = "Fig. 3: Curvas de densidad para distintos núcleos y h",
       subtitle = "La superposición sobre el histograma permite comprender cómo h cambia la 'suavidad' de la estimación.",
       y = "Conteo / 10600 * K(x)", color = "Núcleo 'K'") -> fig3

ggsave("fig3_nucleos_h.png", fig3, width = 7, height = 12)
fig3

## ----implementación de la regresión no paramétrica-----------------------
# El estimador de la esperanza condicional de Y para X = x segun NW es un promedio ponderado de los Y
construir_estimador_NW <- function(X, Y, nombre_nucleo, h) {
  nucleo <- nucleos[[nombre_nucleo]]   
  function(x) {
    pesos <- nucleo((X - x) / h)
    
    return(weighted.mean(Y, pesos))
  }
}


construir_estimador_de_error_NW <- function(X, Y) {
  function(nombre_nucleo, h) {
    # Bajo LOOCV, el estimador de NW para usar en el i-esimo dato, 
    # debe ser entrenado en todos los datos salvo el i-ésimo
    estimador_iesimo <- function(i) { construir_estimador_NW(X[-i], Y[-i], nombre_nucleo, h) }
    estimadores <- map(seq_along(X), estimador_iesimo)
    
    # Las predicciones resultan de aplicar el i-ésimo estimador, al i-ésimo valor de X.
    predicciones <- pmap_dbl(
      list(estimador = estimadores, x = X),
      function(estimador, x) { estimador(x) }
    )
    
    return (mean((Y - predicciones)^2))
    }
}

# Parametrizamos el estimador con los datos de elevación y temp_anual
estimador_de_error_NW <-
  construir_estimador_de_error_NW(X = df$elevacion, Y = df$temp_anual)

## ----ecm_por_nucleo_y_h--------------------------------------------------
cross_df(.l = list(nombre_nucleo = names(nucleos), h = 1:1000)) %>%
  mutate( ecm = pmap_dbl(.l = list(nombre_nucleo, h),
                        .f = estimador_de_error_NW)) %>%
  filter(!is.na(ecm)) -> errores

errores %>%
  ggplot(aes(x = h, y = ecm, color = nombre_nucleo)) +
  geom_line() +
  labs(title = "Fig. 4: Comportamiento comparativo de los núcleos en función de h",
       y = "Error Cuadrático Medio", color = "Núcleo 'K'") -> fig4

ggsave("fig4_ecm_por_nucleo_y_h.png", fig4)
fig4

## ------------------------------------------------------------------------
errores %>%
  group_by(nombre_nucleo) %>%
  filter(rank(ecm) == 1) %>% 
  ungroup() -> minimos_por_nucleo

kable(minimos_por_nucleo, digits = 5)

## ----minimos ecm por nucleo----------------------------------------------
minimos_por_nucleo <- minimos_por_nucleo %>%
  mutate(etiqueta = str_c("ECM: ", round(ecm, 3), " (h = ", h, ")"))

fig4 +
  geom_vline(data = minimos_por_nucleo,
             mapping = aes(xintercept = h, color = nombre_nucleo), lty = "dashed", size = 0.3) +
  geom_text(data = minimos_por_nucleo,
            mapping = aes(label = etiqueta, x = h, y = c(0.775, 0.75, 0.725, 0.7)), size = 2.5, color = "black") +
  coord_cartesian(xlim = c(50,400), ylim = c(0.7, 1.2)) +
  labs(title = "Fig. 5: Comportamiento comparativo, enfocado y anotado") -> fig5

ggsave("fig5_h_optimo_por_nucleo.png", fig5)
fig5

## ----estimaciones por lm-------------------------------------------------

polinomio_dados_coeficientes <- function(coeficientes) {
  function(x) {
    # El polinomio tendra tantos terminos como coeficientes se pasen
    terminos <- vector("numeric", length(coeficientes))
    for (n in seq_along(coeficientes)) {
      # El n-ésimo término es el producto del n-ésimo coeficiente y x elevado a la n-1
      terminos[[n]] <- coeficientes[[n]] * x^(n-1)
    }
    return(sum(terminos))
  }
}

evaluar_modelo_lineal_LOOCV <- function(X, Y, grado) {
  estimador_iesimo <- function(i) {
    # Bajo LOOCV, el estimador de NW para usar en el i-esimo dato, 
    # debe ser entrenado en todos los datos salvo el i-ésimo
    modelo <- lm(Y ~ poly(X , degree = grado, raw = T),
         list(X = X[-i], Y = Y[-i]))
    
    # El estimador no es más que el polinomio de grado `grado` con los coeficientes ajustados
    return ( polinomio_dados_coeficientes(modelo$coefficients) )

  }
  # Al i-ésimo valor de X le corresponde el estimador i-ésimo
  estimadores <- map(seq_along(X), estimador_iesimo)
  
  predicciones <- pmap_dbl(
    list(estimador = estimadores, x = X),
    function(estimador, x) { estimador(x) }
  )

  # Devolvemos el promedio de los length(X) ECMs
  return (mean((Y - predicciones)^2))
}

evaluar_modelo_lineal_en_g <- function(g) {
  evaluar_modelo_lineal_LOOCV(X = df$elevacion, Y = df$temp_anual, grado = g)
}

errores_lineales <- tibble(
  grado = 1L:12L,
  ecm = map_dbl(grado, evaluar_modelo_lineal_en_g))

kable(errores_lineales, digits = 5)

## ------------------------------------------------------------------------
minimo_error_lineal <- filter(errores_lineales, near(ecm, min(ecm)))

errores_lineales %>%
  ggplot(aes(grado, ecm)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 1:12, minor_breaks = NULL) +
  geom_point(minimo_error_lineal, mapping = aes(y = ecm, x = grado), color = "red") +
  geom_text(minimo_error_lineal, mapping = aes(
    y = ecm + 0.2, x = grado, 
    label = str_c("ECM: ", round(ecm, 3), " (g = ", grado, ")"))) +
    labs(title = "Fig. 6: ECM de la regresión lineal según el grado del polinomio",
       y = "Error Cuadrático Medio", x = "Grado 'g' del polinomio") -> fig6

ggsave("fig6_ecm_modelo_lineal.png", fig6)
fig6


## ------------------------------------------------------------------------
minimos_por_nucleo %>% 
  mutate(modelo = str_c("NW c/ núcleo ", nombre_nucleo, " (h = ", h, ")")) %>%
  select(modelo, ecm) -> resumen_no_parametricos

errores_lineales %>%
  filter(grado <= 3) %>%
  mutate(modelo = str_c("Polinomio de grado ", grado)) %>%
  select(modelo, ecm) -> resumen_polinomios

resumen_no_parametricos %>%
  bind_rows(resumen_polinomios) %>%
  arrange(ecm) -> resumen_general

kable(resumen_general, digits = 5)

## ----gráfico encimado de los modelos-------------------------------------

estimador_no_parametrico <- function(nombre_nucleo, h) {
  construir_estimador_NW(df$elevacion, df$temp_anual, nombre_nucleo, h)
}

estimador_polinomico <- function(g){
  modelo <- lm(temp_anual ~ poly(elevacion, g, raw = TRUE), df)
  return(polinomio_dados_coeficientes(modelo$coefficients))
}

estimadores_no_parametricos <- minimos_por_nucleo %>%
  mutate(
    modelo = str_c("NW c/ núcleo ", nombre_nucleo, " (h = ", h, ")"),
    estimador = map2(nombre_nucleo, h, estimador_no_parametrico))
  
estimadores_polinomicos <- tibble(
  grado = 1:3,
  modelo = str_c("Polinomio de grado ", grado),
  estimador = map(1:3, estimador_polinomico))

estimadores <- bind_rows(estimadores_no_parametricos, estimadores_polinomicos)

agregar_x <- function(df, x) {df %>% mutate(x = x)}
aplicar_estimador <- function(estimador, x) { estimador(x) }

curvas_estimadores <-
  map_df(min(elevacion):max(elevacion), ~agregar_x(estimadores, .)) %>%
  mutate(y = map2_dbl(estimador, x, aplicar_estimador))


curvas_estimadores %>%
  filter(grado == 3 | nombre_nucleo == "uniforme" | nombre_nucleo == "epanechnikov") %>%
  ggplot(aes(x, y, color = modelo)) +
  geom_line() +
  geom_point(data = df, mapping = aes(x = elevacion, y = temp_anual), inherit.aes = F) +
  labs(title = "Fig. 7: Datos y curvas para los estimadores seleccionados.",
       x = "Elevación", y = "Temperatura Anual", color = "Modelo") -> fig7

ggsave("fig7_mejores_modelos.png", fig7)
fig7

