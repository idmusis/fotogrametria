# Aplicativo para estimativa de velocidade veicular
# Autor Carlo Ralph De Musis
# Versão: 0.85

# Pacotes ------------------
pacman::p_load(
  # Interface Shiny
  shiny,
  shinydashboard,
  shinyBS, # popovers

  # Manipulação de imagem
  base64enc,
  jpeg,
  png,
  imager, # suavização de imagem

  # Estatística e modelagem ----------
  MASS, # simulação de distribuição normal multivariada (mvrnorm)
  MVN, # testes de normalidade multivariada (ex: Mardia)
  deming, # regressão Deming (erro em x e y)

  # Manipulação de dados -------------
  dplyr,
  tidyr,

  # Visualização ----------------------
  ggplot2,
  grid,
  hexbin, # gráficos de densidade hexagonal
  plotly # gráficos interativos
)

# Objetos ------------------------------------

## Configurações do painel ------
config_plotly <- function(p) {
  plotly::config(p,
    scrollZoom = FALSE,
    locale = "pt-BR",
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("select", "select2d", "lasso2d")
  )
}

input_com_ajuda <- function(input_id, label_text, input_ui, ajuda_id, ajuda_titulo = NULL, ajuda_texto) {
  tagList(
    tags$div(
      style = "display: flex; align-items: center;",
      class = "form-group shiny-input-container",
      tags$label(label_text, class = "control-label"),
      shinyBS::bsButton(ajuda_id, label = NULL, icon = icon("info-circle"), size = "extra-small")
    ),
    input_ui,
    shinyBS::bsPopover(
      ajuda_id, ajuda_titulo, ajuda_texto,
      placement = "right",
      options = list(container = "body", trigger = "hover focus")
    )
  )
}

titulo_com_ajuda <- function(ajuda_id, label_text, ajuda_texto, ajuda_titulo = NULL, centralizar = FALSE) {
  alinhamento <- if (centralizar) "center" else "flex-start"

  tagList(
    tags$div(
      style = paste0("display: flex; align-items: center; justify-content: ", alinhamento, "; gap: 8px;"),
      tags$h4(label_text),
      shinyBS::bsButton(ajuda_id, label = NULL, icon = icon("info-circle"), size = "extra-small")
    ),
    shinyBS::bsPopover(
      id = ajuda_id,
      title = ajuda_titulo,
      content = ajuda_texto,
      placement = "right",
      options = list(container = "body", trigger = "hover focus")
    )
  )
}

## Funções de processamento -----
# Gera repeticoes usando uma distribuicao normal bivariada

gerar_repeticoes <- function(dados_originais, repeticoes) {
  grupos <- unique(dados_originais$ponto) # Identificar os grupos de pontos
  repeticoes_geradas <- list() # Lista para armazenar as repeticoes

  for (grupo in grupos) {
    dados_grupo <- dados_originais[dados_originais$ponto == grupo, ] # Dados do grupo atual
    media <- colMeans(dados_grupo[, c("x", "y")]) # Media do grupo atual
    matriz_cov <- cov(dados_grupo[, c("x", "y")]) # Matriz de covariancia do grupo atual

    repeticoes_grupo <- lapply(1:repeticoes, function(i) {
      repeticao <- MASS::mvrnorm(1, mu = media, Sigma = matriz_cov)
      data.frame(serie = i, ponto = grupo, x = repeticao[1], y = repeticao[2])
    })

    repeticoes_geradas <- c(repeticoes_geradas, repeticoes_grupo)
  }
  # Junta todas as repeticoes em um unico data frame
  repeticoes_geradas <- do.call(rbind, repeticoes_geradas)

  return(repeticoes_geradas)
}

## ______________________________________________________________________________

hex_cinza <- function(hex_colors) {
  # Aplica a conversão a cada cor no vetor usando sapply
  cinza_values <- sapply(hex_colors, function(hex_color) {
    # Converte hexadecimal para RGB
    rgb <- grDevices::col2rgb(hex_color)

    # Converte RGB para cinza usando a fórmula de luminosidade
    cinza <- (rgb[1, ] + rgb[2, ] + rgb[3, ]) / 3

    return(round(cinza))
  }, USE.NAMES = FALSE)

  return(cinza_values)
}

## ______________________________________________________________________________

# Calculo da distancia
distancia <- function(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, l) {
  A <- complex(real = Ax, imaginary = Ay)
  B <- complex(real = Bx, imaginary = By)
  C <- complex(real = Cx, imaginary = Cy)
  D <- complex(real = Dx, imaginary = Dy)

  dAC <- C - A
  dBC <- C - B
  dAD <- D - A
  dBD <- D - B

  k <- (dAC / dBC) / (dAD / dBD)

  d_c <- sqrt(k / (k - 1) * l^2)
  d <- sqrt(Re(d_c)^2 + Im(d_c)^2)

  return(d)
}

## ______________________________________________________________________________

rotacionar_pontos <- function(df, angulo) {
  df_rotacionado <- df
  df_rotacionado$x <- df$x * cos(angulo) - df$y * sin(angulo)
  df_rotacionado$y <- df$x * sin(angulo) + df$y * cos(angulo)
  return(df_rotacionado)
}
