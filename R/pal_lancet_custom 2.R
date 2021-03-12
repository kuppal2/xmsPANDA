pal_lancet_custom <-
function (palette = c("lanonc"), alpha = 1) 
{
  palette = match.arg(palette)
  if (alpha > 1L | alpha <= 0L) 
    stop("alpha must be in (0, 1]")
  raw_cols = c(
    "CongressBlue" = "#00468B", "Red" = "#ED0000",
    "Apple" = "#42B540", "BondiBlue" = "#0099B4",
    "TrendyPink" = "#925E9F", "MonaLisa" = "#FDAF91",
    "Carmine" = "#AD002A", "Edward" = "#ADB6B6",
    "CodGray" = "#1B1919"
  )
  
  raw_cols_rgb = col2rgb(raw_cols)
  alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
                   raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), 
                   maxColorValue = 255L)
  manual_pal(unname(alpha_cols))
}
