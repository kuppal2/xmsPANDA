pal_npg_custom <-
function (palette = c("nrc"), alpha = 1) 
{
  
  
  palette = match.arg(palette)
  if (alpha > 1L | alpha <= 0L) 
    stop("alpha must be in (0, 1]")
  raw_cols = c(
    "Cinnabar" = "#E64B35", "Shakespeare" = "#4DBBD5",
    "PersianGreen" = "#00A087", "Chambray" = "#3C5488",
    "Apricot" = "#F39B7F", "WildBlueYonder" = "#8491B4",
    "MonteCarlo" = "#91D1C2", "Monza" = "#DC0000",
    "RomanCoffee" = "#7E6148", "Sandrift" = "#B09C85"
  )
  raw_cols_rgb = col2rgb(raw_cols)
  alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
                   raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), 
                   maxColorValue = 255L)
  manual_pal(unname(alpha_cols))
}
