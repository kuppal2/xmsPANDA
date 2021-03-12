pal_jama_custom <-
function (palette = c("default"), alpha = 1) 
{
  palette = match.arg(palette)
  if (alpha > 1L | alpha <= 0L) 
    stop("alpha must be in (0, 1]")
  raw_cols = c(
    "Limed Spruce" = "#374E55", "Anzac" = "#DF8F44",
    "Cerulean" = "#00A1D5", "Apple Blossom" = "#B24745",
    "Acapulco" = "#79AF97", "Kimberly" = "#6A6599",
    "Makara" = "#80796B"
  )
  
  raw_cols_rgb = col2rgb(raw_cols)
  alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
                   raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), 
                   maxColorValue = 255L)
  manual_pal(unname(alpha_cols))
}
