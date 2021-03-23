pal_jco_custom <-
function (palette = c("default"), alpha = 1) 
{
  palette = match.arg(palette)
  if (alpha > 1L | alpha <= 0L) 
    stop("alpha must be in (0, 1]")
  raw_cols = c(
    "Lochmara" = "#0073C2", "Corn" = "#EFC000",
    "Gray" = "#868686", "ChestnutRose" = "#CD534C",
    "Danube" = "#7AA6DC", "RegalBlue" = "#003C67",
    "Olive" = "#8F7700", "MineShaft" = "#3B3B3B",
    "WellRead" = "#A73030", "KashmirBlue" = "#4A6990"
  )
  
  raw_cols_rgb = col2rgb(raw_cols)
  alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
                   raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), 
                   maxColorValue = 255L)
  manual_pal(unname(alpha_cols))
}
