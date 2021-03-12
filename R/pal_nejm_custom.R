pal_nejm_custom <-
function (palette = c("default"), alpha = 1) 
{
  palette = match.arg(palette)
  if (alpha > 1L | alpha <= 0L) 
    stop("alpha must be in (0, 1]")
  raw_cols = c(
    "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
    "Zest" = "#E18727", "Eucalyptus" = "#20854E",
    "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
    "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
  )
  raw_cols_rgb = col2rgb(raw_cols)
  alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
                   raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), 
                   maxColorValue = 255L)
  manual_pal(unname(alpha_cols))
}
