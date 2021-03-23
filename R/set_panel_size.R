set_panel_size <-
function (p = NULL, g = ggplot2::ggplotGrob(p), file = NULL, 
          margin = unit(1, "mm"), width = unit(4, "cm"), 
          height = unit(4, "cm")) 
{
  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  g$widths[panel_index_w] <- rep(width, nw)
  g$heights[panel_index_h] <- rep(height, nh)
  if (!is.null(file)) {
    ggplot2::ggsave(file, g, width = grid::convertWidth(sum(g$widths) + 
                                                          margin, unitTo = "in", valueOnly = TRUE), height = grid::convertHeight(sum(g$heights) + 
                                                                                                                                   margin, unitTo = "in", valueOnly = TRUE))
  }
  g
}
