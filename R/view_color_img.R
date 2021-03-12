view_color_img <-
function(color.palette=c("custom1","npg","jama","jco","lancet","nejm"),alpha.col=1){
  
  color.palette=get_hexcolors_for_palettes(color.palette[1],alpha=alpha.col)
  
  show_col(color.palette,ncol=1)
  
}
