get_hexcolors_for_palettes <-
function(color.palette=c("custom1","wong","npg","jama","jco","lancet","nejm"),alpha=1){
  
  if(color.palette[1]=="custom1"){
    color.palette=c("#474A49","#92C147","#F79646","#8064A2","#11BCFF","#0F7BA0")
  }else{
    
    if(color.palette[1]=="wong"){
      
      #https://www.nature.com/articles/nmeth.1618
      #rgb(204,121,167,maxColorValue = 255) Output: #CC79A7
      color.palette=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
      
    }else{
      color.palette=color.palette[1]
      
      if(color.palette=="npg"){
        
        # color.palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
        #                "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
        
        colfunc=pal_npg_custom(alpha=alpha.col)
        color.palette=colfunc(10)
      }else{
        
        if(color.palette=="nejm"){
          
          colfunc=pal_nejm_custom(alpha=alpha.col)
          color.palette=colfunc(8)
          
        }else{
          if(color.palette=="jco"){
            
            colfunc=pal_jco_custom(alpha=alpha.col)
            color.palette=colfunc(10)
            
          }else{
            
            if(color.palette=="jama"){
              
              colfunc=pal_jama_custom(alpha=alpha.col)
              color.palette=colfunc(7)
              
            }else{
              if(color.palette=="lancet"){
                
                colfunc=pal_lancet_custom(alpha=alpha.col)
                color.palette=colfunc(9)
                
              }
              
            }
          }
          
        }
      }
    }
  }
  return(color.palette)
}
