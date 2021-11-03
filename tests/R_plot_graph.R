R_plot_graph = function(FEM, ...){
  nodes <- FEM$FEMbasis$mesh$nodes
  edges <- as.vector(t(FEM$FEMbasis$mesh$edges))
  coeff <- FEM$coeff
  
  p <- jet.col(n=128,alpha=0.8)
  palette(p)
  ncolor <- length(p)
  nplots <- ncol(coeff)
  
  for (i in 1:nplots){
    if (i > 1)
      readline("Press any key for the next plot...")
    
    open3d()
    axes3d()
    rgl.pop("lights")
    light3d(specular="black")
    
    col <- coeff[edges,i]
    col <- (ncolor-1)*(col-min(col))/diff(range(col))+1
    col <- p[col]
    
    rgl.lines(nodes[edges,1], nodes[edges,2],rep(0,dim(nodes)[1]),
              color = col,lwd=2,...)
    rgl.viewpoint(0,0)
    }
}


R_plot_graph.ggplot2<-function( FEM ){
  
  
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    coef=append(coef, rep((FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]/2),times=2) )  
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.) + 
    geom_line(aes(color=coef),size=1) +
    scale_color_gradientn(colours=p ) + 
    labs(x="",y="",color="") + 
    theme(panel.background = element_rect(fill="white",color="gray50"))
  
  
}

