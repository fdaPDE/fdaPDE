R_plot_graph = function(FEM, ...){
  
  nodes <- FEM$FEMbasis$mesh$nodes
  if(FEM$FEMbasis$order==1){
        edges <- as.vector(t(FEM$FEMbasis$mesh$edges))
  }else{
        edges <- cbind(FEM$FEMbasis$mesh$edges[,1],
                      FEM$FEMbasis$mesh$edges[,3],
                      FEM$FEMbasis$mesh$edges[,3],
                      FEM$FEMbasis$mesh$edges[,2])
        edges <- as.vector(t(edges))
  }
  coeff <- FEM$coeff
  
  p <- jet.col(n=128,alpha=0.8)
  palette(p)
  ncolor <- length(p)
  nplots <- ncol(coeff)
  
  for (i in 1:nplots){
    if (i > 1)
      readline("Press any key for the next plot...")
    
    open3d(windowRect=c(10,10,500,500))
    layout3d(matrix(1:2,1,2),c(0.8,0.2),1 ) 
    axes3d()
    rgl.pop("lights")
    light3d(specular="black")
    
    col <- coeff[edges,i]
    col <- (ncolor-1)*(col-min(col))/diff(range(col))+1
    col <- p[col]
    
    rgl.lines(nodes[edges,1], nodes[edges,2],rep(0,dim(nodes)[1]),
              color = col,lwd=1.75,...)
    rgl.viewpoint(0,0,zoom=0.75)  
    next3d()
    bgplot3d({
      plot.new()
      plotrix::color.legend(0.1,0.1,0.9,0.9,
                   rect.col=p,
                   legend = seq(floor(min(FEM$coeff[,i])),
                                max(FEM$coeff[,i]), 
                                0.2),
                   gradient="y",cex=1.)
    })
  
  }
  
}

R_plot_graph.ggplot2<-function( FEM ){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  if(FEM$FEMbasis$order==1){
    coef.f = function(i){(FEM$coeff[mesh$edges[i,1]] + FEM$coeff[mesh$edges[i,2]])/2}
  }else{
    coef.f = function(i){(FEM$coeff[mesh$edges[i,1]] + FEM$coeff[mesh$edges[i,2]] + FEM$coeff[mesh$edges[i,3]])/3 }
  }
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
   # coef=append(coef, rep( (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]])/2 ,times=2) )  
    coef = append(coef, rep(coef.f(e), times=2))
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

R_plot_mesh.1.5D<-function(mesh,alpha=0.){
  if(mesh$order==1){
   x=vector(mode="double")
   y=vector(mode="double")
   grp.nodes=vector(mode="integer")
   num_edges=dim(mesh$edges)[1]
   for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    grp.nodes = append(grp.nodes, rep(e,times=2))
   }
  
   data=data.frame(x,y,grp.nodes)
   ggplot(data=data, aes(x=x,y=y,group=grp.nodes) ) + 
     geom_point(alpha=alpha) + 
     geom_line() + 
     labs(x="",y="",color="") + 
     theme(panel.background = element_rect(fill="white",color="gray50"))
  }else{
    x=vector(mode="double")
    y=vector(mode="double")
    x.middle=vector(mode="double")
    y.middle=vector(mode="double")
    
    grp.nodes=vector(mode="integer")
    num_edges=dim(mesh$edges)[1]
    for(e in 1:num_edges){
      x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
      y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
      grp.nodes = append(grp.nodes, rep(e,times=2))
      
      x.middle = append(x.middle, mesh$nodes[ mesh$edges[e,3],1] )
      y.middle = append(y.middle, mesh$nodes[ mesh$edges[e,3],2] )
      
      }
    
    data=data.frame(x,y,grp.nodes,x.middle,y.middle)
    ggplot(data=data, aes(x=x,y=y,group=grp.nodes) ) + 
      geom_point(alpha=alpha) + 
      geom_line() + 
      geom_point(aes(x=x.middle,y=y.middle),alpha=alpha,color="red") +
      labs(x="",y="",color="") + 
      theme(panel.background = element_rect(fill="white",color="gray50"))
    
  }
  
}


plot.mesh.1.5D.new<-function(x,plot.nodes=FALSE, ...)
{
  if(plot.nodes)
    type="p"
  else
    type="n"
  
  if( x$order == 1 ){
    
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type=type, ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
             x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  }
  else{
    
    plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n",type=type, ...)
    segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
             x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
    points(x$nodes[x$edges[,3],1], x$nodes[x$edges[,3],2],col="red",type=type,...)
  }
}

plot.mesh.2.5D.new<-function(x,plot.nodes=FALSE,...){
  if(plot.nodes) 
    alpha = 1.0 
  else 
    alpha = 0.
  
  nodes <- x$nodes
  edges <- as.vector(t(x$edges))
  
  open3d()
  axes3d()
  rgl.pop("lights")
  light3d(specular="black")
  
  rgl.points(nodes[,1], nodes[,2], nodes[,3], col="black",alpha=alpha, ...)
  rgl.lines(nodes[edges,1], nodes[edges,2], nodes[edges,3], col="black",...)
  
  aspect3d("iso")
  rgl.viewpoint(0,-45)
  
}

plot.mesh.2D.new<-function(x,plot.nodes=FALSE, ...)
{
  if(plot.nodes)
    type="p"
  else
    type="n"
  
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type=type, ...)
  segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  segments(x$nodes[x$segments[,1],1], x$nodes[x$segments[,1],2],
           x$nodes[x$segments[,2],1], x$nodes[x$segments[,2],2], col="red", ...)
}

plot.mesh.3D.new<-function(x,plot.nodes=FALSE,...){
  if(plot.nodes)
    alpha=1.0
  else
    alpha=0.
  
  nodes <- x$nodes
  faces <- as.vector(t(x$faces[x$facesmarkers,]))
  
  aux_mesh <- create.mesh.2.5D(nodes=nodes, triangles=x$faces[x$facesmarkers,], order=1)
  edges <- as.vector(t(aux_mesh$edges))
  
  open3d()
  axes3d()
  rgl.pop("lights")
  light3d(specular="black")
  
  rgl.triangles(nodes[faces,1],nodes[faces,2],nodes[faces,3],col="white",...)
  rgl.points(nodes[,1], nodes[,2], nodes[,3], col="black",alpha=alpha, ...)
  rgl.lines(nodes[edges,1], nodes[edges,2], nodes[edges,3], col="black",...)
  
  aspect3d("iso")
  rgl.viewpoint(0,-45)
  
}
