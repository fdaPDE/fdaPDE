library(ggplot2)
R_plot_graph.ggplot2<-function(FEM){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    coef=append(coef, rep( (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]) /2,times=2) )  
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  MyTheme <- theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
    title = element_text(size=20),
    legend.text = element_text(size=14),
    legend.key.size = unit(1,"cm") 
  )
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.01) + 
    geom_line(aes(color=coef)) +
    scale_color_gradientn(colours=p) + 
    labs(x="",y="",color="") + 
    theme(panel.background = element_rect(fill="white",color="gray50"))+
    MyTheme
  
  
}

R_plot_graph.a.sym.ggplot2<-function(FEM){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    if( abs(FEM$coeff[mesh$edges[e,1]])<= 10*.Machine$double.eps 
        ||
        abs(FEM$coeff[mesh$edges[e,2]]) <= 10*.Machine$double.eps){
      coef = append(coef, c(0.0, 0.0) )
    }else{
      coef= append(coef, rep( (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]) /2,times=2) )  
    }
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  MyTheme <- theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
    title = element_text(size=20),
    legend.text = element_text(size=14),
    legend.key.size = unit(1,"cm") 
  )
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.01) + 
    geom_line(aes(color=coef)) +
    scale_color_gradientn(colours=p) + 
    labs(x="",y="",color="") + 
    theme(panel.background = element_rect(fill="white",color="gray50"))+
    MyTheme
  
  
}


R_plot_graph.ggplot2.2<-function(FEM, color.min=min(FEM$coeff), color.max=max(FEM$coeff)){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    coef=append(coef, rep((FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]])/2,times=2) )  
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.01) + 
    geom_line(aes(color=coef)) +
    scale_color_gradientn(colours=p, limits = c(color.min, color.max)) + 
    labs(x="",y="",color="") + 
    theme(panel.background = element_rect(fill="white",color="gray50"))
  
  
}
