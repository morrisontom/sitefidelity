

#### SOURCE THIS CODE FIRST!
make_clust <- function(data=data, freq_cutoff=1, use_weights=TRUE, name="WB_SE", plot_it=FALSE){
  if(all(names(data)%in%c("from","to","freq"))==FALSE) 
    stop("Your dataframe must contain 3 columns with the names from, to, and freq")
  
  data <- data[data$freq >= freq_cutoff,]
  g <- graph.data.frame(data, directed=F)
  
  if(use_weights==TRUE){
    E(g)$weight <- data$freq
  }
  bc <- edge.betweenness.community(g)
  
  if(plot_it==TRUE){
    plot(bc,g)
  }
  return(data.frame(name=name, AID=names(membership(bc)), cluster=as.numeric(membership(bc))))
}