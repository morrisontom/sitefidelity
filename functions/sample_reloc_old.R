# Function for sampling 1 location per day for a dataset -- can sample >1 by changing line 10 & 13

smp<-function(x) { 
  dat3<-NULL
  for (j in 1:length(unique(x$AID))) { 
    tmp<-x[x$AID==unique(x$AID)[j],]
    for (i in 1:length(unique(tmp$Date))) {
      del<-nrow(tmp[tmp$Date==unique(tmp$Date)[i],])
      if(del>1){
        dat3<-rbind(dat3,tmp[tmp$Date==unique(tmp$Date)[i],][sample(nrow(tmp[tmp$Date==unique(tmp$Date)[i],]),1),])
      } 
      else {
        dat3<-rbind(dat3,tmp[tmp$Date==unique(tmp$Date)[i],])
      }
    }
  }
  return(dat3)
}
