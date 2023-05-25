IRG<-function(midS, scaleS, t, minMax=T){
  numerator<-exp((t+midS)/scaleS)
  demoninator<-(2*scaleS*exp((t+midS)/scaleS))+
    (scaleS*exp((2*t)/scaleS))+     
    (scaleS*exp((2*midS)/scaleS))
  irg<-numerator/demoninator

  irg<-(irg-min(irg, na.rm=T))/(max(irg, na.rm=T)-min(irg, na.rm=T))
  return(irg)
}
