
## Function to identify and return 2 states (sedentary and migratory) 
## for each GPS fix using moveHMM. This ignores transitions probabilities
## and relationships with covariates

# state decoding only seems to work on 1 individual at a time 
# 1=sedentary 
# 2=migratory

hmm2state <- function(x){
  
  library(moveHMM)
  
  # convert everything to km
  x$x <- x$x/1000
  x$y <- x$y/1000
  x<-x[,names(x) %in% c("x","y")]
  
  #prepare data
  suppressWarnings( #suppress warning as this will always give a warning about covariates
    data <- prepData(x,type="UTM",coordNames=c("x","y"))
  )
  ## initial parameters for gamma and von Mises distributions
  mu0 <- c(median(data$step,na.rm=T)-median(data$step,na.rm=T)/2,median(data$step,na.rm=T)) # step mean (two parameters: one for each state)
  sigma0 <- c(sd(data$step,na.rm=T)/2,median(data$step,na.rm=T)) # step SD
  stepPar0 <- c(mu0,sigma0)
  
  if(length(data$step[data$step==0])>1){
    zeromass0 <- c(0.01,0.01)
    stepPar0 <- c(mu0,sigma0,zeromass0)
  }
  angleMean0 <- c(0,0) # angle mean
  kappa0 <- c(0.1,1) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  
  m <- fitHMM(data=data,
              nbStates=2,
              stepPar0=stepPar0,
              anglePar0=anglePar0,
              angleDist = "vm",
              formula=~1)
  
  states <- viterbi(m)
  return(states)
}

