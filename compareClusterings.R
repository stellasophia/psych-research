library(psych)
source("faktorensimulation.R")
source("cmdsolve.R")
source("Vergleichsverfahren.R")

compareClusterings <- function(NLmean,NLsd,phimean,phisd, comparing,toSimulate, nrep=1, addError=F) {
  #  jpeg(paste("/home/andreas/Desktop/Bachelorarbeit/",NL,phi,comparing,".jpg"), width = 600, height = 400)
  
  completecorresult <- 0
  completecorresult <- 0
  completecorcorresult <- 0
  averagecorresult <- 0
  averagecorcorresult <- 0
  kmeansresult <- 0
  kmeansCorResult<-0
  kmeansNewResult<-0
  completecorcornometricresult <-0
  completecornometricresult <-0
  averagecornometricresult <- 0
  averagecorcornometricresult <- 0
  completeCorCorCorResult <- 0
  averageCorCorCorResult <- 0
  kmeanscorcorResult <- 0
  cmdsolveCorResult <- 0
  kMeansOnDistancesCorResult <- 0
  
  
  
  results <- c()
  sumresults <- c(rep(0,length(toSimulate)))
  for(i in 1:nrep) {
    print(paste(NLmean,NLsd,phimean,phisd))
   corM <- setCorrelationMatrix(NLmean,NLsd,phimean,phisd, addError)
 #   cat("corMatrix set! ", NLmean,NLsd,phimean,phisd)
#    completecor <- completeCor(corM)
#    completecorcor <- completeCorCor(corM)
 #   averagecor <- averageCor(corM)
#    averagecorcor <- averageCorCor(corM)
#    kmeans <- cmdsolve(corM)
    
    
    results <- c()
    resultnames <- toSimulate
   for(sim in toSimulate) {
     print(sim)
     if(sim=="averagecor") {
         averagecor <- averageCor(corM)
        
         averagecorresult <- averagecorresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
       results <-append(results,averagecorresult/nrep)
     } else if(sim=="averagecorcor") {
     
         averagecorcor <- averageCorCor(corM)
         averagecorcorresult <- averagecorcorresult  + sum(vergleich(zuordnung.ges,averagecorcor, comparing))

       results <-append(results,averagecorcorresult/nrep)
     } else if(sim=="completecor") {
      
         completecor <- completeCor(corM)
         completecorresult <- completecorresult  + sum(vergleich(zuordnung.ges,completecor, comparing))
       results <-append(results,completecorresult/nrep)
     } else if(sim=="completecorcor") {
      
         completecorcor <- completeCorCor(corM)
         completecorcorresult <- completecorcorresult  + sum(vergleich(zuordnung.ges,completecorcor, comparing))
       results <-append(results,completecorcorresult/nrep)
     }  else if(sim=="kmeansmds") {   
         kmeans <- cmdsolve(corM)
         kmeansresult <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult/nrep)
     } else if(sim=="Dim1") {   
       kmeans <- cmdsolve(corM,dim=1)
       kmeansresult1 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult1/nrep)
     }  else if(sim=="Dim2") {   
       kmeans <- cmdsolve(corM,dim=2)
       kmeansresult2 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult2/nrep)
     }  else if(sim=="Dim3") {   
       kmeans <- cmdsolve(corM,dim=3)
       kmeansresult3 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult3/nrep)
     }  else if(sim=="Dim4") {   
       kmeans <- cmdsolve(corM,dim=4)
       kmeansresult4 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult4/nrep)
     } else if(sim=="Dim10") {   
       kmeans <- cmdsolve(corM,dim=10)
       kmeansresult10 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult10/nrep)
     } else if(sim=="Dim20") {   
       kmeans <- cmdsolve(corM,dim=20)
       kmeansresult20 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult20/nrep)
     } else if(sim=="Dim30") {   
       kmeans <- cmdsolve(corM,dim=30)
       kmeansresult30 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
       results <-append(results,kmeansresult30/nrep)
     }  else if(sim=="kmeanskoord") {
         kmeansNew <- kMeansOnDistances(corM)
         kmeansNewResult <- kmeansNewResult  + sum(vergleich(zuordnung.ges,kmeansNew, comparing))
       results <-append(results,kmeansNewResult/nrep)
     }  else if(sim=="kmeanscor") {
         kmeansCor <- kmeansCor(corM)
         kmeansCorResult <- kmeansCorResult  + sum(vergleich(zuordnung.ges,kmeansCor, comparing))
       results <-append(results,kmeansCorResult/nrep)
     } else if(sim=="completecorcorcor") {
       completeCorCorCor <- completeCorCorCor(corM)
       completeCorCorCorResult <- completeCorCorCorResult  + sum(vergleich(zuordnung.ges,completeCorCorCor, comparing))
       results <-append(results,completeCorCorCorResult/nrep)
     } else if(sim=="averagecorcorcor") {
       averageCorCorCor <- averageCorCorCor(corM)
       averageCorCorCorResult <- averageCorCorCorResult  + sum(vergleich(zuordnung.ges,averageCorCorCor, comparing))
       results <-append(results,averageCorCorCorResult/nrep)
     } else if(sim=="kmeanscorcor") {
       print("kmeanscorcor")
       kmeanscorcor<- kmeansCorCor(corM)
       print(kmeanscorcor)
       kmeanscorcorResult <- kmeanscorcorResult  + sum(vergleich(zuordnung.ges,kmeanscorcor, comparing))
       print(print(kmeanscorcorResult))
       results <-append(results,kmeanscorcorResult/nrep)
     } else if(sim=="kmeansmdscor") {
       cmdsolveCor <- cmdsolveCor(corM)
       cmdsolveCorResult <-  cmdsolveCorResult   + sum(vergleich(zuordnung.ges, cmdsolveCor, comparing))
       results <-append(results,cmdsolveCorResult/nrep)
     }  else if(sim=="kMeansondistancescor ") {
       kMeansOnDistancesCor  <- kMeansOnDistancesCor(corM)
       kMeansOnDistancesCorResult <-  kMeansOnDistancesCorResult   + sum(vergleich(zuordnung.ges,kMeansOnDistancesCor, comparing))
       results <-append(results,kMeansOnDistancesCorResult/nrep)
     }  else if(sim=="averagecornom") {
        averagecor <- averageCorNoMetric(corM)
        averagecornometricresult <- averagecornometricresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
     results <-append(results,averagecornometricresult/nrep)
     } else if(sim=="averageccnom") {
         averagecorcor <- averageCorCorNoMetric(corM)
         averagecorcornometricresult <- averagecorcornometricresult + sum(vergleich(zuordnung.ges,averagecorcor, comparing))
       results <-append(results,averagecorcornometricresult/nrep)
     } else if(sim=="completecornom") {
         completecor <- completeCorNoMetric(corM)
         completecornometricresult <- completecornometricresult  + sum(vergleich(zuordnung.ges,completecor, comparing))
       results <-append(results,completecornometricresult/nrep)
     } else if(sim=="completeccnom") {
         completecorcor <- completeCorCorNoMetric(corM)
         completecorcornometricresult <- completecorcornometricresult  + sum(vergleich(zuordnung.ges,completecorcor, comparing))
       results <-append(results,completecorcornometricresult /nrep)
     } 
   }
   print(results)
    sumresults <- sumresults + results
  }
  
    mainText <- "Rand-Index"
    if(comparing==2) {
      mainText <- "Meila-Heckermann-Measure"
    }
    
    
    
   # completecorresult <- completecorresult + sum(vergleich(zuordnung.ges,completecor, comparing))
  #  completecorcorresult <- completecorcorresult +  sum(vergleich(zuordnung.ges,completecorcor,comparing))
  #  averagecorresult <- averagecorresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
  #  averagecorcorresult <- averagecorcorresult + sum(vergleich(zuordnung.ges,averagecorcor, comparing))
  #  kmeansresult <- kmeansresult + sum(vergleich(zuordnung.ges,kmeans, comparing))

  #results <<- c(completecorresult/nobs, completecorcorresult/nobs,averagecorresult/nobs,averagecorcorresult/nobs,  kmeansresult/nobs)
  
  print(results)
  print(toSimulate)
  names(results) <- toSimulate

  print(results)
  
  results
  
 # barplot(results, beside=T,main="Erkennen von Faktorenstrukturen",
#          xlab="Clusterverfahren bei", ylab=mainText, ylim=c(0,1.0), col=gray.colors(length(results)), sub=optionstext,cex.names=1.3,font.names=2,cex.axis=,cex.main=1.5,cex.lab=1.5,cex.sub=1.5)
  
 # legend(x="topleft",col=c("green","red"), legend= resultnames, fill=gray.colors(5), inset=c(0.01, 0), ce=0.4, bty="n")
}





getComparison <- function(toSimulates, names,distribution,addError=F) {
  counter <- 1
  error <- ""
  if(addError) {
    error <- "Fehler"
  }
  folder <- "/home/andreas/Desktop/Bachelorarbeit2/faktorbasis/"
  for(toSimulate in toSimulates) {
    jpeg(paste(folder,names[counter],"/NLzunehmend",error,".jpeg",sep=""), width = 1600, height = 800)
    
   par(mfrow=c(2,2), oma=c(0,0,0,0))
    
    compareClusterings(0,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.1,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.15,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.2,0,0,0,1,toSimulate, addError=addError)
    
    dev.off()
    
    jpeg(paste(folder,names[counter],"/biggerNLzunehmend",error,".jpeg",sep=""), width = 1600, height = 800)
    
    par(mfrow=c(2,2), oma=c(0,0,0,0))
    
    compareClusterings(0.225,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.25,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.275,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.3,0,0,0,1,toSimulate, addError=addError)

    
    legend(x="topright", legend= toSimulate, fill=gray.colors(length(toSimulate)),inset=c(-0.15,0),  bty="n",cex=1.2)
    
  dev.off()
    
    
   # jpeg(paste(folder,"/biggerNLzunehmendWithVar",error,".jpeg",sep=""), width = 1600, height = 800)
    
  #  par(mfrow=c(2,2), oma=c(0,0,0,0))
    
#    compareClusterings(0.2,0.05,0,0,1,toSimulate, addError=addError)
    
 #   compareClusterings(0.2,0.1,0,0,1,toSimulate, addError=addError)
    
#  compareClusterings(0.2,0.15,0,0,1,toSimulate, addError=addError)
    
#    compareClusterings(0.2,0.2,0,0,1,toSimulate, addError=addError)
    
#    dev.off()
    
 #   jpeg(paste(folder,names[counter],"/biggestNLzunehmend",error,".jpeg",sep=""), width = 1200, height = 800)
    
  #  par(mfrow=c(2,2), oma=c(0,0,0,0))
    
##   compareClusterings(0.3235,0,0,0,1,toSimulate, addError=addError)
    
#    compareClusterings(0.35,0,0,0,1,toSimulate, addError=addError)
    
 #   compareClusterings(0.375,0,0,0,1,toSimulate, addError=addError)
    
  #  compareClusterings(0.4,0,0,0,1,toSimulate, addError=addError)
    
  ##  dev.off()
    
    
  jpeg(paste(folder,names[counter],"/Korrzunehmend",error,".jpeg",sep=""), width = 1600, height = 800)
    par(mfrow=c(2,2), oma=c(0,0,0,0))
    
    compareClusterings(0,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0,0,0.2,0,1,toSimulate,addError=addError)
    
  compareClusterings(0,0,0.3,0,1,toSimulate, addError=addError)
    compareClusterings(0,0,0.4,0,1,toSimulate, addError=addError)
  
    dev.off()
    
    
 jpeg(paste(folder,names[counter],"/biggerKorrzunehmend",error,".jpeg",sep=""), width = 1600, height = 800)
    par(mfrow=c(2,2), oma=c(0,0,0,0))
    
    compareClusterings(0,0,0.5,0,1,toSimulate, addError=addError)
    
    compareClusterings(0,0,0.6,0,1,toSimulate, addError=addError)
    
    compareClusterings(0,0,0.7,0,1,toSimulate, addError=addError)
    
    compareClusterings(0,0,0.8,0,1,toSimulate, addError=addError)
    
    dev.off()
    
    
    jpeg(paste(folder,names[counter],"/KorrundNLzunehmend",error,".jpeg",sep=""), width = 1600, height = 800)
   par(mfrow=c(2,2), oma=c(0,0,0,0))
    
    compareClusterings(0,0,0,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.1,0,0.1,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.15,0,0.15,0,1,toSimulate, addError=addError)
    
    compareClusterings(0.175,0,0.175,0,1,toSimulate, addError=addError)
    
    dev.off()
    
    
    if(distribution) {
      jpeg(paste(folder,names[counter],"/NLnormalverteiltt",error,".jpeg",sep=""), width = 1600, height = 800)
      par(mfrow=c(2,2), oma=c(0,0,0,0))
      
     compareClusterings(0.01,0,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.01,0.1,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.01,0.15,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.01,0.2,0,0,1,toSimulate,nrep=10, addError=addError)
      
      dev.off()
      
      
      jpeg(paste(folder,names[counter],"/NLnormalverteilttbigNL",error,".jpeg",sep=""), width = 1600, height = 800)
      par(mfrow=c(2,2), oma=c(0,0,0,0))
      
      compareClusterings(0.2,0,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.2,0.05,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.2,0.075,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0.2,0.1,0,0,1,toSimulate,nrep=10, addError=addError)
      
      dev.off()
      
      
      jpeg(paste(folder,names[counter],"/Korrnormalverteilt",error,".jpeg",sep=""), width = 1600, height = 800)
      par(mfrow=c(2,2), oma=c(0,0,0,0))
      
      compareClusterings(0,0,0,0,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0,0,0,0.1,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0,0,0,0.3,1,toSimulate,nrep=10, addError=addError)
      
      compareClusterings(0,0,0,0.5,1,toSimulate,nrep=10, addError=addError)
      
      dev.off()
    }
    
    counter <- counter+1
  }
}


