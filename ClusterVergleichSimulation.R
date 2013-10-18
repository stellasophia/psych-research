
simulate <- function(cor.sp,nobs,nrep,k,method,wholecut,compareMethod,data,...) {
  
  print(method)
  v1 <- 1:nrep
  v2 <- 1:nrep
  w1 <- list()
  w2 <- list()
  l1  <-  list()
  l2  <-  list()
  
  cut.sp1 <- c()
  
  vergleich1 <- 0
  n11 <- 0
  n00 <- 0
    if(method=="completecor")  {
      cut.sp1 <- completeCor(cor.sp,k)
    } else if(method=="completecorcor") {
      cut.sp1 <- completeCorCor(cor.sp,k)
    } else if(method=="averagecor") {
      cut.sp1 <- averageCor(cor.sp,k)
    } else if(method=="averagecorcor") {
      cut.sp1 <- averageCorCor(cor.sp,k)
    } else if(method=="kmeansmds") {
      cut.sp1 <- cmdsolve(cor.sp,k)
    }else if(method=="Dim1") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=1)
    }else if(method=="Dim2") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=2)
    }else if(method=="Dim3") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=3)
    }  else if(method=="Dim4") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=4)
    }  else if(method=="Dim10") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=10)
    }  else if(method=="Dim20") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=20)
    }   else if(method=="Dim30") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=30)
    } else if(method=="kmeanskoord") {
      cut.sp1 <- kMeansOnDistances(cor.sp,k)
    } else if(method=="kmeanscor") {
      cut.sp1 <- kmeansCor(cor.sp,k)
    } else if(method=="completecornom")  {
      cut.sp1 <- completeCorNoMetric(cor.sp,k)
    } else if(method=="completecorcornom") {
      cut.sp1 <- completeCorCorNoMetric(cor.sp,k)
    } else if(method=="averagecornom") {
      cut.sp1 <- averageCorNoMetric(cor.sp,k)
    } else if(method=="averagecorcornom") {
      cut.sp1 <- averageCorCorNoMetric(cor.sp,k)
    }  else if(method=="totalclust") {
      cut.sp1  <- kmeans(t(data),centers=5,nstart=100)$cluster
    } else if(method=="averagecorcorcor") {
      cut.sp1  <- averageCorCorCor(cor.sp,k)
    } else if(method=="completecorcorcor") {
      cut.sp1  <- completeCorCorCor(cor.sp,k)
    }  else if(method=="kmeansmdscor") {
      cut.sp1  <- cmdsolveCor(cor.sp)
    } else if(method=="kmeanscorcor") {
      cut.sp1  <- kmeansCorCor(cor.sp,k)
    } else if(method=="kmeansneucor") {
      cut.sp1  <- kMeansOnDistancesCor(cor.sp,k)
    } else if(method=="faclust") {
      cut.sp1 <- fclustering(cor.sp,k)
    }

  print(cut.sp1)
  print(wholecut)
    n00 <- n00 + vergleich(wholecut,cut.sp1, compareMethod)[1]
    n11 <- n11 + vergleich(wholecut,cut.sp1, compareMethod)[2]
print(n00)

  result <- c(n00/nrep, n11/nrep)
  print(result)
  result
}


simulateClusterSamplesComparison <- function(facs,toSimulate,compareWith, allnobs, nrep,numbercluster,compareMethod ) {
  simulationresults <- c()
  counter <- 1
for(nobs in allnobs) {
 
  #numbercluster4
sumvalues <- c(rep(0,2*length(toSimulate)))
resultnames <- c()
resultvalues <- c()
  for (i in 1:nrep){
    
    daten.sp <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    resultvalues <- c()

    cor.sp <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
   for(sim in toSimulate) {
    newresult <-  simulate(cor.sp,nobs,nrep, numbercluster,method=sim,compareWith[[sim]], compareMethod,daten.sp)
    print(newresult)
    resultvalues <- append(resultvalues,newresult)
    resultnames <- sim
  }
    sumvalues  <- sumvalues + resultvalues
  }
  resultvalues <- sumvalues
  print(resultvalues)
  #resultnames <- c("Average,  Korrelation", "Complete, Korrelation",
   #                "Averag, Korrelation der Korrelation", "Complete, Korrelation der Korrelation", "kmeans cmd")
  
  simulationinfo <- list()
  
  simulationinfo[["resultvalues"]] <- resultvalues
  
  simulationinfo[["resultnames"]] <- resultnames
  
  simulationinfo[["samplesize"]] <- nobs
  
  simulationinfo[["nreps"]] <- nrep
  
#  simulationinfo[["numbercluster"]] <- numbercluster
  
  simulationresults[[counter]] <-  simulationinfo

  counter <- counter + 1
}

  print(simulationresults)

##bringt diese Daten in ein richtiges Format für die grafische Darstellung
dataasvector11 <- c()
dataasvector00 <- c()
samplesizes <- c()


for(i in 1:length(simulationresults)) {
  print(simulationresults[[i]][["resultvalues"]])
  for(j in 1:length(simulationresults[[i]][["resultvalues"]])) { 
    dataasvector00 <- as.vector(na.omit(append(dataasvector00, simulationresults[[i]][["resultvalues"]][2*j-1])))
    samplesizes <- append(samplesizes,paste(simulationresults[[i]][["samplesize"]]) )
  }
}


for(i in 1:length(simulationresults)) {
  print(simulationresults[[i]][["resultvalues"]])
  for(j in 1:length(simulationresults[[i]][["resultvalues"]])) { 
    dataasvector11 <- as.vector(na.omit(append(dataasvector11, simulationresults[[i]][["resultvalues"]][2*j])))
    samplesizes <- append(samplesizes,paste(simulationresults[[i]][["samplesize"]]) )
  }
}


dataasmatrix00 <- matrix(c(dataasvector00), ncol=length(toSimulate), byrow=T)
dataasmatrix11 <- matrix(c(dataasvector11), ncol=length(toSimulate), byrow=T)
dataasmatrixall <- matrix(c(dataasvector00+dataasvector11), ncol=length(toSimulate), byrow=T)


#relevant für den Rand-Index in dem nur die wahr positiven oder wahr negativen
#und für den Rand-Index Spezifität und Sensitivität
#ansonsten immer dataasmatrixall verwenden 
dataasmatrix <-  dataasmatrixall
mainText <- "Rand-Index mit Clusteranzahl=5"
#dataasmatrix <- dataasmatrix11
#mainText <- "Rand-Index Sensitivität"
#dataasmatrix <- dataasmatrix00
#mainText <- "Rand-Index Spezifität"
if(compareMethod==2) {
  mainText <- "Meila-Heckerman-Measure"
} 


colnames(dataasmatrix) <- toSimulate
rownames(dataasmatrix) <- allnobs

print(t(dataasmatrix))

#dies ist die graphische Darstellung des Ergebnisses
par(mfrow=c(1,1), oma=c(0,0,0,0),xpd=TRUE,mar=c(5.1, 4.1, 4.1, 8.1))

barplot(as.table(t(dataasmatrix)), beside=T,main="Strukturerhaltungsanalyse von Klassifizierungsalgorithmen",
        xlab="Stichprobengrössen", ylab=mainText, ylim=c(0,1.0), ,cex.names=1.5,font.names=2,cex.axis=1.5,cex.main=1.5,cex.lab=1.5)

legend(x="topright", legend= toSimulate, fill=gray.colors(length(toSimulate)),inset=c(-0.1,0),  bty="n",cex=1.2)
}





#Anzahl der Wiederholungen der Simulation
nrep <- 20
#die Stichprobengrößen
allnobs <- c(100,150,250,500,1000)

#die gewünschten Clusteranzahlen. Wenn Wert größer als 2 wird es als die Clusteranzahl interpretiert,
#bei Werten kleiner als 2 als der Cutoff-wert h, in diesem Fall wird die Clusteranzahl für
#den Stichprobendatensatz dann dynamisch bestimmt. Bei der Faktoranalyse muss dafür variableClusterNumber
#auf True gesetzt werden
#average Korrelation
k <- 5
#complete Korrelation


#welcher Vergleichsalgorithmus gewählt wird:
#1: normaler Rand-Index
#2: Meila-Heckermann-Measure
#3: Rand Index für Spezifität und Sensitivität
compareMethod <- 2


corM <- cor(facs, use="pairwise.complete.obs", method="pearson")
#corcorM <- cor(corM, use="pairwise.complete.obs", method="pearson")

averagecor <- averageCor(corM,k) 
averagecorcor <- averageCorCor(corM,k) 
completecor <- completeCor(corM,k)
completecorcor <- completeCorCor(corM,k)
kmeansmds <- cmdsolve(corM,k)
kmeansmds1 <- cmdsolve(corM,k,dim=1)
kmeansmds2 <- cmdsolve(corM,k,dim=2)
kmeansmds3 <- cmdsolve(corM,k,dim=3)

kmeansmds4 <- cmdsolve(corM,k,dim=4)
kmeansmds10 <- cmdsolve(corM,k,dim=10)
kmeansmds20 <- cmdsolve(corM,k,dim=20)
kmeansmds30 <- cmdsolve(corM,k,dim=30)
kmeansnew <- kMeansOnDistances(corM,k)
kmeanscor <- kmeansCor(corM,k)
averagecornometric <- averageCorNoMetric(corM,k) 
averagecorcornometric <- averageCorCorNoMetric(corM,k) 
completecornometric <- completeCorNoMetric(corM,k)
completecorcornometric <- completeCorCorNoMetric(corM,k)
completecorcornometric <- completeCorCorNoMetric(corM,k)
totalclust <- kmeans(t(facs),centers=5,nstart=100)$cluster
fclust <- fclustering(corM,k)
averagecorcorcor <- averageCorCorCor(corM)
completecorcorcor <- completeCorCorCor(corM)
kmeansmdscor <- cmdsolveCor(corM,k)
kmeanscorcor <- kmeansCorCor(corM,k)
#kmeansneucor <- kMeansOnDistancesCor(corM,k)




#hier findet die tatsächliche Simulation statt, in die Ähnlichkeit der Clusterung für verschiedene 
#Stichprobengrößen vergelichen mit der Clusterung mit allen Daten (mit der gleichen Methode) für
#jeweils alle Methoden vergleichen wird

#toSimulate <- c("averagecor","completecor", "averagecorcor", "completecorcor", "kmeanscmd","kmeanscor")
#toSimulate <- c("averagecor")


#compares if there's a difference 





compareWith <- list("averagecor" =averagecor, "completecor" = completecor, "averagecorcor" = averagecorcor,
                    "completecorcor"= completecorcor,"kmeansmds" = kmeansmds, "kmeanskoord"=kmeansnew,
                    "kmeanscor"=kmeanscor, 
                    "averagecornom" = averagecornometric, "averagecorcornom" = averagecorcornometric, 
                    "completecornom"=completecornometric,"completecorcornom"=completecorcornometric
                    , "Dim4" = kmeansmds4, "Dim10" = kmeansmds10, "Dim20" = kmeansmds20, "Dim30" = kmeansmds30
                    , "Dim1" = kmeansmds1,"Dim2" = kmeansmds2,
                    "Dim3" = kmeansmds3, "totalclust"=totalclust, "averagecorcorcor"=averagecorcorcor,
                    "completecorcorcor"=completecorcorcor,"kmeansmdscor"=kmeansmdscor,"kmeanscorcor"=kmeanscorcor,"faclust"=fclust)



hierarchical <- c("averagecor", "averagecornom", "completecor","completecornom", "averagecorcor", "averagecorcornom",
                  "completecorcor", 
                  "completecorcornom")


#nonhierarchical <- c("totalclust","kmeanscmd","kmeanscor","kmeansneu")


allsmall <- c("averagecor","completecor", "averagecorcor", "completecorcor","kmeansmds","kmeanskoord","kmeanscor","faclust")


cmds <- c("Dim1","Dim2","Dim3","Dim4","Dim10","Dim20","Dim30","kmeansmds")


additionalCor <- c("averagecorcor", "averagecorcorcor","completecorcor",
                   "completecorcorcor", "kmeansmds", "kmeansmdscor", "kmeanscor","kmeanscorcor")

toSimulates <- c()

toSimulates[[1]] <- allsmall

#toSimulates[[1]] <- hierarchical 

#toSimulates[[2]] <- cmds

distribution <- F

names <- c("allsmall")#,  "allsmall")
#toSimulate <- hierarchical
#simulateClusterSamplesComparison(facs,toSimulate,compareWith, allnobs, 20,numbercluster=5,1 )



folder <- "/home/andreas/Desktop/Bachelorarbeit/strukturerhaltung/"
for(i in 1:length(toSimulates)) {
  jpeg(paste(folder,names[i],"/strukturerhaltung.jpeg",sep=""), width = 1200, height = 800)
  
  par(mfrow=c(2,2), oma=c(0,0,0,0), xpd=TRUE)
simulateClusterSamplesComparison(facs,toSimulates[[i]],compareWith, allnobs, nrep=100,numbercluster=5,compareMethod=1 )
  dev.off()
}





toSimulates[[1]] <- additionalCor

#toSimulates[[2]] <- cmds

distribution <- F

names <- c("additionalCor")#,  "allsmall")
#simulateClusterSamplesComparison(facs,toSimulate,compareWith, allnobs, 20,numbercluster=5,1 )




folder <- "/home/andreas/Desktop/Bachelorarbeit/strukturerhaltung/"
for(i in 1:length(toSimulates)) {
  jpeg(paste(folder,names[i],"/strukturerhaltung.jpeg",sep=""), width = 1200, height = 800)

  par(mfrow=c(2,2), oma=c(0,0,0,0), xpd=TRUE)
  simulateClusterSamplesComparison(facs,toSimulates[[i]],compareWith, allnobs, nrep=50,numbercluster=5,compareMethod=1 )
  dev.off()
}





