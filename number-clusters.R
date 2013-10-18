library(NbClust)

#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluc <- function (data,nobs,nrep,method,h) {
  
  v <- 1:nrep
  for (i in 1:nrep){
    daten.sp <- data[sample((x=1:nrow(data)), size=nobs, replace=T),]
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    d <- as.dist((1 - (cor)))
    hca <- hclust(d, method, members=NULL)
    v[i] <- max(cutree(hca,h=h))
  }
  v
}


#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numclukmeans <- function (data,nobs,nrep,method,h) {
  
  v <- c()
  for (i in 1:nrep){
    daten.sp <- data[sample((x=1:nrow(data)), size=nobs, replace=T),]
    
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    result <- clValid(obj=points, nClust=2:15,clMethods="kmeans", validation=c("internal","stability"))
    number.cluster <- as.numeric(as.character(optimalScores(result)[,3]))
    print(measures(result))http://www.sueddeutsche.de/medien/ruhestand-im-rotlichtviertel-auf-arte-stellung-halten-1.1797499
    v[[i]] <- number.cluster 
  }
  v
}



#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
#Dieses Mal über die Korrelation der Korrelationen
numclucc <- function (data,nobs,nrep,method,h) {
  
  v <- 1:nrep
  for (i in 1:nrep){
    daten.sp <- data[sample((x=1:nrow(data)), size=nobs, replace=T),]
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    corcor <- cor(cor, use="pairwise.complete.obs", method="pearson")
    d <- as.dist((1 - (corcor)))
    hca <- hclust(d, method, members=NULL)
    v[i] <- max(cutree(hca,h=h))
  }
  v
}

#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
#Dieses Mal über die Korrelation der Korrelationen
# Dieses Mal Anzahl der Cluster ?ber H?he der Striche
numcluh <- function (data,nobs,nrep,method) {
  
  v <- 1:nrep
  for (i in 1:nrep){
    daten.sp <- data[sample((x=1:nrow(data)), size=nobs, replace=T),]
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    d <- as.dist((1 - (cor)))
    hca <- hclust(d, method, members=NULL)
    # plot(hca,hang=-1)
    h1 <- hca$height[-1]
    h2 <- c(hca$height[-length(hca$height)])
    
    v[i] <-length(hca$height)+3-which.max(h1-h2)
    
  }
  v
}

#berechnet für  mehrere Clusterungen bei bestimmten Cutoff-h
#und gibt erstellt ein Diagramm, in dem für verschiedene Stichprobengrößen
#die verschiedenen Clusteranzahlen die für dieses h sich ergeben miteinander vergleichen werden.
drawNumberComparison <- function(facs, nrep, method, h, cc) {
  bins=seq(0,20,by=1)
  
  if(method=="average") {
    if(!cc) {
      y<-13
    } else {
      y<-11
    }
  } else {
    if(!cc) {
      y<-4
    } else {
      y<-8
    }
  }
  
  
  nobs <- c(100,150,250,500,1000)
  varvector <- c()
  vector <- c()
  #erstellt den Plot in dem die verschiedenen Cluster-anzahlen berechnet werden!
  if(!cc)  {
    for(nob in nobs) {
      vector <- as.vector(numcluc(facs, nob, nrep, method, h))
      var <- mean(vector)
      varvector <- append(varvector, var)
      print(vector)
      if(nob==nobs[1]) {
        drawBarplot(vector,ylab=paste(method, " Korr"),nob=nob, cex.lab=1.5)
      } else {
        drawBarplot(vector, ylab="relative Haeufigkeit",nob=nob)
      }
      drawExpectedLine(y)
    }
  } else {
    for(nob in nobs) {
      vector <- as.vector(numclucc(facs, nob, nrep, method, h))
      var <- mean(vector)
      varvector <- append(varvector, var)
      if(nob==nobs[1]) {
        drawBarplot(vector,ylab=paste(method, "KorrKorr"),nob=nob,cex.lab=1.5)
      } else {
        drawBarplot(vector, ylab="relative Haeufigkeit",nob=nob)
      }
      drawExpectedLine(y)
    }
  }
  varvector
}

drawNumberComparisonFactor <- function(facs,nrep) {
  y<-11
  varvector <- c()
  nobs <- c(100,150,250,500,1000)
  
  for(nob in nobs) {
    vector <- as.vector(numclufactor(facs, nob, nrep))
    var <- mean(vector)
    varvector <- append(varvector, var)
    if(nob==nobs[1]) {
      drawBarplot(vector,ylab=paste("Faktorenanalyse"),nob=nob, cex.lab=1.5)
    } else {
      drawBarplot(vector, ylab="relative Haeufigkeit", nob=nob)
    }
    drawExpectedLine(y)
  }
  
  
}




drawNumberComparisonKmeans <- function(facs,nrep) {
  y<-11
  
  nobs <- c(100,150,250,500,1000)
  vector <- c()
  varvector <- c()
  for(nob in nobs) {
    vector <- as.vector(numclukmeans(facs, nob, nrep))
    var <- mean(vector)
    varvector <- append(varvector, var)
    if(nob==nobs[1]) {
      drawBarplot(vector,ylab=paste("Faktorenanalyse"),nob=nob, cex.lab=1.5)
    } else {
      drawBarplot(vector, ylab="relative Haeufigkeit", nob=nob)
    }
    drawExpectedLine(y)
  }
  
  
}



drawExpectedLine <-function(y) {
  abline(v=y-0.5, col="red", lwd=2)
}

drawBarplot <- function(data,ylab,nob,original,...) {
  print("drawPlot!")
  rel.table <- table(data)/length(data)
  barplot(rel.table, main=paste("ganz=",original),  xlab="Anzahl Cluster", ylab=ylab, ylim=c(0,1),...)
  
}



averageCorrelation <- drawNumberComparison(facs, 30, "average",h=0.8, F)

completeCorrelation <- drawNumberComparison(facs, 30, "complete",h=1, F)

averageCorrelationCorrelation <- drawNumberComparison(facs, 30, "average",h=0.8, T)

completeCCorrelationCorrelation <- drawNumberComparison(facs, 30, "complete",h=1, T)

drawKmeans <- function(nrep=5) {
  dim <- dim(corM)[1] - 1
  dist <- getDist(corM, F)
  fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
  points <- fit$points
  result <- clValid(obj=points, nClust=2:15,clMethods="kmeans", validation=c("internal","stability"))
  number.cluster.original <- as.numeric(as.character(optimalScores(result)[,3]))
  
  nobs <- c(100,150,250,500,1000)
  for(nob in nobs) {
    sample.values <- numclukmeans(facs, nob, nrep)
    
    par(mfrow=c(1,6))
    names <-  measNames(result)
    
    for(i in 1:length(names)) {
      name <- names[i]
      hist <- c()
      for(j in 1:nrep) {
        hist <- append(hist, sample.values[[j]][i])
      }
      drawBarplot(data=hist, ylab=name,nob=nob, original=number.cluster.original[i])
    }  
    
    
  }
}



completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs, 1)
