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
numclukmeans <- function (data,nobs,nrep,type="kmeans") {
  
  v <- c()
  for (i in 1:nrep){
    samps <- sample((x=1:nrow(data)), size=nobs, replace=T)
    print(samps)
    daten.sp <- data[samps,]
    
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    print(cor[1,2])
    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    result <- getClusterNumbers(points=points, type=type)
    number.cluster <- as.numeric(as.character(optimalScores(result)[,3]))
    print( number.cluster)
    v[[i]] <- number.cluster 
  }
  v
}

getClusterNumbers <- function(points, type="kmeans") {
  if(type=="kmeans") {
    result <- clValid(obj=points, nClust=2:15,clMethods="kmeans", validation=c("internal","stability"))
  } else if(type=="complete") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
  } else if(type=="average") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="average")
  }
  result
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


drawNumberComparisonKmeans <- function(facs,nrep, type="kmeans") {
  y<-11
  
  nobs <- c(100)
  vector <- c()
  
  for(nob in nobs) {
    
    vector <- as.vector(numclukmeans(data=facs, nobs=nob, nrep=nrep, type=type))
    print(vector)
    
    
    m <- matrix( ncol=length(vector[[1]]), nrow=length(vector))
    
    m <- matrix( ncol=length(vector), nrow=length(vector[[1]]))
    
    for(i in 1:length(vector)) {
      for(j in 1:length(vector[[1]])) {
        m[j,i] <- vector[[i]][j]
      }
    }
    
    par(mfrow=c(3,3))
    ###aufteilen auf vektoren der einzelnen Methoden, die dann geplottet werden
    for(i in 1:dim(m)[1]) {
    if(nob==nobs[1]) {
      drawBarplot(m[i,],ylab=paste("Faktorenanalyse"),nob=nob,type=type, cex.lab=1.5, method= measNames(result)[i])
    } else {
      drawBarplot(m[i,], ylab="relative Haeufigkeit",type=type, nob=nob, method=measNames(result)[i])
    }
   # drawExpectedLine(y)
  }
}

}



drawExpectedLine <-function(y) {
  abline(v=y-0.5, col="red", lwd=2)
}


drawBarplot <- function(data,ylab,nob,type,method,...) {

  data.rel <- table(data)/length(data)
  
  barplot(data.rel, main=paste("n=",nob,"clustertype= ", type, "method= ", method),  xlab="Anzahl Cluster", ylab=ylab, ylim=c(0,1))
}


completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs,nrep=30, type="average")

t <- drawKmeans(1,type="average") 

test <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
