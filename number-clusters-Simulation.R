library(gridExtra)

source("faktorensimulation.R")
source("cmdsolve.R")
source("Vergleichsverfahren.R")


#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluadvanced.simulation <- function (cor,nrep,type="kmeans") {
  
  v <- c()
  for (i in 1:nrep){

    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    result <<- getClusterNumbers(points=points, type=type)
    number.cluster <- as.numeric(as.character(optimalScores(result)[,3]))
    print( number.cluster)
    v[[i]] <- number.cluster 
  }
  v
}





getClusterNumbers.simulation <- function(points, type="kmeans") {
  if(type=="kmeans") {
    result <- clValid(obj=points, nClust=2:15,clMethods="kmeans", validation=c("internal","stability"))
  } else if(type=="complete") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
  } else if(type=="average") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="average")
  }
  result
}


drawNumberClusterAdvanced.simulation <- function(cor,nrep, type="kmeans") {
  
  whole.cluster.number <- c()

    whole.cluster.number <- c(5,5,5,5,5,5,5)
  
  
  result.names <- c("whole", "Bias")  
  nobs <- c(500)
  
  
  resultsmatrix <- matrix(nrow=length(result.names), ncol=length(method.names))
  rownames(resultsmatrix) <- result.names
  colnames(resultsmatrix) <- method.names
  
  for(o in 1:length(nobs)) {
    vector <- as.vector(numcluadvanced.simulation(cor, nrep=nrep, type=type))
    
    m <- matrix( ncol=length(vector), nrow=length(vector[[1]]))
    
    for(i in 1:length(vector)) {
      for(j in 1:length(vector[[1]])) {
        m[j,i] <- vector[[i]][j] - whole.cluster.number[j]
      }
    }
    
    
    
    # par(mfrow=c(3,3))
    ###aufteilen auf vektoren der einzelnen Methoden, die dann geplottet werden
    for(i in 1:dim(m)[1]) {
      #    drawBarplot(m[i,],ylab=paste("Faktorenanalyse"),nob=nob,type=type, cex.lab=1.5, method= measNames(result)[i])
      method= measNames(result)[i]
      method.var <- var(m[i,])
      method.bias <- mean(m[i,])
      method.whole <-  whole.cluster.number[i]
      
      resultsmatrix[1,i] <- method.whole
      resultsmatrix[2,i] <- method.bias
  
    }
    
  }
  
  resultsmatrix
  
}

method.names <<- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")

getClusterNumberBias.simulation <- function(NL.mus,  Kor.mus,  type="kmeans") {
  
  r.names <- c()
  descriptions <- ""
  for(i in 1:length(NL.mus)) {
    r.names[i] <- paste0("constr ", i)
    descriptions <- paste0(descriptions,  " und " , r.names[i] , " mit NL von ",
                           NL.mus[i], " und Faktorkorrelation von ", Kor.mus[i], " \n  ")
  }
  
 
  
  rs <- matrix(nrow=length(NL.mus), ncol=length(method.names))
  rownames(rs) <- r.names
  colnames(rs) <- method.names
  for(i in 1:length(NL.mus)) {
    corM1 <- setCorrelationMatrix(NL.mus[i],0,Kor.mus[i],0, F)
    r1 <- drawNumberClusterAdvanced.simulation(corM1,1,type="kmeans")
    rs[i,] <- r1[2,]
  }
  paintTable(rs, "Clusteranzahlsabweichung bei EFA-Similation", 
             paste0( "type=", type, " \n ", descriptions)) 
}






getClusterNumberBias.simulation.methods <- function(method, fa.ges) {
  
  r.names <- c()
  descriptions <- ""
 
  
  descriptions<- "bei allen NL gleich und entsprechen Kommunalität (NL.equal)"
  
  if(method==2) {
    descriptions<- "bei einer NL und entsprechen Kommunalität (NL.one)"
  } else if(method==3) {
    descriptions<- "bei zwei NL und entsprechen Kommunalität (NL.two)"
  }
  
  
  
  loads <- NL.equal(fa.ges$loadings)
  
  # loads <- NL.fixed(fa.ges$loadings, 0.2)
  
  if(method==2) {
    loads <- NL.one(fa.ges$loadings)
  } else if(method == 3) {
    loads <- NL.two(fa.ges$loadings)
  }
  
  zuordnung.ges <- apply(loads,1,function(x) which.max(abs(x)))
  
  Phi <- Phi.fixed(fa.ges$Phi, 0)
  corM <- sim.structure(fx=loads,Phi=Phi,n=0)$model
  
  types <- c("kmeans", "average", "complete")
  
  rs <- matrix(nrow=length(types), ncol=length(method.names))
  rownames(rs) <- types
  colnames(rs) <- method.names

  
  for(i in 1:length(types)) {
  r <- drawNumberClusterAdvanced.simulation(corM,1,types[i])
  rs[i,] <- r[2,]
  }

  paintTable(rs, "Clusteranzahlsabweichung bei EFA-Similation mit 5 Faktoren", 
             paste0(" \n ", descriptions)) 
}

#NL.mus <- c(0,0.1,0.2,0.1)
#Kor.mus <- c(0,0,0,0.4)

#getClusterNumberBias.simulation(NL.mus,Kor.mus,type="kmeans")


