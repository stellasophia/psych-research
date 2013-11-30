#ruft die Vergleichsmethode auf, die zu verwenden ist
vergleich <- function(x,y,compareMethod) {
  result <- 0
  if(compareMethod == 1) {
   result  <- randIndex(x,y)
  # result  <- meilaHeckermann(x,y)
  }
  
  if(compareMethod==2) {
    result  <- meilaHeckermann(x,y)
  }
  
  
  if(compareMethod==3) {
    result  <- randIndexSensitivitySpecifity(x,y, result00)
  }
  result 
}




## erstellt für 2 clusterungen je eine Paarmatrix 
## und berechnet die Ähnlichkeiten dieser Matrizen (prozentuale Übereinstimmung der gleichen Felder)
randIndex <- function(x,y){
  g1 <- createPairMatrix(x)
  
  g2 <- createPairMatrix(y)
  
  counter1 <- 0
  counter2 <- 0
  for(i in 1:nrow(g1)) {
    for(j in 1:ncol(g2)) {
      if(i!=j) {
      if(g1[i,j] == g2[i,j]) {
        if(g1[i,j]  == 0 ) {
          counter1 <- counter1 + 1
        } else {
          counter2 <- counter2 +1
        }
      }
      }
    }
  }
  
  leng <- nrow(g1)
  
  counter1 <- counter1/(leng*(leng-1))
  
  counter2 <- counter2/(leng*(leng-1))
  
  result <- c(counter1,counter2)
}


#das Meila-Heckermann-Ähnlichkeitsmaß
#die Idee ist es zu jedem Cluster aus Clusterung C1 das ähnlichste Cluster aus C2 findet und die Anzahl der gleichen Elemente zählt
#der preozentuale Anteil dieser Elemente entspricht dann dem Meila-Heckermann-Ähnlichkeitsmaß
meilaHeckermann <- function(x,y) {
  
  t <- as.matrix(table(x,y)) ##Konfusionsmatrix wird gebildet
  minCluster <- min(max(dim(t[0,])),max(dim(t[,0]))) ##Die Anzahl der Cluster der Clusterung mit weniger Clustern wird ermittelt
  ClMatch <- 0 ##Übereinstimmungen in der Clusterung
  i <- 1
  while(i < minCluster){       
    w <- which(t==max(t), arr.ind=TRUE)##Zeile und Spalte des größten Eintrags in t ermitteln
    
    ClMatch <- ClMatch + max(t)
    
    j <- w[1,1]
    k <- w[1,2]
    t <- t[-j,-k] ##entsprechende Zeile / Spalte löschen
    
    i <- i+1
  }
  
  ClMatch <- ClMatch + max(t) ##Ein weiteres Mal den höchsten Wert addieren
  MM <- c(ClMatch/(2*length(x)-2), ClMatch/(2*length(x)-2)) ##Anzahl "richtige" Einträge/Gesamteinträge

  MM
}



##CLusterveähnlichkeitsvergleichsmethode3: 
##erstellt die Paarmatrix genau wie beim Rand-Index-Verfahren
##zählt den Anteil der 0en die in der zweiten Paarmatrix auftreten
#wenn sie auch schon in der ersten Paarmatrix 0en waren
##analog für 1en wenn resullt00 gleich false ist
randIndexSensitivitySpecifity <- function(x,y, result00){
  g1 <- createPairMatrix(x)
  
  g2 <-createPairMatrix(y)
  
  counter0 <- 0
  counter1 <- 0
  counter00 <- 0
  counter11 <- 0
  
  
  for(i in 1:nrow(g1)) {
    for(j in 1:ncol(g2)) {
      if(g1[i,j]  ==0 ) {
        counter0 <- counter0 + 1
        if(g1[i,j] == g2[i,j]) {
          counter00 <- counter00 + 1 
        }
      } else {
        counter1 <- counter1 +1
        if(g1[i,j] == g2[i,j]) {
          counter11 <- counter11 + 1 
        }
      }
    }
  }
  
  
  result <- c(counter00/counter0, counter11/counter1)
  
  result
}



## clusternumbercomparison  erstellt eine Matrix für jede Clusterung mit 1 an den Stellen wo 2 Items im gleichen Cluster liegen,sonst 0
createPairMatrix <- function(x){
  x <- as.numeric(x)
  g <- matrix(nrow=length(x),ncol=length(x))
  i <- 1
  j <- 1
  while(i <= length(x)){
    j <- 1
    while(j <= length(x)){
      
      if(x[i]==x[j]){
        g[i,j] <- 1
      }
      else{
        g[i,j] <- 0
      }
      
      j <- j+1
    }
    i <- i+1
  }
  g
}
