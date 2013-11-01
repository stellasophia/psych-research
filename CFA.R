library(lavaan)

getClustering <- function(cor.sp, k, method) {
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
  cut.sp1
}


getCFASimiliarity <- function(facs, nrep, nobs,method) {
  
  for(i in 1:nrep) {
    daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
    
    daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp2 <- cor(as.matrix(daten.sp2), use="pairwise.complete.obs", method="pearson")
    
    k <- 5
    
    clustering <- getClustering(cor.sp1, k, method)
    
    latent.zuweisung <- ''
    
    for(i in 1:k) {
      clustername = as.character(paste0("class",i))
      cluster.name <- names(which(clustering==i))
      latent.zuweisung <- paste(latent.zuweisung, clustername, " =~ " )
      
      for(name in cluster.name) {
        latent.zuweisung <- paste(latent.zuweisung, name, "+ ")
      }
      
      #removes the , at the end
      
      latent.zuweisung <- substring(latent.zuweisung, 1, nchar(latent.zuweisung) - 2)
      
      latent.zuweisung <- paste(latent.zuweisung,  " ",  sep='\n')
      
    }
    
    frame <- as.data.frame(facs)
    fit <- cfa(latent.zuweisung,data = facs)
    
  }
}