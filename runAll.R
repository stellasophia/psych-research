source("Daten_einlesen.R")
source("Vergleichsverfahren.R")
source("GesamtDatenAnalysen.R")
source("cmdsolve.R")
source("faktorensimulation.R")

source("compareClusterings.R")
source("ClusterVergleichSimulation.R")
source("koorddendrograms.R")





allsmall <- c("averagecor","completecor", "averagecorcor", "completecorcor","kmeansmds","kmeanskoord","kmeanscor")


toSimulates <- c()


toSimulates[[1]] <- kmeanscmds

#toSimulates[[1]] <- hierarchical

#toSimulates[[1]] <- kmeanscmds
distribution <- T

names <- c("kmeanscmds")

#names <- c("hierarchical")

#names <- c("kmeanscmd")


distribution <- F
getComparison(toSimulates,names,distribution, addError=F)
distribution <- F
getComparison(toSimulates,names,distribution, addError=T)




allsmall <- c("averagecor","completecor", "averagecorcor", "completecorcor","kmeansmds","kmeanskoord","kmeanscor")



toSimulates <- c()

toSimulates[[1]] <- kmeanscmds

#toSimulates[[1]] <- hierarchical 

#toSimulates[[2]] <- cmds

distribution <- F

names <- c("kmeanscmds")#,  "allsmall")
#toSimulate <- hierarchical
simulateClusterSamplesComparison(facs,toSimulate,compareWith, allnobs, 20,numbercluster=5,1 )






folder <- "/home/andreas/Desktop/Bachelorarbeit2/strukturerhaltung/"
for(i in 1:length(toSimulates)) {
  jpeg(paste(folder,names[i],"/strukturerhaltung.jpeg",sep=""), width = 1200, height = 800)
  
  par(mfrow=c(2,2), oma=c(0,0,0,0), xpd=TRUE)
  simulateClusterSamplesComparison(facs,toSimulates[[i]],compareWith, allnobs, nrep=30,numbercluster=5,compareMethod=1 )
  dev.off()
}


