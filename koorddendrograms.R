addError=F
error <- ""
#erstellt die Dendrogramme als Grafiken
if(addError) {
  error <-"Fehler" 
}
folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/dendrogramComplete",error,".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
completeCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0.2,0,0,0,addError)
completeCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.25,0,0,0,addError)
completeCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.3,0,0,0,addError)
completeCor(corM,k=5,show=T)
dev.off()


folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/dendrogramAverage",error,".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
averageCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0.1,0,0,0,addError)
averageCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.15,0,0,0,addError)
averageCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.2,0,0,0,addError)
averageCor(corM,k=5,show=T)
dev.off()



folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/dendrogramAverageCorCor",error,".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
averageCorCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0.1,0,0,0,addError)
averageCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.15,0,0,0,addError)
averageCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.2,0,0,0,addError)
averageCorCor(corM,k=5,show=T)
dev.off()


folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/dendrogramCompleteCorCor",error,".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
completeCorCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0.1,0,0,0,addError)
completeCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.15,0,0,0,addError)
completeCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0.2,0,0,0,addError)
completeCorCor(corM,k=5,show=T)
dev.off()




addError <- F

folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/KorrelationzunehmenddendrogramComplete",".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
completeCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0,0,0.5,0,addError)
completeCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.6,0,addError)
completeCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.7,0,addError)
completeCor(corM,k=5,show=T)
dev.off()



folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/KorrelationzunehmenddendrogramAverage",".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
averageCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0,0,0.5,0,addError)
averageCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.6,0,addError)
averageCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.7,0,addError)
averageCor(corM,k=5,show=T)
dev.off()






folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/KorrelationzunehmenddendrogramCompleteCor",".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
completeCorCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0,0,0.5,0,addError)
completeCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.6,0,addError)
completeCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.7,0,addError)
completeCorCor(corM,k=5,show=T)
dev.off()



folder <- "/home/andreas/Desktop/Bachelorarbeit/"
jpeg(paste(folder,"/KorrelationzunehmenddendrogramAverageCor",".jpeg",sep=""), width = 1200, height = 800)
par(mfrow=c(2,2), oma=c(0,0,0,0), cex=1.1)
corM <-setCorrelationMatrix (0,0,0,0,addError)
averageCorCor(corM,k=5,show=T)
corM <-setCorrelationMatrix (0,0,0.5,0,addError)
averageCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.6,0,addError)
averageCorCor(corM,k=5,show=T)
corM <- setCorrelationMatrix(0,0,0.7,0,addError)
averageCorCor(corM,k=5,show=T)
dev.off()




