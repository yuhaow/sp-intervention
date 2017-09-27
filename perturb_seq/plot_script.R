rm(list=ls())

load('greedy_sp_summary.rda')
third_batch = do.call(cbind, list(false_positives[30,], true_positives[30,]))
colnames(third_batch) = c('fp', 'tp')

load('gies_summary.rda')
third_gies = do.call(cbind, list(false_positives[30,], true_positives[30,]))
colnames(third_gies) = c('fp', 'tp')

load('kernel_summary.rda')
kernel_batch = do.call(cbind, list(kernel_false_positives[30,], kernel_true_positives[30,]))
colnames(kernel_batch) = c('fp', 'tp')
tpplot(third_batch, third_gies, kernel_batch, 'plot', row.names(third_batch), colnames(kernel_true_positives))

tpplot <- function(gauss.plot, gies.plot, kernel.plot, title, alphas, alphas_k){
  gauss.plot = data.frame(gauss.plot)
  gies.plot = data.frame(gies.plot)
  kernel.plot = data.frame(kernel.plot)
  maxy = max(c(gauss.plot$tp, gies.plot$tp, kernel.plot$tp))
  maxx = max(c(gauss.plot$fp, gies.plot$fp, kernel.plot$fp))
  png(paste(title, ".png", sep=""))
  par(mar=c(5,5.5,3,1))
  plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), xlab="number of false positives", ylab="number of true positives", col="blue", type="p", pch="O", axes=FALSE, cex.lab=2)
  plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), axes=FALSE,  xlab="number of false positives", ylab="number of true positives", 
       col="blue", type="p", pch="O", cex.lab=2, lwd=3)
  axis(1, las=1, cex.axis=1.4)
  axis(2, las=1, cex.axis=1.4)
  box()
  text(gauss.plot$fp, gauss.plot$tp, labels=sapply(alphas, toString), cex=0.7, pos=3)
  text(kernel.plot$fp, kernel.plot$tp, labels=sapply(alphas_k, toString), cex=0.7, pos=3)
  lines(x=gies.plot$fp, y=gies.plot$tp, col="purple", type="p", pch=3)
  lines(x=kernel.plot$fp, y=kernel.plot$tp, col="red", type="p", pch=20)
  abline(0,0.218)
  legend("bottomright", legend=c("IGSP", "k-IGSP", "GIES"), col=c("blue", "red", "purple"), pch=c(1, 20, 3), cex=2)
  dev.off()
}
