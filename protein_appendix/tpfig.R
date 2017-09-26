library(graph)
library(pcalg)

#ground truth dag
ground.truth <- list(c(1,2), c(2,6), c(3,4), c(3,9), c(4,9), c(5,3), c(5,4), c(5,7), c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), c(9,1), c(9,2), c(9,10), c(9,11))
#ground.truth <- list(c(1,2), c(2,6), c(3,4), c(3,9), c(4,9), c(5,3), c(5,4), c(5,7), c(8,1), c(8,2), c(8,6), c(8,7), c(8,10), c(8,11), c(9,1), c(9,2), c(9,10), c(9,11), c(9,8), c(6,7))
dag <- matrix(0, nrow=11, ncol=11)
for (t in ground.truth)
	dag[t[1], t[2]] <- 1

#get true positive rate
tprate <- function(essgraph){
	essmat <- as(essgraph, "matrix")
	edgenum <- sum(essmat | t(essmat)) / 2
	tp <- sum(essmat & dag)
	skel.tp <- sum((essmat | t(essmat)) & dag)
	return(list(tp=tp, fp=edgenum-tp, skel.tp=skel.tp, skel.fp=edgenum-skel.tp))
}

#plot new figures
tpplot <- function(gauss.plot, gies.plot, kernel.plot, title){
	maxy = max(c(gauss.plot$tp, gies.plot$tp, kernel.plot$tp))
	maxx = max(c(gauss.plot$fp, gies.plot$fp, kernel.plot$fp))
	png(paste("figure/safe_choice_", title, ".png", sep=""))
	par(mar=c(5,5.5,3,1))
	plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), xlab="number of false positives", ylab="number of true positives", col="blue", type="p", pch="O", axes=FALSE, cex.lab=2)
	plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), axes=FALSE,  xlab="number of false positives", ylab="number of true positives", 
		col="blue", type="p", pch="O", cex.lab=2, lwd=3)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	text(gauss.plot$fp, gauss.plot$tp, labels=sapply(alphas, toString), cex=0.7, pos=3)
	text(kernel.plot$fp, kernel.plot$tp, labels=sapply(alphas, toString), cex=0.7, pos=3)
	lines(x=gies.plot$fp, y=gies.plot$tp, col="purple", type="p", pch=3)
	lines(x=kernel.plot$fp, y=kernel.plot$tp, col="red", type="p", pch=20)
	legend("bottomright", legend=c("IGSP", "k-IGSP", "GIES"), col=c("blue", "red", "purple"), pch=c(1, 20, 3), cex=2)
	dev.off()
}

skel.tpplot <- function(gauss.plot, gies.plot, kernel.plot, title){
	maxy = max(c(gauss.plot$tp, gies.plot$tp, kernel.plot$tp))
	maxx = max(c(gauss.plot$fp, gies.plot$fp, kernel.plot$fp))
	png(paste("figure/safe_choice_", title, ".png", sep=""))
	par(mar=c(5,5.5,3,1))
	plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), xlab="number of false positives", ylab="number of true positives", col="blue", type="p", pch="O", axes=FALSE, cex.lab=2)
	plot(x=gauss.plot$fp, y=gauss.plot$tp, xlim=c(0, maxx), ylim=c(0, maxy), axes=FALSE,  xlab="number of false positives", ylab="number of true positives", 
		col="blue", type="p", pch="O", cex.lab=2, lwd=3)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	text(gauss.plot$fp, gauss.plot$tp, labels=sapply(alphas, toString), cex=0.7, pos=3)
	text(kernel.plot$fp, kernel.plot$tp, labels=sapply(alphas, toString), cex=0.7, pos=3)
	lines(x=gies.plot$fp, y=gies.plot$tp, col="purple", type="p", pch=3)
	lines(x=kernel.plot$fp, y=kernel.plot$tp, col="red", type="p", pch=20)
	coef <- 15 / (11 * 5)
	lines(x=0:33, y=coef * (0:33), type="l")
	legend("bottomright", legend=c("IGSP", "k-IGSP", "GIES"), col=c("blue", "red", "purple"), pch=c(1, 20, 3), cex=2)
	dev.off()
}

#get results for Gaussian CI tests
load("sp/result/safe_choice.rda.rda")
sp.list <- lapply(essgraph.list, tprate)

#get results for kernel CI tests
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
ksp.list <- lapply(alphas, function(alpha){
	load(paste("sp/result/safe_choice.rda_alpha_", toString(alpha), "_method_hsic.gamma.rda", sep=""))
	tprate(dag2essgraph(grspdag, targets=t.list))})

#plot the results, first plot the directed rates
gauss.plot <- list(tp=sapply(sp.list, function(t) t$tp), fp=sapply(sp.list, function(t) t$fp))
gies.plot <- list(tp=c(0,0, rep(1,6), rep(2,6), 3,3, rep(4,5), rep(5,19), 6,6,5,5,6,6), fp=c(0,1, 1:6, 6:11, 11,12, 12:16, 16:24, 26:35, 35,36,37,38,38,39))
kernel.plot <- list(tp=sapply(ksp.list, function(t) t$tp), fp=sapply(ksp.list, function(t) t$fp))
tpplot(gauss.plot, gies.plot, kernel.plot, "directed")

#plot the results, first plot the skeleton
gauss.plot <- list(tp=sapply(sp.list, function(t) t$skel.tp), fp=sapply(sp.list, function(t) t$skel.fp))
gies.plot <- list(tp=c(0,1,2,2,3,3,4,4, rep(5,5), 6,7,7, rep(8,5), 9,9,9,10,10, rep(11, 15), 12,12,12,13,14,13), 
	fp=c(0,0,0,1,1,2,2,3,3, 4:6, 7,7,7, 8,8,9,10,11,12,12,13,14,14,15, 15:29, 29,30,31,31,31,32))
kernel.plot <- list(tp=sapply(ksp.list, function(t) t$skel.tp), fp=sapply(ksp.list, function(t) t$skel.fp))
skel.tpplot(gauss.plot, gies.plot, kernel.plot, "skeleton")
