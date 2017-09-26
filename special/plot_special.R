library(pcalg)

#get results for SP algorithm
sp.consis <- function(prefix){
	load(paste(prefix, sep=""))
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag2essgraph(dag.list[[i]], t.list[[i]]))))
	sum(sapply(1:length(dag.list), function(i) {
		t1 <- as(grspdag.list[[i]], "matrix")
		t2 <- as(dag2essgraph(dag.list[[i]], t.list[[i]]), "matrix")
		sum(t1 != t2) == 0})) / length(dag.list)
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag.list[[i]])))
}

ps <- c(7)
ns <- c(100000, 10000, 1000)
lbs <- c(0.25, 0.1, 0.5)
dagnum <- 100
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
cols <- rainbow(length(ns) * 2)

for(p in ps) for (lb in lbs) {

	#proportion of consistently estimated DAGs
	fig.prefix <- paste("figure/p_", toString(p), "_lb_", toString(lb), sep="")
	png(paste(fig.prefix, "_consist.png", sep=""))
	par(mar=c(5,5.5,3,1))
	#get the data for each setting
	#greedy SP with restart
	a <- do.call(rbind, lapply(ns, function(n){ 
		prefix <- paste("p_", toString(p), "_n_", toString(n), "_lb_", toString(lb), ".rda", sep="")
		sapply(alphas, function(alpha) sp.consis(paste("sp/result/", prefix, "_alpha_", toString(alpha), "_restart.rda", sep="")))
		}))
	plot(x=log(alphas), y=a[1,], ylim=c(0, 1), axes=FALSE,  xlab="log(alpha)", ylab="proportion of consistently estimated DAGs", 
		col=cols[1], type="o", cex.lab=1.9, lwd=3, pch=NA)
	axis(1, las=1, cex.axis=1.4)
	axis(2, las=1, cex.axis=1.4)
	box()
	for(i in 2:nrow(a))
		lines(x=log(alphas), y=a[i,], col=cols[i], type="o", lwd=3, pch=NA)
	a <- do.call(rbind, lapply(ns, function(n){
		prefix <- paste("p_", toString(p), "_n_", toString(n), "_lb_", toString(lb), ".rda", sep="")
		rep(sp.consis(paste("ges/result/", prefix, ".rda", sep="")), length(alphas))
		}))
	for(i in 1:nrow(a))
		lines(x=log(alphas), y=a[i,], col=cols[i+length(ns)], type="o", lwd=3, pch=NA)
	legend("bottomright", legend=c(paste("IGSP n = ", sapply(ns, toString), sep=""), paste("GIES n = ", sapply(ns, toString), sep="")), col=cols, lty=rep(1, length(cols)), cex=1.8)
	dev.off()
}
