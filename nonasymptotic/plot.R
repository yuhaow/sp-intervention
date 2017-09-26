library(pcalg)

#get results for SP algorithm
sp.consis <- function(prefix){
	load(paste(prefix, sep=""))
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag2essgraph(dag.list[[i]], t.list[[i]]))))
	sum(sapply(1:length(dag.list), function(i) {
		t1 <- wgtMatrix(grspdag.list[[i]])
		t2 <- wgtMatrix(dag2essgraph(dag.list[[i]], t.list[[i]]))
		sum(t1 != t2) == 0})) / length(dag.list)
	#mean(sapply(1:length(dag.list), function(i) shd(grspdag.list[[i]], dag.list[[i]])))
}

ps <- c(20)
neighs <- c(1.5)
ns <- c(100000, 10000, 1000)
dagnum <- 100
#ks <- c(3,5)
ks <- c(1, 2, 3)
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
cols <- rainbow(length(ns) * 2)

for(p in ps) for (k in ks) for(neigh in neighs){

	#proportion of consistently estimated DAGs
	fig.prefix <- paste("figure/p_", toString(p), "_neigh_", toString(neigh), "_k_", toString(k), sep="")
	png(paste(fig.prefix, "_consist.png", sep=""))
	par(mar=c(5,5.5,3,1))
	#get the data for each setting
	#greedy SP with restart
	a <- do.call(rbind, lapply(ns, function(n){ 
		prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), ".rda", sep="")
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
		prefix <- paste("p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), ".rda", sep="")
		rep(sp.consis(paste("ges/result/", prefix, ".rda", sep="")), length(alphas))
		}))
	for(i in 1:nrow(a))
		lines(x=log(alphas), y=a[i,], col=cols[i+length(ns)], type="o", lwd=3, pch=NA)
	if(p == 10 && k == 1){
		legend("bottom", legend=c(paste("IGSP n = ", sapply(ns, toString), sep=""), paste("GIES n = ", sapply(ns, toString), sep="")), col=cols, lty=rep(1, length(cols)), cex=1.8)
	}
	dev.off()
}
