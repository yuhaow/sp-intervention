library(pcalg)
library(graph)

set.seed(1)

args <- commandArgs()
print(args[4])
load(args[4])

#do joint estimation given single data
ges.alg <- function(data, targets){
	#train the DAG
	tidx <- unlist(lapply(1:length(targets), function(t) rep(t, nrow(data[[t]]))))
	data <- do.call(rbind, data)
	l0score <- new("GaussL0penIntScore", data = data, targets=targets, target.index=tidx, intercept = FALSE, use.cpp = TRUE)
	gies.fit <- gies(l0score)
	return(gies.fit)
}

gies.fit <- ges.alg(data.list, t.list)

#get the output
p <- ncol(data.list[[1]])
wmat <- matrix(0, nrow=p, ncol=p)
for(i in 1:p){
	wmat[i, gies.fit$repr$.in.edges[[i]]] <- abs(gies.fit$repr$.params[[i]][3:length(gies.fit$repr$.params[[i]])])
	wmat[gies.fit$repr$.in.edges[[i]], i] <- wmat[i, gies.fit$repr$.in.edges[[i]]]
}

essgraph <- as(gies.fit$essgraph, "graphNEL")

save(essgraph, wmat, t.list, file=paste("result/", basename(args[4]), ".rda", sep=""))
