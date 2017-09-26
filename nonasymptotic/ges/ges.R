library(pcalg)
library(graph)

set.seed(1)

args <- commandArgs()
print(args[4])
load(args[4])

#do joint estimation given single data
ges.alg <- function(data, targets){
	uniqtargets <- sapply(1:length(targets), function(i) !(i > 1 && sum(sapply(1:(i-1), function(j) setequal(targets[[i]], targets[[j]])))))
	tidx <- unlist(lapply(1:length(targets[uniqtargets]), function(t) rep(t, nrow(data[[t]]))))
	data <- do.call(rbind, data[uniqtargets])
	targets <- targets[uniqtargets]
	l0score <- new("GaussL0penIntScore", data = data, targets=targets, target.index=tidx, intercept = FALSE, use.cpp = TRUE)
	ges.fit <- gies(l0score)
	return(as(ges.fit$essgraph, "graphNEL"))
}

grspdag.list <- lapply(1:dagnum, function(i) ges.alg(data.list[[i]], t.list[[i]]))
save(grspdag.list, dag.list, t.list, file=paste("result/", basename(args[4]), ".rda", sep=""))
