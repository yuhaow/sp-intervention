library(pcalg)
library(graph)

set.seed(1)

args <- commandArgs()
print(args[4])
load(args[4])
k1 = as.numeric(args[5])
k2 = as.numeric(args[6])
k3 = as.numeric(args[7])

#do joint estimation given single data
ges.alg <- function(data, targets, k1, k2, k3){
	#set up the scoring function
	tidx <- unlist(lapply(1:length(targets), function(t) rep(t, nrow(data[[t]]))))
	data <- do.call(rbind, data)
	l0score <- new("GaussL0penIntScore", data = data, targets=targets, target.index=tidx, intercept = FALSE, use.cpp = TRUE)
	p <- ncol(data)

	#estimate essential graph
	essgraph <- new("EssGraph", nodes=sapply(1:p, toString), targets=targets, score=l0score)
	#start the forward phase
	essgraph$caus.inf(algorithm="GIES-F", maxSteps=k1)
	#start the backward phase
	essgraph$caus.inf(algorithm="GIES-B", maxSteps=k2)
	#start the turning phase
	essgraph$caus.inf(algorithm="GIES-T", maxSteps=k3)
	return(essgraph)
}

essgraph <- ges.alg(data.list, t.list, k1, k2, k3)

in_edges = essgraph$.in.edges
# obtain DAG matrix, where (i,j) = 1 indicates causal edge from i to j
grspdag = matrix(0, 24, 24)
idx = 1;
for (i in in_edges){
  for (j in i){
    grspdag[j, idx] = grspdag[j, idx] + 1
  }
  idx = idx+1
}
essgraph = as(essgraph, 'graphNEL')
save(essgraph, grspdag, t.list, file=paste("result/", basename(args[4]), "_", k1, "_", k2, "_", k3, ".rda", sep=""))
