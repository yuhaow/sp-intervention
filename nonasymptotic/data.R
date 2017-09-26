library(pcalg)
library(graph)
library(MASS)

#parameter used for generating data
ps <- c(20)
neighs <- c(1.5)
ns <- c(1000, 10000, 100000)
dagnum <- 100
ks <- c(1)

#set random seed
set.seed(1)

#function for generating random dags
randgraph <- function(p, neigh){
	dag <- randomDAG(p, prob = neigh / (p - 1), lB = 0.1, uB = 1)
	dag <- as(dag, "matrix") * matrix(sample(c(-1, 1), p * p, replace=TRUE), nrow=p)
	dag <- lapply(1:p, function(i) list(edges=which(dag[i,] != 0), weights=dag[i, dag[i,] != 0]))
	names(dag) <- 1:p
	dag <- new("graphNEL", node=sapply(1:p, toString), edgeL=dag, edgemode="directed")
	return(dag)
}

#generate random graph for each specific K
randngraph <- function(p, dag, nodes){
	dag <- as(dag, "matrix")
	dag[, nodes] <- 0
	dag <- lapply(1:p, function(i) list(edges=which(dag[i,] != 0), weights=dag[i, dag[i,] != 0]))
	names(dag) <- 1:p
	dag <- new("graphNEL", node=sapply(1:p, toString), edgeL=dag, edgemode="directed")
	return(dag)
}

#generate random graphs
for(p in ps){
for(neigh in neighs){

	#get list of observational dags
	ori.dag.list <- lapply(1:dagnum, function(i) randgraph(p, neigh))
	#get random permutation for each dag
	permut <- lapply(1:dagnum, function(i) sample(p, p))
	#get error variances for the data-generating dag

	#start generating data and interventions
	for(n in ns){
	for(k in ks){

		#get intervented nodes
		vpert <- 10
		nodes.list <- lapply(1:dagnum, function(i) sample(1:p, vpert * k))
		t.list <- lapply(nodes.list, function(nodes) lapply(1:k, function(t) nodes[(vpert * (t-1) + 1) : (vpert * t)]))

		#get all K dags
		dag.list <- lapply(1:dagnum, function(i) c(list(ori.dag.list[[i]]), lapply(1:k, function(t) randngraph(p, ori.dag.list[[i]], t.list[[i]][[t]]))))

		#get data for all K dags under interventional case
		data.list <- lapply(1:dagnum, function(i) c(list(mvrnorm(n, mu=rep(0,p), Sigma=trueCov(dag.list[[i]][[1]]))), 
			lapply(1:k, function(t) {
				err.vec <- rep(1, p)
				err.vec[t.list[[i]][[t]]] <- 0.2
				sigmat <- solve(diag(p) - wgtMatrix(dag.list[[i]][[t+1]]))
				sigmat <- sigmat %*% diag(err.vec) %*% t(sigmat)
				mvrnorm(n, mu=rep(0,p), Sigma=sigmat)})))

		#inference with permutations 
		data.list <- lapply(1:dagnum, function(i) lapply(1:(k+1), function(j) {
			t <- data.list[[i]][[j]][, permut[[i]]]
			colnames(t) <- 1:p
			return(t)}))
		dag.list <- lapply(1:dagnum, function(i) {
			dagmat <- as(ori.dag.list[[i]], "matrix")
			dagmat <- dagmat[permut[[i]], permut[[i]]]
			colnames(dagmat) <- 1:p
			rownames(dagmat) <- 1:p
			return(as(abs(dagmat), "graphNEL"))})
		t.list <- lapply(1:dagnum, function(i) c(list(integer(0)), lapply(t.list[[i]], function(nodes) as.integer(sapply(nodes, function(t) which(permut[[i]] == t))))))

		title <- paste("data/p_", toString(p), "_neigh_", toString(neigh), "_n_", toString(n), "_k_", toString(k), ".rda", sep="")
		intnum <- k
		save(data.list, dag.list, dagnum, intnum, t.list, title, file=title)
		print(title)
	}}
}}
