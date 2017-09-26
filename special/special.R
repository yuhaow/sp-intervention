library(pcalg)
library(graph)
library(MASS)

#parameter used for generating data
ps <- c(7)
ns <- c(1000, 10000, 100000)
lbs <- c(0.25, 0.1, 0.5)
dagnum <- 100

#set random seed
set.seed(1)

randomGeneReg <- function (n_TF, n_genes, n_levels, prob, lB, uB){
    stopifnot(is.numeric(lB), is.numeric(uB), lB <= uB)
    n=n_TF+n_genes
    V <- as.character(1:n)
    edL <- as.list(V)
    names(edL) <- V
    genWgt <- function(n, min, max) {
      w <- runif(n, min = min, max = max)
      sgn <- sample(c(-1,1), n, replace = TRUE)
      w*sgn
    }
    if (n_levels==1) {
	levels <- c(n_TF, n_genes)
    }
    else {
	levels <- tabulate(c(rep(1,n_TF), c(2:(n_levels+1)), sample(2:(n_levels+1), size = n_genes-n_levels, replace = TRUE)))
    }
    nodes <- list()
    nodes[[1]] <- 1:n_TF
    for (i in 2:(n_levels+1)) {
	nodes[[i]] <- (nodes[[i-1]][levels[i-1]]+1):(nodes[[i-1]][levels[i-1]]+levels[i])
    }
    parents <- matrix(rep(0,n*n), nrow=n)
    jj <- n_TF+1
    for (j in 2:(n_levels+1)){
    	for (k in 1:levels[j]){
    		parents[jj,nodes[[j-1]]] = rbinom(levels[j-1],1,prob)
			jj = jj+1
		}
	}
    for (i in 1:n){
	if (sum(parents[,i])==0) {
		edL[[i]] <- list(edges = integer(0), weights = numeric(0))	
	}
	else {
		edgeList <- which(parents[,i]==1)
		weightList <- genWgt(length(edgeList), min = lB, max = uB)
		edL[[i]] <- list(edges = edgeList, weights = weightList)
	}
    }
    res <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    res
}

#function for generating random dags
randgraph <- function(lb){
	dag <- matrix(0, nrow=p, ncol=p)
	dag[1, 4] <- 1
	dag[2, 4] <- 1
	dag[3, 4] <- 1
	dag[3, 7] <- 1
	dag[4, 6] <- 1
	dag[5, 7] <- 1
	dag <- dag * (matrix(runif(p * p, min = lb, max=1), nrow=p) * matrix(sample(c(-1, 1), p * p, replace=TRUE), nrow=p))
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
for(p in ps) for(lb in lbs){

	#get list of observational dags
	#ori.dag.list <- lapply(1:dagnum, function(i) randomGeneReg(5,15,5,0.5,0.1,1))
	ori.dag.list <- lapply(1:dagnum, function(i) randgraph(lb))
	#get random permutation for each dag
	permut <- lapply(1:dagnum, function(i) sample(p, p))
	#get error variances for the data-generating dag
	err.list <- lapply(1:dagnum, function(i) runif(p, min=1, max=1.5))

	#start generating data and interventions
	for(n in ns){

		#get intervented nodes
		t.list <- lapply(1:dagnum, function(t) list(list(c(4), c(7)), list(c(4)), list(c(7)))[[sample(1:3, 1)]])

		#get all K dags
		dag.list <- lapply(1:dagnum, function(i) c(list(ori.dag.list[[i]]), lapply(t.list[[i]], function(t) randngraph(p, ori.dag.list[[i]], t))))

		#get data for all K dags under interventional case
		data.list <- lapply(1:dagnum, function(i) c(list(mvrnorm(n, mu=rep(0,p), Sigma=trueCov(dag.list[[i]][[1]]))), 
			lapply(1:k, function(t) {
				err.vec <- rep(1, p)
				err.vec[t.list[[i]][[t]]] <- 0.2
				sigmat <- solve(diag(p) - wgtMatrix(dag.list[[i]][[t+1]]))
				sigmat <- sigmat %*% diag(err.vec) %*% t(sigmat)
				mvrnorm(n, mu=rep(0,p), Sigma=sigmat)})))

		#inference with permutations 
		data.list <- lapply(1:dagnum, function(i) lapply(1:(length(t.list[[i]])+1), function(j) {
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

		title <- paste("data/p_", toString(p), "_n_", toString(n), "_lb_", toString(lb), ".rda", sep="")
		save(data.list, dag.list, dagnum, t.list, title, file=title)
		print(title)
	}
}
