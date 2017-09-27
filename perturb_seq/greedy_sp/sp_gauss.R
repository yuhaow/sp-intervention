library(pcalg)
library(graph)

#input to Greedy SP with random restarts: mat -- precision matrix, order -- initial permutation
sp.restart.alg <- function(suffstat, intsuffstat, inttargets, alpha){
	#set up the initial parameters for all functions
	p <- ncol(suffstat$C)

	#I-contradicting edges
	contCItest <- function(i, j, s){
		k1 <- which(sapply(inttargets, function(a) i %in% a))
		k2 <- which(sapply(inttargets, function(a) j %in% a))
		k11 <- setdiff(k1, k2)
		k22 <- setdiff(k2, k1)
		return((length(k11) > 0 || length(k22) > 0) && ifelse(length(k11) > 0, !sum(sapply(k11, function(t) gaussCItest(i, j, c(), intsuffstat[[t]]) < alpha)), 1) 
			&& ifelse(length(k22) > 0, !sum(sapply(k22, function(t) gaussCItest(i, j, c(), intsuffstat[[t]]) > alpha)), 1))
	}

	#I-covered edges
	icovtest <- function(i, j, s){
		k1 <- which(sapply(inttargets, function(a) i %in% a))
		k2 <- which(sapply(inttargets, function(a) j %in% a))
		k <- setdiff(k1, k2)
		length(k) > 0 && sum(sapply(k, function(t) gaussCItest(i, j, s, intsuffstat[[t]]) < alpha)) 
	}

	#get new dag based on edge flip
	get.newdag <- function(dag, contdag, order, edge, vorders){


		#get the new orders
		a <- which(order == edge[1])
		b <- which(order == edge[2])
		order <- order[c(0:(a-1), b, a:(b-1), (b+1):(p+1))[2:(p+1)]]
		#check if the new order has been visited
		if(list(order) %in% vorders) return(NULL)
		#if it has not been visited, check if this edge is an I-covered edge
		revorder <- sapply(1:p, function(t) which(order==t))
		par <- subset(1:p, dag[,edge[1]] == 1)
		if(icovtest(edge[1], edge[2], order[c(1:(revorder[edge[1]]-1))[-revorder[edge[2]]]])) return(NULL)
		#then you can continue
		dag[edge[1], edge[2]] <- 0
		dag[edge[2], edge[1]] <- 1
		contdag[edge[1], edge[2]] <- 0
		contdag[edge[2], edge[1]] <- contCItest(edge[2], edge[1], par)
		#parent set of the flipped components
		if(length(par) != 0){
			dag[par, edge[1]] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[1], c(par[-i], edge[2]), suffstat) < alpha)
			dag[par, edge[2]] <- sapply(1:length(par), function(i) gaussCItest(par[i], edge[2], par[-i], suffstat) < alpha)
			contdag[par, edge[1]] <- sapply(1:length(par), function(i) if(dag[par[i], edge[1]] != 0) contCItest(par[i], edge[1], c(par[-i], edge[2])) else 0)
			contdag[par, edge[2]] <- sapply(1:length(par), function(i) if(dag[par[i], edge[2]] != 0) contCItest(par[i], edge[2], par[-i]) else 0)
		}
		#get updates of the number of contradicting edges
		return(list(dag=dag, contdag=contdag, order=order))
	}

	#get the initial dag
	init.dag <- function(order){
		revorder <- sapply(1:p, function(t) which(order==t))
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(revorder[i] < revorder[j]) gaussCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]], suffstat) < alpha else 0)))
	}

	#get the initial dag
	init.contdag <- function(dag, order){
		revorder <- sapply(1:p, function(t) which(order==t))
		return(sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i, j] != 0) contCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]]) else 0)))
	}

	#the stack for visited orders
	sing.restart <- function(order){
		vorders <- list()
		vtrace <- list()
		vdags <- list()
		dag <- init.dag(order)
		contdag <- init.contdag(dag, order)
		mindag <- list(dag=dag, n=sum(contdag != 0))
		while(TRUE){
			#get the list of covered edges
			cov.edge <- which(dag != 0, arr.ind = TRUE)
			cov.edge <- data.frame(subset(cov.edge, apply(cov.edge, 1, function(x) all.equal(c(dag[-x[1], x[1]]), c(dag[-x[1], x[2]])) == TRUE)))
			#get the list of DAGs after I-covered edge reversals
			rdags <- if(nrow(cov.edge) > 0) apply(cov.edge, 1, function(edge) get.newdag(dag, contdag, order, edge, vorders)) else list()
			if(length(rdags) > 0) rdags <- subset(rdags, sapply(rdags, function(t) !is.null(t)))
			select <- which(sapply(rdags, function(rdag) sum(rdag$dag != 0) < sum(dag != 0)) == TRUE)
			#start the searching
			if((length(rdags) > 0 && length(vtrace) != 3) || length(select) != 0){
				if(length(select) != 0){
					vorders <- list()
					vtrace <- list()
					vdags <- list()
					order <- rdags[[select[1]]]$order
					dag <- rdags[[select[1]]]$dag
					mindag <- list(dag=dag, n=sum(rdags[[select[1]]]$contdag != 0))
				}else{
					vorders <- append(vorders, list(order))
					vtrace <- append(vtrace, list(order))
					vdags <- append(vdags, list(dag))
					order <- rdags[[1]]$order
					dag <- rdags[[1]]$dag
					if(sum(rdags[[1]]$contdag != 0) < mindag$n) mindag <- list(dag=dag, n=sum(rdags[[1]]$contdag != 0))
				}
			}else{
				if(length(vtrace) == 0)
					break
				vorders <- append(vorders, list(order))
				order <- tail(vtrace, 1)[[1]]
				vtrace <- head(vtrace, -1)
				dag <- tail(vdags, 1)[[1]]
				vdags <- head(vdags, -1)
			}
		}
		print(order)
		return(mindag)
	}

	#main part of the algorithm
	start.order <- lapply(1:10, function(x) sample(1:p, p, replace=F))
	dag.list <- lapply(start.order, function(x) sing.restart(x))
	edgenum.list <- sapply(dag.list, function(dag) sum(dag$dag != 0))
	minidx <- which(edgenum.list == min(edgenum.list))
	contedgenum.list <- sapply(dag.list, function(dag) dag$n)
	minidx <- minidx[which.min(contedgenum.list[minidx])]
	return(dag.list[[minidx]]$dag)
}

#get the test significance for estimated edges
conf.mat <- function(dag, order, suffstat){
	p <- ncol(suffstat$C)
	revorder <- sapply(1:p, function(t) which(order==t))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) 1 - gaussCItest(i, j, order[c(1:(revorder[j]-1))[-revorder[i]]], suffstat) else 0))
	t <- sapply(1:p, function(j) sapply(1:p, function(i) if(dag[i,j] != 0) t[i, j] else t[j, i]))
}

#get data for simulation
args <- commandArgs()
load(args[4])
print(title)

set.seed(1)

#get result for each dag
suffstat <- list(C=cor(data.list[[1]]), n=nrow(data.list[[1]]))
intsuffstat <- lapply(2:length(t.list), function(t) list(C=cor(data.list[[t]]), n=nrow(data.list[[t]])))
inttargets <- t.list[2:length(t.list)]

#get the estiated graph
#alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)
alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
grspdag.list <- lapply(alphas, function(alpha) sp.restart.alg(suffstat, intsuffstat, inttargets, alpha))
essgraph.list <- lapply(grspdag.list, function(grspdag) dag2essgraph(as(grspdag, "graphNEL"), targets=t.list))
grspdag.list <- lapply(grspdag.list, function(grspdag) as(grspdag, "graphNEL"))
save(alphas, essgraph.list, grspdag.list, t.list, file=paste("result/", basename(args[4]), ".rda", sep=""))
