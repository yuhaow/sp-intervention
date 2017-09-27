rm(list=ls())

true_descendents = read.table('log_pvals', header = F)
perturbations = read.csv('perturbations', stringsAsFactors = F, header = F)
perturbations = perturbations$V1
colnames(true_descendents) = perturbations
true_descendents = abs(true_descendents)

# delete control column
true_descendents = true_descendents[c(1,3:9)]
perturbations = perturbations[c(1,3:9)]


gene_names = read.csv('genes', stringsAsFactors = F, header = F)
gene_names = gene_names$V1

files = list.files('./result')

#Greedy sp
DAG_list = list()
DAG_names = c()
for (file in files) {
  load(paste('./result/', file, sep = ''))
  for (i in c(1:length(alphas))){
    DAG_list[[length(DAG_list)+1]] = as(grspdag.list[[i]], "matrix")
    DAG_names[length(DAG_names)+1] = paste(file,alphas[i],sep="_")
  }
}

get_children <-function(DAG, int_node) {
  DAG = t(DAG)
  ans = (DAG[,int_node] == 1)+0
  names(ans) = NULL
  return(ans)
}

get_children_skeleton <- function(DAG, int_node){
  DAG = t(DAG)
  ans = ((DAG[,int_node] == 1) | (DAG[int_node,] == 1))+0
  names(ans) = NULL
  return(ans)
}

# determine true positives, false positives, for given p-val cutoff, DAG, and interventional set
get_tp <- function(true_descendents, p_threshold, predicted_children){
  true_descendents = true_descendents>p_threshold
  true_positives = (predicted_children == 1) & (true_descendents == 1)
  return (sum(true_positives))
}
get_fp <- function(true_descendents, p_threshold, predicted_children){
  true_descendents = true_descendents>p_threshold
  false_positives = (predicted_children == 1) & (true_descendents == 0)
  return (sum(false_positives))
}
#,
for (sig in alphas){
print(sig)
predicted_children = true_descendents
#added for skeleton
#predicted_children_skeleton = true_descendents

# for each perturbation, find DAG name with same perturbation, get children
  for (i in 1:length(perturbations)) {
    perturb = perturbations[i]
    print(perturb)
    dag_idx = which(unlist(lapply(DAG_names, function(x) grepl(paste(perturb, '', sep = '_'), x))))
    sig_idx = which(unlist(lapply(DAG_names, function(x) grepl(gsub("\\.", "\\\\.", sig), x))))
    dag_idx = intersect(dag_idx, sig_idx)
    print(dag_idx)
    predicted_children[,i] = get_children(DAG_list[[dag_idx]], which(gene_names == perturb))
    # added for skeleton
    #predicted_children_skeleton[,i] = get_children_skeleton(DAG_list[[dag_idx]], which(gene_names == perturb))
  }
  
# changed for skeleton
true_positives = unlist(lapply(1:50/10, function(t) get_tp(true_descendents, t, predicted_children)))
#true_positives = unlist(lapply(1:50/10, function(t) get_tp(true_descendents, t, predicted_children_skeleton)))
false_positives = unlist(lapply(1:50/10, function(t) get_fp(true_descendents, t, predicted_children)))

results = data.frame(1, true_positives, false_positives)
save(results, file=paste('./accuracy_results/', paste(sig,'accuracy_results.rda', sep = '_'), sep = ''))

}


