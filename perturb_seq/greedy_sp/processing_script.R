rm(list=ls())
files = list.files('./accuracy_results')

greedy_sp_pvals = as.numeric(gsub("_accuracy_results\\.rda", "", files))

sorted = sort.int(greedy_sp_pvals, index.return = T)
greedy_sp_pvals = sorted$x
sorting_idx = sorted$ix

# reorder files
files = files[sorting_idx]

true_positives = list()
for (file in files) {
  load(paste('./accuracy_results/', file, sep = ''))
  true_positives[[length(true_positives)+1]] = results$true_positives
}

true_positives = do.call(cbind, true_positives)
rownames(true_positives) = 1:50/10
colnames(true_positives) = greedy_sp_pvals

false_positives = list()
for (file in files) {
  load(paste('./accuracy_results/', file, sep = ''))
  false_positives[[length(false_positives)+1]] = results$false_positives
}
false_positives = do.call(cbind, false_positives)
rownames(false_positives) = 1:50/10
colnames(false_positives) = greedy_sp_pvals

load('proportion_true.rda')

save(proportion_true, true_positives, false_positives, file = 'greedy_sp_summary.rda')
