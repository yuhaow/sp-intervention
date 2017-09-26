#set random seed
set.seed(1)

#get the data
samples <- read.csv('sachs.csv', sep = '\t');

#GIES data
experiments <- c(3, 4, 5, 6, 7)
t.list <- list(integer(0), c(7), c(9), c(4), c(2), c(5))
data.list <- lapply(experiments, function(t) log(samples[which(samples[['experiment']]==t), 1:11]))
data.list <- c(list(log(samples[which(samples[['experiment']]==1 | samples[['experiment']]==2), 1:11])), data.list)
title <- "safe_choice.rda"
save(data.list, t.list, title, file=title)
