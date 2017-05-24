# Calculates the stemminess of a tree as the proportion of the total tree length accounted for by internal branches.

stemmy <- function(tre) sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))], na.rm = T) / sum(tre$edge.length, na.rm = T)