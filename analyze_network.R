rm(list=ls())
library(igraph)

# functions

# do majarity cut and calculate correlation matrix
# otu, OTU table
# keep, number of samples an OTU has to occur in to be included in corrlation calculation
maj_cut_cor = function(otu, keep) {
  maj_id = apply(otu, 1, FUN=function(x){
    if (sum(x[]!=0) >= keep)
      return(1)
    else
      return(0)
  })
  maj = otu[which(maj_id[] == 1),]
  
  sp = cor(t(maj), method="spearman")
  
  #	pearson = cor(t(maj), method = "pearson")
  
  #	clr = apply((maj + 1), 2, function(xc){
  #            log(xc, 2) - sum(log(xc, 2))/length(xc)
  #	})
  #	clr_pearson = cor(t(clr), method="pearson")
}

# fit log power law
fit_power_law = function(graph) {
  # calculate degree
  d = centr_degree(graph)$res
  dd = degree_distribution(graph)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  # print(paste("Alpha =", round(alpha, 3)))
  return(round(R.square, 3))
  # plot
  # plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", col = 1, main = "Degree Distribution")
  # curve(power.law.fit, col = "red", add = T, n = length(d))
}


setwd("/Users/maggieyuan/Desktop/")

# input OTU table
otutable = read.table("sample.txt", header=T, row.names=1)

# majority cut and correlation calculation
correlation_matrix = maj_cut_cor(otu=otutable, keep=10)

cutoff = 0.6 # correlation cutoff, decided by RMT

# calculate adjacency matrix
transition_matrix = correlation_matrix
transition_matrix[which(abs(transition_matrix) <= cutoff)] <- 0
  
adjacancy_matrix = transition_matrix
adjacancy_matrix[adjacancy_matrix != 0] <- 1
  
# construct graph
my_graph = graph_from_adjacency_matrix(as.matrix(adjacancy_matrix), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)
  
# global properties of the graph
conn_nodes = gorder(my_graph) # nodes
links = ecount(my_graph) # links
r2 = fit_power_law(my_graph) # r2 of power law
avgK = mean(centr_degree(my_graph)$res) # avrerage degree
avgCC = transitivity(my_graph, type="average", isolates="zero") # average connectivity
GD = mean_distance(my_graph, directed=F, unconnected=T) # geodesic distance
gd = cluster_fast_greedy(my_graph) # greedy module separation
modules = length(gd_no_iso) # number of modules
M = modularity(gd_no_iso) # modularity