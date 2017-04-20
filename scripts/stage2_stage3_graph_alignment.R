library(GraphAlignment)
library(igraph)
setwd("/Users/BC/Documents/spinglass/results/matching_communities/")


#temp = list.files(pattern="*.csv")
#myfiles = lapply(temp, read.delim)

stage2 <- read.csv("/Users/BC/Desktop/GeNets_Analysis_Stage2/Stage_2_WithCandidates/NetworkEdges.csv", header=T)

stage3 <- read.csv("/Users/BC/Desktop/GeNets_Analysis_Stage3/Stage_3_WithCandidates/NetworkEdges3.csv", header=T)

y2=as.matrix(stage2[,1:2])
y3=as.matrix(stage3[,1:2])


g2 <- graph.edgelist(y2, directed=FALSE)
g3 <- graph.edgelist(y3, directed=FALSE)

g_2 <- get.adjacency(g2)
g_3 <- get.adjacency(g3)

g_2.1 <- graph.adjacency(g_2)
g_3.1 <- graph.adjacency(g_3)

#Calculate the node similarity between two graphs
g_sim <- graph.intersection(g_2.1,g_3.1, byname = "auto", keep.all.vertices = FALSE)

adj_sim <- get.adjacency(g_sim, type="both")

g<- graph.adjacency(adj_sim,mode = "undirected",weighted = T)

edge_ara <- get.data.frame(g,what = "edges")

edge_ara<- edge_ara[order(abs(edge_ara$weight),decreasing = T),]

write.csv(edge_ara, "comm2_node_sim.csv")

#calculate all possible edges,



