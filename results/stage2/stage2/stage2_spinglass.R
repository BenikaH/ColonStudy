library(igraph)
library(GraphAlignment)
setwd("/Users/BC/Desktop/GeNets_Analysis_Stage2/stage2/")
genes <- read.csv("stage2_connected.csv", header=T)
g=as.matrix(genes)
y <- graph.edgelist(g, directed=FALSE)#, weights=E(y)$weight)
#snps <- read.table("~/Documents/NeXTProject/combined_MTLasso_mEQTL_eqtls.txt", header=TRUE, sep="\t")
#snps <- as.data.frame(snps)
#snps <- snps$SNP
#snps
y <- simplify(y)

bad.vs<-V(y)[degree(y) == 1] 

# remove isolated nodes
y <-delete.vertices(y, bad.vs)
is.connected(y)
diameter(y, directed = FALSE, unconnected = FALSE, weights = NULL)
get.diameter(y, directed = FALSE, unconnected = FALSE, weights = NULL)

sg1 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg2 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg3 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg4 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg5 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg6 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg7 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg8 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg1 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)
sg10 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)




edgelist <- tapply(seq_along(membership(sg1)), membership(sg1), function(xx) xx)
comList <- tapply(membership(sg1), membership(sg1), names)

length(comList)                         ## number of communities
comsize <- sapply(comList, length)
comsize
bigComIndex <- which(comsize>7)
bigComIndex
#bigComIndex <- which(neighborhood(graph=y,order=7,nodes=snps))
h <- comList[bigComIndex]
h
coms3 <- edgelist[bigComIndex[[2]]]    ## pull out first community
coms3
#coms4 <- edgelist[bigComIndex[[3]]]    ## pull out third community
comGraphBig1 <- induced.subgraph(y, coms3[[1]])
#ly1 <- layout[unlist(coms3[[2]]),] 
cc <-get.edgelist(comGraphBig1)
cc
write.csv(cc, "results_default/stage2_comm2.csv")

