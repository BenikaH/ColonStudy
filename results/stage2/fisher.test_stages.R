
#Test the hypothesis that the stage 2 and stage 3 networks are independent of each other.
#If we obtain a low p-value (<<.05), we can reject the null hypothesis.

#Load in the data
stage2 <- as.data.frame(read.csv("/Users/BC/Desktop/GeNets_Analysis_Stage2/Stage_2_WithCandidates/NetworkEdges.csv", header=T))
stage3 <- as.data.frame(read.csv("/Users/BC/Desktop/GeNets_Analysis_Stage3/Stage_3_WithCandidates/NetworkEdges3.csv", header=T))
stage2 <- stage2[,1:2]
stage3 <- stage3[,1:2]

#Combine both data frames
all_edges <- rbind(stage2,stage3)

#Find all possible edge combinations for the universe
a <- all_edges[,1]
b <- all_edges[,2]
df = expand.grid(a = a, b = b)
sorted <- df[order(df$a), ]
num_all_possible <- nrow(sorted)

#Count the number of edges in common for both stages
all_edges$count <- 1
t <- aggregate(count ~ ., all_edges, FUN = sum)
common_edges <- subset(t,count=="2")
num_in_both <- nrow(common_edges) 

#Find the edges in stage3, but not in stage 2
x <- rbind(stage3, stage2)
not_in_stage2 <- x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(stage3), ] 
stage3_unique <- nrow(not_in_stage2)

#Find the edges in stage 2, but not in stage 3
x2 <- rbind(stage2, stage3)
not_in_stage3 <- x2[! duplicated(x2, fromLast=TRUE) & seq(nrow(x2)) <= nrow(stage2), ] 
stage2_unique <- nrow(not_in_stage3)

#Find the edges that are not found in both
#not_in_either <- sorted[! duplicated(sorted, fromLast=TRUE) & seq(nrow(sorted)) <= nrow(rbind(stage2,stage3)), ]
#not_in_either <- sorted[! duplicated(all_edges), ]
not_in_either <- num_all_possible - stage2_unique - stage3_unique

#Perform Fisher Exact test
c <- matrix(c(97,202,439,696584), nrow=2)
f <- fisher.test(c)
