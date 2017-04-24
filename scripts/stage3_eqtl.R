library(MatrixEQTL)
library(ggplot2)

setwd("/Users/BC/Documents/clinical/eqtl/results/test/stage3/")

miR <- as.data.frame(read.table("stage3_mirna_matrix.txt", sep="\t", header=TRUE, row.names=1))
rna <- as.data.frame(read.table("stage3_rna_matrix.txt", sep="\t", header=TRUE,row.names=1))

mi <- as.matrix(miR);
mi2 <- mi[which(mi[,ncol(mi)]!='NA'),];
#x <- mi2[,1:(ncol(mi2))];
x <- cbind(as.matrix(rownames(mi)),mi);

r <- as.matrix(rna);
r2 <- r[which(r[,ncol(r)]!='NA'),];
y <- cbind(as.matrix(rownames(r)),r);
#y <- r2[,1:(ncol(r2))];


x1 <- NULL;
for(i in 1:nrow(x)){
  cat(i,',',sep="");
  index <- as.numeric(x[i,(2:ncol(x))]);
  if(sum(index!=0)){x1 <- rbind(x1,x[i,]);}
}
x2 <- NULL;
for(i in 1:nrow(x1)){
  cat(i,',',sep="");
  index <- as.numeric(x1[i,(2:ncol(x1))]);
  if(sum(index!=0) > 47){ x2 <- rbind(x2,x1[i,]);}
}

y1 <- NULL;
for(i in 1:nrow(y)){
  cat(i,',',sep="");
  index <- as.numeric(y[i,(2:ncol(y))]);
  if(sum(index!=0)){y1 <- rbind(y1,y[i,]);}
}
y2 <- NULL;
for(i in 1:nrow(y1)){
  cat(i,',',sep="");
  index <- as.numeric(y1[i,(2:ncol(y1))]);
  if(sum(index!=0) > 47){ y2 <- rbind(y2,y1[i,]);}
}
#x3 <- cbind(as.matrix(rownames(mi)),x2)
write.table(x2, file="stage3_mirna_80.txt",row.names=T, sep="\t")
write.table(y2, file="stage3_rna_80.txt",row.names=T, sep="\t")

miR <- as.data.frame(read.table("stage3_mirna_80.txt", sep="\t", header=TRUE, row.names=1))
rna <- as.data.frame(read.table("stage3_rna_80.txt", sep="\t", header=TRUE,row.names=1))

z.mirna <- scale(miR, center = TRUE, scale = TRUE)
z.rna <- scale(rna, center = TRUE, scale = TRUE)

write.table(z.mirna, file="stage3_mirna_sd.txt",row.names=T, sep="\t")
write.table(z.rna, file="stage3_rna_sd.txt",row.names=T, sep="\t")

library(preprocessCore)

#Quantile normalization of miRNA data
m <- as.matrix(read.table("stage3_mirna_sd.txt", header=TRUE, row.names = 1))
#dimnames(m) <- NULL
m <- as.matrix(read.table("stage3_rna_sd.txt", header=TRUE, row.names = 1))
#miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/hg19/
#Matrix must be numeric
new_m <- matrix(as.numeric(m))
#mode(new_m)

#Normalize the matrix
n <- normalize.quantiles(m,copy=F)
n<- as.matrix(n)
#Write to file
write.table(n, file = "quantile_stage3_rna.txt", sep = "\t")



# inverse_quantile_norm <- function(filename)
# {
#   library(MatrixEQTL)
#   
#   gene = SlicedData$new()$LoadFile(filename);
#   
#   
#   for( sl in 1:length(gene) ) {
#     mat = gene[[sl]];
#     mat = t(apply(mat, 1, rank, ties.method = "average"));
#     mat = qnorm(mat / (ncol(gene)+1));
#     gene[[sl]] = mat;
#   }
#   rm(sl, mat);
#   #
#   gene <- as.matrix(gene)
#   gene <- cbind(rownames(gene), gene)
#   colnames(gene)[1] = "ID"
#   write.table(gene,paste0(filename,'.qnorm.txt'),col.names=T,row.names=F,quote=F,sep="\t")
# }
# 
# inverse_quantile_norm('stage3_rna.txt')

miR <- as.data.frame(read.table("quantile_stage3_mirna.txt", sep="\t", header=TRUE, row.names=1))
rna <- as.data.frame(read.table("quantile_stage3_rna.txt", sep="\t", header=TRUE,row.names=1))

#miR <- log(miR2, base=exp(2))
old_colnames <- colnames(rna)

mirna_pc <- princomp(~., data=miR, center=F, scale=F, na.action=na.omit)
pc1_5 <- as.data.frame(mirna_pc$loadings[,1:7])
#rownames(pc1_5) <- oldsnpcolnames
colnames(pc1_5) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7")
pc1_5 <- t(pc1_5)
pc1_5<- cbind(rownames(pc1_5), pc1_5)
colnames(pc1_5)[1] = "ID"
write.table(pc1_5,paste0('stage3_mirna_pcs.txt'),col.names=T,row.names=F,quote=F,sep="\t")

pdf(paste0('stage3_mirna_screeplot.pdf'))
screeplot(mirna_pc, type='lines')
dev.off()
p1p2 <- cbind(mirna_pc$loadings[,1],mirna_pc$loadings[,2],mirna_pc$loadings[,3],mirna_pc$loadings[,4],mirna_pc$loadings[,5],mirna_pc$loadings[,6], mirna_pc$loadings[,7])
#
race <- read.table('race.txt',header=T,sep="\t",row.names=1,na.strings='NA',stringsAsFactors=T)
rownames(race) <- paste('X', rownames(race),sep='')
pc1pc2<-cbind.data.frame(p1p2,race$race[match(rownames(p1p2),rownames(race))])
colnames(pc1pc2) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7", 'Race')

pcaplot<-ggplot(pc1pc2,aes(x=PC1,y=PC2,colour=factor(Race)))+geom_point()+labs(x="PC1",y="PC2",colour="Race",title=paste0("stage 2 miRNA PC1 vs PC2"))
ggsave(paste0('stage3_race_pc1_pc2.pdf'),plot=pcaplot)

gender <- read.table(paste0('gender.txt'),header=T,sep="\t",row.names=1,check.names=F,stringsAsFactors=T)
#gender <- t(gender)
rownames(gender) <- paste('X', rownames(gender),sep='')
pc1pc2.gender<-cbind.data.frame(p1p2,gender$gender[match(rownames(p1p2),rownames(gender))])
colnames(pc1pc2.gender) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7",'Gender')
#
pcaplot<-ggplot(pc1pc2.gender,aes(x=PC1,y=PC2,colour=factor(Gender)))+geom_point()+labs(x="PC1",y="PC2",colour="Gender",title=paste0("stage 2 miRNA PC1 vs PC2"))
ggsave(paste0('stage3_mirna_pc1_pc2-gender.pdf'),plot=pcaplot)

#####PC Analysis gene expression#####
gene_pc <- princomp(~., data=rna, center=F, scale=F, na.action=na.omit)
gene_pc1_5 <- as.data.frame(gene_pc$loadings[,1:7])
rownames(gene_pc1_5) <- old_colnames
colnames(gene_pc1_5) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7")
gene_pc1_5 <- t(gene_pc1_5)
gene_pc1_5<- cbind(rownames(gene_pc1_5), gene_pc1_5)
colnames(gene_pc1_5)[1] = "ID"
write.table(gene_pc1_5,paste0('stage3_gene_pcs.txt'),col.names=T,row.names=F,quote=F,sep="\t")

pdf(paste0('stage3_gene_expression-screeplot.pdf'))
screeplot(gene_pc, type='lines')
dev.off()
gene_p1p2 <- cbind(gene_pc$loadings[,1],gene_pc$loadings[,2],gene_pc$loadings[,3],gene_pc$loadings[,4],gene_pc$loadings[,5],gene_pc$loadings[,6],gene_pc$loadings[,7])
#
gene_race <- read.table('race.txt',header=T,sep="\t",row.names=1,na.strings='NA',stringsAsFactors=T)
rownames(gene_race) <- paste('X', rownames(gene_race),sep='')
gene_pc1pc2<-cbind.data.frame(gene_p1p2,gene_race$race[match(rownames(gene_p1p2),rownames(gene_race))])
colnames(gene_pc1pc2) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7", 'Race')
                   
gene_pcaplot<-ggplot(gene_pc1pc2,aes(x=PC1,y=PC2,colour=factor(Race)))+geom_point()+labs(x="PC1",y="PC2",colour="Race",title=paste0("stage 2 gene PC1 vs PC2"))
ggsave(paste0('stage3_gene_pc1_pc2.pdf'),plot=gene_pcaplot)
                   
gene_gender <- read.table(paste0('gender.txt'),header=T,sep="\t", row.names=1,na.strings='NA',check.names=F,stringsAsFactors=T)                   #gene_gender <- t(gene_gender)
rownames(gene_gender) <- paste('X', rownames(gene_gender),sep='')

pc1pc2<-cbind.data.frame(p1p2,gene_gender$gender[match(rownames(p1p2),rownames(gene_gender))])
colnames(pc1pc2) <- c("PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7",'Gender')
#
gene_pcaplot<-ggplot(pc1pc2,aes(x=PC1,y=PC2,colour=factor(Gender)))+geom_point()+labs(x="PC1",y="PC2",colour="Gender",title=paste0(" stage 2 gene PC1 vs PC2"))
ggsave(paste0('stage3_gene_pc1_pc2-gender.pdf'),plot=gene_pcaplot)

p1p2 <- cbind(gene_pc$loadings[,1],gene_pc$loadings[,2])
pc1pc2<-cbind.data.frame(p1p2,race$race[match(rownames(p1p2),rownames(race))])
                   
pdf('stage3_gene-iqnorm-pc1pc2.pdf')
plot(pc1pc2[,1],pc1pc2[,2],col=pc1pc2[,3])
legend('bottomleft',legend=unique(pc1pc2[,3]),col=as.numeric(unique(pc1pc2[,3])),pch=1)
dev.off()

#PCA Covariates for Inverse Quantile Normed Expression
expr_pcs <- gene_pc$loadings[,1:10]
expr_pcs <- t(expr_pcs)
colnames(expr_pcs) <- old_colnames
                   
gender <- read.table('gender.txt',header=T,sep="\t", row.names=1, check.names=F,stringsAsFactors=T)
#gender <- t(gender)
rownames(gender) <- paste('x', rownames(gender),sep='')
pc1pc2<-cbind.data.frame(p1p2,gender[,1][match(rownames(p1p2),rownames(gender))])
                   
#cov_race <-cbind(gene_p1p2,gene_race$race[match(rownames(gene_p1p2),rownames(gene_race))])

library(preprocessCore)

#Quantile normalization of miRNA data
m <- as.matrix(read.table("stage3_mirna.txt", header=TRUE, row.names = 1))
#dimnames(m) <- NULL
m <- as.matrix(read.table("stage3_rna.txt", header=TRUE, row.names = 1))
#miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/hg19/
#Matrix must be numeric
new_m <- matrix(as.numeric(m))
#mode(new_m)

#Normalize the matrix
n <- normalize.quantiles(m,copy=F)
n<- as.matrix(n)
#Write to file
write.table(n, file = "quantile_stage3_rna_.txt", sep = "\t")




                   
useModel = modelLINEAR;
#cov <- as.matrix(read.table("stage3_covariates.txt",header=T,row.names=1, sep="\t"))

snps = SlicedData$new()$LoadFile("quantile_stage3_mirna.txt")#stage3_mirna_6.txt.qnorm.txt");                   #covariates_file = cov
cvrt = SlicedData$new()$LoadFile('stage3_cov.txt')
gene = SlicedData$new()$LoadFile("quantile_stage3_rna.txt")#stage3_rna_6.txt.qnorm.txt");
mirna_pos <- read.table("miRNA_positions.txt", header=T, sep="\t")
gene_pos <- read.table("gene_positions.txt", header=T, sep="\t")

snpspos <- read.table('miRNA_positions.txt', header = TRUE, stringsAsFactors = FALSE, "\t");
genepos <- read.table('gene_positions.txt', header = TRUE, stringsAsFactors = FALSE, "\t");


# p <- read.table("qnorm_stage_3-7_cis_meqtl_3_res.txt", header=T, sep="\t")
# p <- p[,5]
# fdr <- round(p.adjust(0.0162311323, "fdr"), 3)                   
                   
pvOutputThreshold_trans= 0;
pvOutputThreshold_cis = .05;
errorCovariance= numeric();
#pvOutputThreshold = pvOutputThreshold.cis
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = "quantile_trans_stage3.txt",
  output_file_name.cis = "quantile_cis_stage3_res.txt",
  #output_file_name = output_file_name.cis,  
  pvOutputThreshold =  pvOutputThreshold_trans,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvOutputThreshold.cis =  pvOutputThreshold_cis,
  snpspos <- snpspos,
  genepos <- genepos,
  cisDist <- 1e6,
  pvalue.hist <- "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);
                   
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected ',me$cis$neqtls,' local eQTLs:', '\n');
cat('Detected ',me$trans$neqtls,' distant eQTLs:', '\n');
                   
plot(me, pch = 16, cex = 0.7)
                   
                   