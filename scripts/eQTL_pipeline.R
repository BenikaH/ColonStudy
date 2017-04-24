library(MatrixEQTL)
library(ggplot2)

setwd("/Users/BC/Documents/clinical/eqtl/")

#Preprocessing
#####Load in the data######

mi <- read.table("miRNA_matrix_condensed.txt",header=T, sep="\t",row.names=1)

mi <- as.matrix(mi);
mi2 <- mi[which(mi[,ncol(mi)]!='NA'),];
x <- cbind(as.matrix(rownames(mi)),mi);
x <- mi2[,2:(ncol(mi2))];

r <- read.table("RNAseq_matrix_condensed.txt",header=T,sep="\t", row.names=1);#mirna2_matrix.txt
r <- as.matrix(r);
r2 <- r[which(r[,ncol(r)]!='NA'),];
y <- cbind(as.matrix(rownames(r)),r);
#y <- r2[,2:(ncol(r2))];
#####Filter out those with 90% or more zeros#####
#x1[i,(2:ncol(x))], and combine x1 <- rbind(x1,x[i,])
#x <- rbind(t(as.matrix(dimnames(x)[[2]])),x);

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
  if(sum(index!=0) > 14){ x2 <- rbind(x2,x1[i,]);}
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
  if(sum(index!=0) > 14){ y2 <- rbind(y2,y1[i,]);}
}
#####Write matrix to file#####
write.table(rna, file="quantile_rna.txt",row.names=T, sep="\t")
write.table(miR, file="quantile_mirna.txt",row.names=T, sep="\t")

#####Load in the new data######
rna.1 <- as.matrix(read.table("quantile_rna.txt", sep="\t", header=TRUE, row.names=1))
miR.1 <- as.matrix(read.table("quantile_mirna.txt", sep="\t", header=TRUE,row.names=1))
miR.1 <-as.data.frame(miR.1)
rna.1 <- as.data.frame(rna.1)
#####Normalize the data#####
for(sl in 1:nrow(miR)) {
  mat = as.matrix(miR[sl,]);
  mat = t(apply(mat, 2, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(miR)+1));
  miR[sl,] = mat;
}
#rm(sl, mat);
gene <- as.matrix(gene)
gene <- cbind(rownames(gene), gene)
colnames(gene)[1] = "ID"
write.table(gene,paste0(filename,'.qnorm'),col.names=T,row.names=F,quote=F,sep="\t")

for(sl2 in 1:nrow(rna)) {
  mat = as.matrix(rna[sl2,]);
  mat = t(apply(mat, 2, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(rna)+1));
  rna[sl2,] = mat;
}
gene <- as.matrix(rna)
gene <- cbind(rownames(gene), gene)
#	colnames(gene)[1] = "ID"
#	write.table(gene,paste0(filename,'.qnorm'),col.names=T,row.names=F,quote=F,sep="\t")
rna.1 <- as.matrix(read.table("rna.txt", sep="\t", header=TRUE, row.names=1))
miR.1 <- as.matrix(read.table("mirna.txt", sep="\t", header=TRUE,row.names=1))
#####PCA Analysis miRNA#####
old_colnames <- colnames(rna.1)
mirna_pc <- princomp(~., data=miR.1, center=F, scale=F, na.action=na.omit)
pc1_5 <- as.data.frame(mirna_pc$loadings[,1:5])
#rownames(pc1_5) <- oldsnpcolnames
colnames(pc1_5) <- c("PC1","PC2", "PC3", "PC4", "PC5")
pc1_5 <- t(pc1_5)
pc1_5<- cbind(rownames(pc1_5), pc1_5)
colnames(pc1_5)[1] = "ID"
write.table(pc1_5,paste0('data/mirna_pc_covariates_gender'),col.names=T,row.names=F,quote=F,sep="\t")

pdf(paste0('mirna_genotype-screeplot.pdf'))
screeplot(mirna_pc, type='lines')
dev.off()
p1p2 <- cbind(mirna_pc$loadings[,1],mirna_pc$loadings[,2],mirna_pc$loadings[,3],mirna_pc$loadings[,4],mirna_pc$loadings[,5])
#
race <- read.table('data/covar_race.txt',header=T,sep="\t",row.names=1,na.strings='NA',stringsAsFactors=T)
rownames(race) <- paste('X', rownames(race),sep='')
pc1pc2<-cbind.data.frame(p1p2,race$race[match(rownames(p1p2),rownames(race))])
colnames(pc1pc2) <- c('PC1','PC2','PC3', 'PC4', 'PC5', 'Race')

pcaplot<-ggplot(pc1pc2,aes(x=PC1,y=PC2,colour=factor(Race)))+geom_point()+labs(x="PC2",y="PC3",colour="Race",title=paste0(" PC2 vs PC3"))
ggsave(paste0('_pc2_pc3.pdf'),plot=pcaplot)

gender <- read.table(paste0('data/covariates'),header=T,check.names=F,stringsAsFactors=T)
gender <- t(gender)
rownames(gender) <- paste('X', rownames(gender),sep='')
pc1pc2.gender<-cbind.data.frame(p1p2,gender[,1][match(rownames(p1p2),rownames(gender))])
colnames(pc1pc2.gender) <- c('PC1','PC2', 'PC3', 'PC4', 'PC5','Gender')
#
pcaplot<-ggplot(pc1pc2.gender,aes(x=PC2,y=PC3,colour=factor(Gender)))+geom_point()+labs(x="PC2",y="PC3",colour="Gender",title=paste0("miRNA PC2 vs PC3"))
ggsave(paste0('mirna_pc2_pc3-gender.pdf'),plot=pcaplot)

#####PC Analysis gene expression#####
gene_pc <- princomp(~., data=rna.1, center=F, scale=F, na.action=na.omit)
gene_pc1_5 <- as.data.frame(gene_pc$loadings[,1:5])
#  rownames(pc1_2) <- oldsnpcolnames
colnames(gene_pc1_5) <- c("PC1","PC2", "PC3", "PC4", "PC5")
gene_pc1_5 <- t(gene_pc1_5)
gene_pc1_5<- cbind(rownames(gene_pc1_5), gene_pc1_5)
colnames(gene_pc1_5)[1] = "ID"
write.table(gene_pc1_5,paste0('data/gene_pc_covariates'),col.names=T,row.names=F,quote=F,sep="\t")

pdf(paste0('gene_expression-screeplot.pdf'))
screeplot(gene_pc, type='lines')
dev.off()
gene_p1p2 <- cbind(gene_pc$loadings[,1],gene_pc$loadings[,2],gene_pc$loadings[,3],gene_pc$loadings[,4],gene_pc$loadings[,5])
#
gene_race <- read.table('data/covar_race.txt',header=T,sep="\t",row.names=1,na.strings='NA',stringsAsFactors=T)
rownames(gene_race) <- paste('X', rownames(gene_race),sep='')
gene_pc1pc2<-cbind.data.frame(gene_p1p2,gene_race$race[match(rownames(gene_p1p2),rownames(gene_race))])
colnames(gene_pc1pc2) <- c('PC1','PC2','PC3', 'PC4', 'PC5', 'Race')

gene_pcaplot<-ggplot(gene_pc1pc2,aes(x=PC2,y=PC3,colour=factor(Race)))+geom_point()+labs(x="PC2",y="PC3",colour="Race",title=paste0(" PC2 vs PC3"))
ggsave(paste0('gene_pc2_pc3.pdf'),plot=gene_pcaplot)

gene_gender <- read.table(paste0('data/covariates'),header=T,check.names=F,stringsAsFactors=T)
gene_gender <- t(gene_gender)
rownames(gene_gender) <- paste('X', rownames(gene_gender),sep='')
pc1pc2<-cbind.data.frame(p1p2,gene_gender[,1][match(rownames(p1p2),rownames(gene_gender))])
colnames(pc1pc2) <- c('PC1','PC2', 'PC3', 'PC4', 'PC5','Gender')
#
gene_pcaplot<-ggplot(pc1pc2,aes(x=PC2,y=PC3,colour=factor(Gender)))+geom_point()+labs(x="PC2",y="PC3",colour="Gender",title=paste0(" PC2 vs PC3"))
ggsave(paste0('gene_pc2_pc3-gender.pdf'),plot=gene_pcaplot)

p1p2 <- cbind(gene_pc$loadings[,1],gene_pc$loadings[,2])
pc1pc2<-cbind.data.frame(p1p2,race$race[match(rownames(p1p2),rownames(race))])

pdf('gene-iqnorm-pc1pc2.pdf')
plot(pc1pc2[,1],pc1pc2[,2],col=pc1pc2[,3])
legend('bottomleft',legend=unique(pc1pc2[,3]),col=as.numeric(unique(pc1pc2[,3])),pch=1)
dev.off()

#PCA Covariates for Inverse Quantile Normed Expression
expr_pcs <- gene_pc$loadings[,1:10]
expr_pcs <- t(expr_pcs)
colnames(expr_pcs) <- old_colnames

gender <- read.table('data/covariates',header=T,check.names=F,stringsAsFactors=T)
gender <- t(gender)
rownames(gender) <- paste('x', rownames(gender),sep='')
pc1pc2<-cbind.data.frame(p1p2,gender[,1][match(rownames(p1p2),rownames(gender))])


useModel = modelLINEAR;
cov <- as.matrix(read.table("data/gene_pc_covariates",header=T,row.names=1, sep="\t"))

snps = SlicedData$new()$LoadFile("mirna.txt");
covariates_file = cov
cvrt = SlicedData$new()$LoadFile('data/gene_pc_covariates')
gene = SlicedData$new()$LoadFile("rna.txt");
mirna_pos <- read.table("miRNA_positions.txt", header=T, sep="\t")
gene_pos <- read.table("data/gene_positions.txt", header=T, sep="\t")

snpspos <- read.table('miRNA_positions.txt', header = TRUE, stringsAsFactors = FALSE, "\t");
genepos <- read.table('data/gene_positions.txt', header = TRUE, stringsAsFactors = FALSE, "\t");



## Run the analysis
## Load miRNA data
getmiRNA <- function(sep, missing, header, rownames, mirna_filename)
{
  if(!file.exists('mirna.txt')){
    stop(cat(mirna_filename, 'no such file', "\n"));
  }
  mirna = SlicedData$new();
  mirna$fileDelimiter = "\t";      # the TAB character
  mirna$fileOmitCharacters = "NA"; # denote missing values;
  mirna$fileSkipRows = 1;          # one row of column labels
  mirna$fileSkipColumns = 1;       # one column of row labels
  mirna$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  mirna$LoadFile('mirna.txt');
}


## Load gene expression data
getExpr <- function(sep, missing, header, rownames, expr_filename)
{
  if (!file.exists('rna.txt'))
  {
    stop(cat('rna.txt', "no such file", "\n"))
  }
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile('rna.txt');
}

getCovariates <- function(sep, missing, header, rownames, cvrt_filename)
{
  skipRows <- toSkip(header);
  skipCols <- toSkip(rownames);
  cvrt <- SlicedData$new();
  cvrt$fileDelimiter <- sep;
  cvrt$fileOmitCharacters <- missing;
  cvrt$fileSkipRows <- skipRows;
  if (length('data/covariates') == 0)
  {
    return(cvrt);
  }
  if (!file.exists('data/covariates'))
  {
    stop(cat('data/covariates',' does not exist','\n'));
  }
  cvrt$LoadFile('data/covariates');
  return(cvrt)
}


output_file_name_cis <- cis_output_file
if (trans_output_file!="")
{
  output_file_name_tra <- trans_output_file;
}
else
{
  output_file_name_tra <- "onlyTRANSresults.txt"
}

pvOutputThreshold_trans= 0;
pvOutputThreshold_cis = 5e-2;
errorCovariance= numeric();
#pvOutputThreshold = pvOutputThreshold.cis
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name <- "trans_only.txt",
  output_file_name.cis = "cis_meqtl_res.txt",
  #output_file_name = output_file_name.cis,  
  pvOutputThreshold =  pvOutputThreshold_trans,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvOutputThreshold.cis =  pvOutputThreshold_cis,
  snpspos <- snpspos,
  genepos <- genepos,
  cisDist <- 100000,
  pvalue.hist <- "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected ',me$cis$neqtls,' local eQTLs:', '\n');
cat('Detected ',me$trans$neqtls,' distant eQTLs:', '\n');


    
    
    ## Plot the Q-Q plot of local and distant p-values
if (qq!="")
{
  pdf(qq)
  plot(me)
  dev.off()
    }
  else
    {
  plot(me)
    }
    
  return(me) #nocov end
  }
