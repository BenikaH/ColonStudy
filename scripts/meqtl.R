library(MatrixEQTL)

#Preprocessing
#####Load in the data######
mi <- read.table("mirna2_matrix.txt",header=T, row.names=1);#mirna2_matrix.txt
mi <- as.matrix(mi);
mi2 <- mi[which(mi[,ncol(mi)]!='NA'),];
x <- cbind(as.matrix(rownames(mi)),mi);
x <- mi2[,2:(ncol(mi2))];

r <- read.table("rna_matrix.txt",header=T, row.names=1);#mirna2_matrix.txt
r <- as.matrix(r);
r2 <- r[which(r[,ncol(r)]!='NA'),];
#y <- cbind(as.matrix(rownames(r)),r);
y <- r2[,2:(ncol(r2))];
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
  if(sum(index!=0) > 22){ x2 <- rbind(x2,x1[i,]);}
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
  if(sum(index!=0) > 22){ y2 <- rbind(y2,y1[i,]);}
}
#####Write matrix to file#####
write.table(x2, file="mirna_norm_matrix.txt",row.names=T, sep="\t")
write.table(y2, file="norm_rnas.txt",row.names=T, sep="\t")

#####Load in the new data######
rna <- as.matrix(read.table("norm_rnas.txt", sep="\t", header=TRUE, row.names=1))
miR <- as.matrix(read.table("rna_norm_matrix.txt", sep="\t", header=TRUE,row.names=1))

#####Normalize the data#####
for(sl in 1:nrow(miR)) {
  mat = as.matrix(miR[sl,]);
  mat = t(apply(mat, 2, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(miR)+1));
  miR[sl,] = mat;
}
#rm(sl, mat);
for(sl2 in 1:nrow(rna)) {
  mat = as.matrix(rna[sl2,]);
  mat = t(apply(mat, 2, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(rna)+1));
  rna[sl2,] = mat;
}

phe <- read.table("clinical.txt", header=T)






useModel = modelLINEAR;

## Load miRNA data

mirna = SlicedData$new();
mirna$fileDelimiter = "\t";      # the TAB character
mirna$fileOmitCharacters = "NA"; # denote missing values;
mirna$fileSkipRows = 1;          # one row of column labels
mirna$fileSkipColumns = 1;       # one column of row labels
mirna$fileSliceSize = 2000;      # read file in slices of 2,000 rows
mirna$LoadFile(mirna_file_name);



## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);



## Run the analysis

me = Matrix_eQTL_engine(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
