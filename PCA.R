# Generate fake data set with 10 samples which have 100 genes in each
data.matrix <- matrix(nrow=100, ncol=10)

# Label the samples 
# (wt=wild type=normal,every day ),(ko=knock-out=missing a gene)
colnames(data.matrix) <- c(
  paste("wt",1:5,sep=""),
  paste("ko",1:5,sep=""))
rownames(data.matrix) <- paste("gene",1:100,sep="")

# Give fake genes to fake read counts
for(i in 1:100){
  wt.values <- rpois(5, lambda=sample(x=10:1000,size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000,size=1))
  
  data.matrix[i,] <- c(wt.values,ko.values);
}

# Take a look at fake data
head(data.matrix)

# Do the PCA (prcomp returns x, standart deviation and rotation)
pca <- prcomp(t(data.matrix),scale=TRUE)

# x are principal components, let's look at first and second ones
plot(pca$x[,1],pca$x[,2])

# To make clear these clusters
# Calculate variance of the data
pca.var <- pca$sdev^2

#Calculate percentages of variance
pca.var.per <- round(pca.var/sum(pca.var)*100,1)

# Plot percentages
barplot(pca.var.per,main="Scree Plot",xlab="Principal Component",
        ylab="Percent Variation")
# Comment: Since PC1 has almost entire variation, clusters have a big difference

# Load ggplot2 library for examining PCA
library(ggplot2)

# Format the data
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
# Look the data, one row per sample and each row has sample ID and coordinates
pca.data

# Now, make a graph with ggplot
ggplot(data=pca.data,aes(x=X,y=Y,label=Sample))+
  geom_text()+
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep=""))+
  theme_bw()+
  ggtitle("PCA graph")

# To determine which genes have the largest effect on where plotted samples in PCA plot, 
# extract loading scores
loading_scores <- pca$rotation[,1]

# left side is negative, right side is positive but we are interested in magnitude
gene_scores <- abs(loading_scores)

# Sort scores and get top 10 genes
gene_scores_ranked <- sort(gene_scores,decreasing=TRUE)
top_10_genes <- names(gene_scores_ranked[1:10])
top_10_genes

# Show sign of these genes
pca$rotation[top_10_genes,1]

# end.
