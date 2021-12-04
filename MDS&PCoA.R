# Load ggplot2 for graphs
library(ggplot2)

# Generate fake data
data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))
rownames(data.matrix) <- paste("gene", 1:100, sep="")
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
head(data.matrix)
dim(data.matrix)

# For reference, draw a PCA plot
pca <- prcomp(t(data.matrix), scale=TRUE, center=TRUE) 

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PCA Graph")

# Create MDS/PCoA plot
#Step 1: Create a distance matrix
distance.matrix <- dist(scale(t(data.matrix), center=TRUE, scale=TRUE),
                        method="euclidean")

#Step 2: Perform multi-dimensional scaling on dm
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

#Step 3: Calculate variation using eigen values
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

# Format the data for ggplot
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

# Plot
ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot using Euclidean distance")
# Commment: PCA and MDS are same since we use Euclidean.

# Choose a different metric (average of absolute values of log fold change)
# Calculate log2 values
log2.data.matrix <- log2(data.matrix)

# Create an empty distance matrix
log2.distance.matrix <- matrix(0,
                               nrow=ncol(log2.data.matrix),
                               ncol=ncol(log2.data.matrix),
                               dimnames=list(colnames(log2.data.matrix),
                                             colnames(log2.data.matrix)))

log2.distance.matrix

# Fill the matrix with average of the absolute values of log fold changes
for(i in 1:ncol(log2.distance.matrix)) {
  for(j in 1:i) {
    log2.distance.matrix[i, j] <-
      mean(abs(log2.data.matrix[,i] - log2.data.matrix[,j]))
  }
}
log2.distance.matrix
# Comment: Because fold matrix is symmetrical, use low triangle

# Perfom scaling
mds.stuff <- cmdscale(as.dist(log2.distance.matrix),
                      eig=TRUE,
                      x.ret=TRUE)

# Calculate the percentage of variation that each MDS axis accounts for
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

# Format the data
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

# Plot
ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot using avg(logFC) as the distance")

# end.