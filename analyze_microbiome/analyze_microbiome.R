## Load the data into a R 
ibd_taxa <- read.csv('ibd_taxa.xls', row.names=1) # taxon x sample matrix
ibd_lineages <- read.csv('ibd_lineages.xls', row.names=1) # lineage x taxon matrix
ibd_metadata <- read.csv('ibd_metadata.xls', row.names=1) # metadata x sample data matrix

## Viewing the data
head(ibd_taxa)

head(ibd_lineages)

head(ibd_metadata)

## Generating a toy count matrix
A <- rmultinom(5, size=100, prob=c(0.2, 0.4, 0.05, 0.02, 0.15, 0.13, 0.01, 0.04))

B <- rmultinom(5, size=100, prob=c(0.6, 0.25, 0, 0.04, 0.02, 0.06, 0.02, 0))

counts <- cbind(A, B)

groups <- c(rep('A', 5), rep('B', 5))

## bar plot
barplot(counts)

counts <- counts[order(rowMeans(counts)),]
barplot(counts)

## Make a bar plot part 1
barplot(as.matrix(ibd_taxa))

# generate number of colours equal to number of phyla
colours <- rainbow(length(unique(ibd_lineages$Phylum)))
color.indices <- match(ibd_lineages$Phylum, unique(ibd_lineages$Phylum))

# generate a vector with color for each taxon
colvec <- colours[color.indices]

# number of phyla
length(unique((ibd_lineages$Phylum)))

## Make a barplot part 2
barplot(as.matrix(ibd_taxa), col=colvec)

legend('topright', fill=colours,legend=unique(ibd_lineages$Phylum))

par(xpd=TRUE, mfrow = c(2,1), mar=c(1, 1, 1, 1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topright", fill=colours,legend=unique(ibd_lineages$Phylum), cex=.8)
barplot(as.matrix(ibd_taxa), col=colvec,  xaxt='n', yaxt='n', ylim=c(0, 1))

dev.off() # restore default plotting area

## make the labels with species coloured by their taxonomic class
length(unique(ibd_lineages$Class))

colours = rainbow(length(unique(ibd_lineages$Class)))

color.indices = match(ibd_lineages$Class, unique(ibd_lineages$Class))

colvec = colours[color.indices]

par(xpd=TRUE, mfrow = c(2,1), mar=c(1, 1, 1, 1))

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

legend("topright", fill=colours,legend=unique(ibd_lineages$Class), cex=.8)

barplot(as.matrix(ibd_taxa), col=colvec, xaxt='n', yaxt='n', ylim=c(0, 1))

dev.off()

## Making nicer plots
library(ggplot2)
library(reshape2)
long_data <- cbind(ibd_taxa, ibd_lineages)
long_data <- reshape2::melt(long_data)

ggplot(data=long_data, aes(x=variable, y=value, fill=Phylum))+
geom_bar(position='stack', stat='identity')+
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x='Sample', y='Relative abundance')
  
## Calculating richness
toy_richness <- colSums(counts != 0)

# Is the richness significantly different?
# first define a vector with group names A and B per sample
groups <- c(rep('A', 5), rep('B', 5))
wilcox.test(toy_richness~groups)

# The Chao 1 estimator
# calculate number of singletons
singletons <- colSums(counts == 1)
doubletons <- colSums(counts == 2)
rares <- singletons / (2*doubletons)
rares[doubletons == 0] <- 0

chao1 <- toy_richness + rares


## The Shannon index
# convert the toy matrix to a relative abundance matrix, implement normalization
sums <- apply(counts, 2, sum) # get the column sums (sample sums)
norm_counts <- counts # generate new data frame for normalized data
for (i in 1:ncol(counts)){ # divide per sample by sample total
  norm_counts [,i] <- counts[,i]/sums[i]
}


shannon_div <- function(vector){
  vector <- vector*log(vector)
  # the vector has NA values if the species proportions are 0
  vectorsum <- sum(na.omit(vector))
  return((-1)*vectorsum)
}

shannon_diversities <- apply(norm_counts, 2, shannon_div)

#  p-value for a Wilcoxon ranked-sum test on the Shannon diversity of the simulated data
wilcox.test(shannon_diversities~c(rep('A', 5), rep('B', 5)))
            
## The Simpson index           
simpson_indices <- colSums(norm_counts^2)
inverse_simpson_indices <- 1 / simpson_indices    


## rarefy count data
rarefaction <- function(count.vector=c(),depth=100){
  probab=count.vector/sum(count.vector)
  return(rmultinom(1,size=depth,prob=probab))
}
# we rarefy an example count vector
count.vector=rmultinom(n=1,size=150,prob=c(0.5,rep(0.1,5)))
rarefaction(count.vector)

## Diversity analysis
plot(shannon~Age, data=ibd_metadata, ylab='Shannon')

cor(ibd_metadata$shannon, ibd_metadata$Age, method=c("spearman"))
cor.test(ibd_metadata$shannon, ibd_metadata$Age, method=c("spearman"))

## Calculating the Jaccard distance
jaccard_sets <- function(a, b){
  taxa <- as.character(c(1:length(a)))
  aset <- taxa[a != 0]
  bset <- taxa[b != 0]
  return(1 - length(intersect(aset, bset))/length(union(aset, bset)))
}
jaccard_sets(counts[,1], counts[,10])


## calculate the Bray-Curtis dissimilarity
1-(2 * sum(apply(counts,1,min)))/ sum(counts)


## Carrying out the ordination
library(vegan)
bray <- vegdist(t(norm_counts))

pcoa.res.sim <- capscale(t(norm_counts)~1, distance='bray', na.action='na.omit')

eigenvalues <- eigenvals(pcoa.res.sim)

# What percentage of variation is contained in the first eigenvector?
eigenvalues[[1]]/sum(eigenvalues) * 100

## Plotting the axes
plot(pcoa.res.sim)

plot(pcoa.res.sim$CA$u[,c(1,2)],
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

## colour the samples
colvec = c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))

plot(pcoa.res.sim$CA$u[,c(1,2)], col=colvec,xlab='PCoA1', ylab='PCoA2')


## metaMDS
nmda.res.sim <- metaMDS(t(norm_counts), distance="bray", k=2, trymax=50)
stressplot(nmda.res.sim)

## Plotting the nMDS analysis
plot(nmda.res.sim)

## Colouring the nMDS
colvec = c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))
ordiplot(nmda.res.sim,type="n")
orditorp(nmda.res.sim, display="sites", col=colvec)
points(nmda.res.sim, display="species", col='green')

## PCoA for IBD data
pcoa.res.ibd <- capscale(t(ibd_taxa)~1, distance='bray', na.action='na.omit')
eigenvalues <- eigenvals(pcoa.res.ibd)
plot(pcoa.res.ibd$CA$u[,c(1,2)],
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

## Colouring the PCoA
colvec = rainbow(length(unique(ibd_metadata$Diagnosis)))[as.numeric(as.factor(ibd_metadata$Diagnosis))]
plot(pcoa.res.ibd$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))),
     ylim=c(-0.25, 0.35))
legend("topleft", legend=unique(ibd_metadata$Diagnosis), col=unique(colvec), cex=1, pch=1)

col_pal = colorRampPalette(c('red','blue'))(length(unique(ibd_metadata$Age)))
colvec = col_pal[ibd_metadata$Age]
plot(pcoa.res.ibd$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

## CnMDS for IBD data
colvec = rainbow(length(unique(ibd_metadata$Diagnosis)))[as.numeric(as.factor(ibd_metadata$Diagnosis))]
nmda.res.ibd <- metaMDS(t(ibd_taxa), distance="bray", k=2, trymax=50)
ordiplot(nmda.res.ibd,type="n")
points(nmda.res.ibd, display="sites", col=colvec)

## Envfit
# run envfit on simulation
groupvec <- c(rep("A", 5), rep("B", 5))
sim_meta <- data.frame(group=groupvec)
pcoa.res.sim <- capscale(t(norm_counts)~1, distance='bray', na.action='na.omit')
ef.sim = envfit(pcoa.res.sim, sim_meta, perm=1000, choices=c(1,2), na.rm=TRUE)

# Summarizing envfit results
print(ef.sim)

## Plotting envfit
plot(pcoa.res.sim)
plot(ef.sim)

## Visualizing envfit results
ef.sim$factors$centroids

colvec <- c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))
sim_ev <- eigenvals(pcoa.res.sim)
plot(pcoa.res.sim$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(sim_ev [[1]],2))), ylab=paste(c('PCoA2', round(sim_ev [[2]],2))))

x <- ef.sim$factors$centroids[,1]
y <- ef.sim$factors$centroids[,2]
labels <- rownames(ef.sim$factors$centroids)
text(x=x, y=y, labels=labels)

x <- ef.sim$factors$centroids[,1]*0.6
y <- ef.sim$factors$centroids[,2]*0.6

colvec <- c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))
sim_ev <- eigenvals(pcoa.res.sim)
plot(pcoa.res.sim$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(sim_ev [[1]],2))), ylab=paste(c('PCoA2', round(sim_ev [[2]],2))))
text(x=ef.sim$factors$centroids[,1]*0.6, y=ef.sim$factors$centroids[,2]*0.6, labels=rownames(ef.sim$factors$centroids))

x <- ef.sim$factors$centroids[,1]*0.6
y <- ef.sim$factors$centroids[,2]*0.6
x0 <- rep(0, length(x))
y0 <- rep(0, length(y))
arrows(x1=x, y1=y, x0=x0, y0=y0)
text(x=x, y=y, labels=labels, pos=3)

## IBD envfit
pcoa.res.ibd = capscale(t(ibd_taxa)~1, distance='bray', na.action='na.omit')
ibd_metadata2=ibd_metadata[,2:ncol(ibd_metadata)]
ef.ibd = envfit(pcoa.res.ibd, ibd_metadata2, perm=1000, choices=c(1,2), na.rm=TRUE)
colvec = rainbow(length(unique(ibd_metadata$Diagnosis)))[as.numeric(as.factor(ibd_metadata$Diagnosis))]
ibd_ev = eigenvals(pcoa.res.ibd)
plot(pcoa.res.ibd$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(ibd_ev[[1]],2))), ylab=paste(c('PCoA2', round(ibd_ev[[2]],2))),
     ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))
legend("topleft", legend=unique(ibd_metadata$Diagnosis), col=unique(colvec), cex=0.8, pch=1)

# Plot the centroids
x = ef.ibd$factors$centroids[,1]*0.3
y = ef.ibd$factors$centroids[,2]*0.3
labels = rownames(ef.ibd$factors$centroids)
text(x=x, y=y, labels=labels)

# Plot the arrows
x = ef.ibd$vectors$arrows[,1]*0.2
y = ef.ibd$vectors$arrows[,2]*0.2
x0 = rep(0, length(x))
y0 = rep(0, length(y))
labels = rownames(ef.ibd$vectors$arrows)
arrows(x1=x, y1=y, x0=x0, y0=y0)
text(x=x, y=y, labels=labels, pos=3)

## Plotting species
# generate the base PCoA figure
colvec <- c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))
sim_ev <- eigenvals(pcoa.res.sim)
plot(pcoa.res.sim$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(sim_ev [[1]],2))), ylab=paste(c('PCoA2', round(sim_ev [[2]],2))))
 
# transpose the normalized counts and get the sample number
Y=t(norm_counts)
n <- nrow(Y)

# standardize the first and second eigenvectors.
ev.stand <- scale(pcoa.res.sim$CA$u[,c(1,2)])

# calculate covariance between taxon abundances and the standardized eigenvectors
S <- cov(Y, ev.stand)

# scale the covariance matrix by the eigenvalues
U <- S %*% diag((pcoa.res.sim$CA$eig[c(1,2)]/(n-1))^(-0.5))

# add as arrows to the plot
x <- U[,1]*0.4
y <- U[,2]*0.4
x0 <- rep(0, length(x))
y0 <- rep(0, length(y))
labels <- seq(8)
arrows(x1=x, y1=y, x0=x0, y0=y0)
text(x=x, y=y, labels=labels, pos=3)

# calculate the length of the vectors
norm_func <- function(x){
  return(sqrt(sum(x^2)))
}

norms=apply(U,1,norm_func)
sorted=sort(norms,index.return=TRUE,decreasing=TRUE)
U.selected <- U[sorted$ix[1:5],]# can be used similarly to plot arrows

##  biplot for the IBD taxa
pcoa.res.ibd = capscale(t(ibd_taxa)~1, distance='bray', na.action='na.omit')
colvec = rainbow(length(unique(ibd_metadata$Diagnosis)))[as.numeric(as.factor(ibd_metadata$Diagnosis))]
eigenvalus = eigenval(pcoa.res.ibd)
plot(pcoa.res.ibd$CA$u[,c(1,2)], col=colvec,
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))),
     ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))
legend("topleft", legend=unique(ibd_metadata$Diagnosis), col=unique(colvec), cex=0.8, pch=1)

# Get the scaled covariance matrix
Y = t(ibd_taxa)
n = nrow(Y)
ev.stand = scale(pcoa.res.ibd$CA$u[,c(1,2)])
S = cov(Y, ev.stand)
U = S %*% diag((pcoa.res.ibd$CA$eig[c(1,2)]/(n-1))^(-0.5))

# Get only the top taxa
norms = apply(U,1,norm_func)
sorted = sort(norms,index.return=TRUE, decreasing=TRUE)
U.selected = U[sorted$ix[1:5],]

# Add vectors for top taxa
x = U.selected[,1]*0.4
y = U.selected[,2]*0.4
x0 = rep(0, length(x))
y0 = rep(0, length(y))
labels = rownames(U.selected)
arrows(x1=x, y1=y, x0=x0, y0=y0)
text(x=x, y=y, labels=labels, pos=3)

## Triplots in ggplot2
# construct the PCoA object, carry out envfit and get taxon covariances
pcoa.res <- capscale(t(norm_counts)~1, distance='bray', na.action='na.omit')
ef <- envfit(pcoa.res,sim_meta,perm=1000, choices=c(1,2), na.rm=TRUE)
Y=t(norm_counts)
n <- nrow(Y)
ev.stand <- scale(pcoa.res$CA$u[,c(1,2)])
S <- cov(Y, ev.stand)
U <- S %*% diag((pcoa.res$CA$eig[c(1,2)]/(n-1))^(-0.5))
U <- U*0.4 # scale U so vectors fit better

# construct the data frames from these objects
sample_df <- data.frame(PCoA1=pcoa.res$CA$u[,1], PCoA2=pcoa.res$CA$u[,2], group=sim_meta)
centroid_df <- data.frame(PCoA1=ef$factors$centroids[,1], PCoA2=ef$factors$centroids[,2],
                          centroid=rownames(ef$factors$centroids))
arrow_df <- data.frame(PCoA1=U[,1], PCoA2=U[,2],
                       taxon=seq(8), x0=rep(0, nrow(U)), y0=rep(0, nrow(U)))


library(ggplot2)

gplot <- ggplot() + geom_point(data=sample_df, aes(x=PCoA1, y=PCoA2, color=group)) +
  geom_text(data=centroid_df, aes(x=PCoA1, y=PCoA2, label=centroid)) +
  geom_segment(data=arrow_df, aes(x=x0, y=y0, xend=PCoA1, yend=PCoA2),arrow=arrow()) +
  geom_text(data=arrow_df, aes(x=PCoA1, y=PCoA2, label=taxon), nudge_x=-0.05, nudge_y=0.05)
gplot
gplot + xlim(-0.6, 1) + theme_minimal() + scale_color_brewer(palette='Dark2') + labs(color='Sample\ngroup')

## build a network
## Step 1: Data preprocessing
ibd_taxa = read.csv("ibd_taxa.csv", row.names=1)

# filter taxa in the IBD data that are present in less than 30 samples
min.prevalence=30
incidence=ibd_taxa
incidence[incidence>0]=1
ibd_taxa_filtered=ibd_taxa[which(rowSums(incidence)>=min.prevalence),]

## Step 2: Computing pairwise associations
# compute the Pearson correlation of Bacteroides thetaiotaomicron and Roseburia hominis 
index.bt=which(rownames(ibd_taxa)=="Bacteroides_thetaiotaomicron")
index.rh=which(rownames(ibd_taxa)=="Roseburia_hominis")
cor(t(ibd_taxa[index.bt,]),t(ibd_taxa[index.rh,]))

## Step 3: Filtering pairwise associations
# How many edges does the Pearson network of the filtered IBD taxon matrix have when only considering correlations below or above 0.6?
ibd_correls=cor(t(ibd_taxa_filtered))
length(ibd_correls[abs(ibd_correls)>0.6])
(length(ibd_correls[abs(ibd_correls)>0.6])-nrow(ibd_correls))/2

## Step 4: Plotting the network
library(igraph)

ibd_correls=cor(t(ibd_taxa_filtered), method="spearman")
ibd_correls[abs(ibd_correls)<0.55]=0

microbe.network=graph_from_adjacency_matrix(ibd_correls,mode="lower",weighted=TRUE, diag=FALSE)

plot(microbe.network)

microbe.network=delete.vertices(microbe.network,degree(microbe.network)==0)

plot(microbe.network, layout=layout.kamada.kawai(microbe.network))

plot(microbe.network, edge.width=E(microbe.network)$weight*10, layout=layout.kamada.kawai(microbe.network))

tkplot(microbe.network, edge.width=E(microbe.network)$weight*10, layout=layout.kamada.kawai(microbe.network))

length(E(microbe.network))

## Step 5: Saving the network
write.graph(graph=microbe.network, file='microbe_network.graphml', format='graphml')

## Interpreting a network
## Mapping information onto nodes
ibd_metadata = read.csv("ibd_metadata.xls", row.names=1)

# compute the average relative abundance of species
# first, we find out the indices of species that made it into the network
node.indices=match(V(microbe.network)$name,rownames(ibd_taxa_filtered))

# then we identify the indices of control and UC/CD samples
control.indices=which(ibd_metadata$Diagnosis=="Control")
ibd.indices=setdiff(c(1:ncol(ibd_taxa_filtered)),control.indices)

# finally, we compute the mean of species in control samples and in UC/CD samples
mean.control=apply(ibd_taxa_filtered[node.indices,control.indices],1,mean)
mean.ibd=apply(ibd_taxa_filtered[node.indices,ibd.indices],1,mean)

# log-ratio of abundance in healthy and IBD samples
logratio=log(mean.control/mean.ibd)

# species that do not change much across conditions will be orange
colors=rep("orange",length(logratio))

# species more abundant in IBD samples will be red
colors[which(logratio < (-1))]="red"

# species more abundant in healthy samples will be green
colors[which(logratio > 1)]="green"

# assign these colors to the nodes in the network
V(microbe.network)$color=colors

# tkplot allows manually moving nodes around
tkplot(microbe.network, edge.width=E(microbe.network)$weight*10, layout=layout.kamada.kawai(microbe.network))

## Clustering the network
clusters=cluster_walktrap(microbe.network)

plot(clusters,microbe.network,col=colors)

table(clusters$membership)

label.colors=rep("black",length(clusters$membership))
label.colors[clusters$membership==1]="darkgreen"
label.colors[clusters$membership==2]="blue"
label.colors[clusters$membership==4]="red"
tkplot(microbe.network, edge.width=E(microbe.network)$weight*10, layout=layout.kamada.kawai(microbe.network), vertex.label.color=label.colors)












## end.