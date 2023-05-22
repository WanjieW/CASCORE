#setwd("~/Project/Proj-Social Network/New")
library(Matrix)

A = as.matrix(read.table("connection.txt", header=F))
# X = as.matrix(read.table("Nfreq.txt", header=T))
load("Nfreq.Rdata")
colnames(A) = NULL
n = dim(A)[1]; p = dim(X)[2]

degree = colSums(A)

newind = which(degree >= 50)
length(newind) #326

require(igraph)
graph1 = graph.adjacency(A[newind, newind], mode = "max")
summary(graph1)
set.seed(2015)
L = layout_nicely(graph1, dim = 2)


# choose the number of communities
#Generate a figure about the eigenvalues of A
Aeig = eig(A)
plot(Aeig[1:6])
K = 3 #we used to try K  from 3 to 6

source("CASCORE.R")
d = colSums(A); alpha = mean(d)/2;
est_cascore = matrix(0, n, nrow = 4, ncol = n); 
for(K in 3:6){  
est_cascore[K-2,] = CASCORE(A, X, K, alpha = alpha, startn = 50)
}

#Match the labels for different K
est_cascore_large = est_cascore[,newind]; 
for(kk in 4:6){
  est_cascore_large = est_cascore[,newind]; 
  tmp = table(est_cascore_large[kk-3,], est_cascore_large[kk-2, ])
  tmplabel = est_cascore[kk - 2,]; 
  recordj = NULL;
  for(ii in 1:nrow(tmp)){
    j = which.max(tmp[ii,]);  
    est_cascore[kk - 2, tmplabel == j] = ii; 
    recordj = c(recordj, j); 
  }
  j = setdiff(1:kk, recordj); 
  est_cascore[kk-2, tmplabel == j] = kk;
}

# Plot 4 figures with K = 3:6 (Figure in Supplementary Material)
#pdf(file = "CitationK.pdf", height = 8, width =8)
setEPS()  
postscript("CitationK.eps", height = 8, width =8)  
par(mfrow = c(2,2), mai = c(0.1, 0.1, 0.2, 0.1))
for(kk in 3:6){
plot(graph1, layout = L, vertex.color = est_cascore[kk - 2, newind], vertex.label = NA, 
     vertex.size = 5, main = paste("CASCORE, K = ", kk))
}
dev.off()


#Explore the case K = 5
#Find the hot corpus
K = 5;
for(kk in 1:5){
  subclass = which(est_cascore[K - 2, ] == kk);
  print(names(sort(apply(X[subclass,], 2, sum), decreasing = TRUE)[1:10]))
}

#Plot
#pdf(file = "Citation.pdf", height = 5, width =7)
setEPS()  
postscript("Citation.eps", height = 3, width =6)  
par(mfrow = c(1,1), mai = c(0.1, 0.1, 0.2, 0.1))
kk = 5; 
csize = summary(as.factor(est_cascore[kk-2,]));
newlabel = rank(csize)
plot(graph1, layout = L, xlim = c(-1, 2), vertex.color = newlabel[est_cascore[kk - 2, newind]], vertex.label = NA, 
     vertex.size = 5)
legendnames = c("Variable Selection \n (Regression)\n", "Large-Scale \n Multiple Testing", "Biostatistics", 
                "Bayesian\n", "Variable Selection \n (Semiparametric)\n");
legend(1.2, 0.8, legend = legendnames, col = categorical_pal(kk), pch = 16, 
       cex = 1, y.intersp = 0.5)
dev.off()
### End of Plot ####

#find the most popular papers in each community
paperList = read.table("paperList.txt", sep=",", stringsAsFactors=F, header=T, 
                       na.strings = "")
for(ii in 1:kk){
  dcomm = degree*0;
  dcomm[est_cascore[kk-2, ] == ii] = degree[est_cascore[kk-2, ] == ii]; 
  indcomm = order(dcomm, decreasing = TRUE);
  print(paperList[indcomm[1:10], 3])
}

#number of nodes in each community
allsize = summary(as.factor(newlabel[est_cascore[kk-2,]]))
print(allsize)

ix = components(graph.adjacency(A))
componentLabel = ix$membership
GiantLabel = which(componentLabel == which.max(ix$csize))
length(GiantLabel) #2179
Giant = A[GiantLabel, GiantLabel]
giantsize = summary(as.factor(newlabel[est_cascore[kk-2, GiantLabel]]))
print(giantsize)

print(allsize - giantsize)
## Consider an example isolated paper
nodes = c(2893, 2481);
degree[nodes]

paperList[2893, 3]
#Bayesian pseudo-empirical-likelihood intervals for complex surveys
legendnames[newlabel[est_cascore[kk-2, 2893]]]
#Bayesian

paperList[2481, 3]
#Testing dependence among serially correlated multicategory variables
legendnames[newlabel[est_cascore[kk-2, 2481]]]
#Large-Scale Multiple Testing

write.table(newlabel[est_cascore[kk-2, ]], file = "Result.txt", row.names = paperList[,3])
save.image(file = "Citation.Rdata");
