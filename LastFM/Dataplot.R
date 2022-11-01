###Select the corresponding nodes 
ind.select = which(label %in% as.vector(class_select))
Aselect = A[ind.select, ind.select]; Xselect = X[ind.select,]; 
labelselect = label[ind.select];

n = dim(Aselect)[1]; p = dim(Xselect)[2]

# Select the regional popular artists##
dartist = colSums(Xselect);
prop = dartist/dartist_all;
prob = 1 - round(min(n/2, 600)/p, 4);
Xselect = Xselect[, prop > quantile(prop, probs = prob, na.rm = TRUE)]

## Remove the users who like no artists
dnode = rowSums(Xselect)
like = which(dnode > 0)
Xselect = Xselect[like, ]
Aselect = Aselect[like, like]
labelselect = labelselect[like]

n = dim(Aselect)[1]; p = dim(Xselect)[2]

#eigenvectors
K = length(class_select);



Y1 = Xselect%*%t(Xselect); 
eig1 = eigen(Y1)$vectors;
eig1n = eig1[,1:K]/apply(eig1[,1:K], 1, Norm);

M = Aselect%*%Xselect;
Y2 = M%*%t(M)
eig2 = eigen(Y2)$vectors;
eig2n = eig2[,1:K]/apply(eig2[,1:K], 1, Norm);

alpha = 12;
d = rowSums(Aselect);
diagcomp = cbind(1, median(d)/(d+1));
lambda = apply(diagcomp, 1, min)
Y3 = M + alpha*diag(lambda)%*%Xselect;
Y3 = Y3%*%t(Y3);
eig3 = eigen(Y3)$vectors;
eig3n = eig3[,1:K]/apply(eig3[,1:K], 1, Norm);


##Make the plot
K = length(class_select);
pshape = rep(0, length(labelselect))
shapepara = c(1,2,3,4,5,6,8)
for(i in 1:K){
  pshape[labelselect == class_select[i]] = shapepara[i];
}

pdf(file = "Comparison.pdf", width = 9, height = 4)
par(mfrow = c(1,3))

##2d plot
data = eig1n[,1:2];
data[,1] = data[,1]*sign(data[1,1])
data[,2] = data[,2]*sign(data[1,2])
plot(data[,1], data[,2],  col = pshape, pch = pshape, 
              xlab = "First Singular Vector (Normalized)", ylab = "Second Singular Vector (Normalized)", 
     cex = 0.8, 
              cex.lab=0.9, cex.axis = 0.8, mar = c(5,4,1,2)+0.1)
legend("bottomleft", legend = paste("Country", class_select, sep = " "), 
       col = shapepara[1:K], pch = shapepara[1:K], cex = 1)

data = eig2n[,1:2];
data[,1] = data[,1]*sign(data[1,1])
data[,2] = data[,2]*sign(data[1,2])*(-1)
plot(data[,1], data[,2],  col = pshape, pch = pshape, 
     xlab = "First Singular Vector (Normalized)", ylab = "", 
     cex = 0.8, 
     cex.lab=0.9, cex.axis = 0.8, mar = c(5,4,1,2)+0.1)
legend("bottomleft", legend = paste("Country", class_select, sep = " "), 
       col = shapepara[1:K], pch = shapepara[1:K], cex = 1)

data = eig3n[,1:2];
data[,1] = data[,1]/sign(data[1,1])
data[,2] = data[,2]/sign(data[1,2])
plot(data[,1], data[,2],  col = pshape, pch = pshape, 
     xlab = "First Singular Vector (Normalized)", ylab = "", 
     cex = 0.8, 
     cex.lab=0.9, cex.axis = 0.8, mar = c(5,4,1,2)+0.1)
legend(0.15, -0.5, legend = paste("Country", class_select, sep = " "), 
       col = shapepara[1:K], pch = shapepara[1:K], cex = 1)
dev.off()

#Remark: Compare the results
##Compare the methods 
est_eig1 = kmeans(eig1n, K)$cluster
cluster(table(est_eig1, labelselect))$error
NMI(est_eig1, labelselect)

# est_eig2 = kmeans(eig2n, K)$cluster
# cluster(table(est_eig2, labelselect))$error
# NMI(est_eig2, labelselect)

est_eig3 = kmeans(eig3n, K)$cluster
cluster(table(est_eig3, labelselect))$error
NMI(est_eig3, labelselect)

