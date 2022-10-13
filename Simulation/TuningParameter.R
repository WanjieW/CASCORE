#setwd("~/Community Detection")
#setwd("~/Project/Proj-Social Network/Final Code")
library(Matrix)
library(igraph)
library(pracma)
source("Mtable.R")

##Adjust the original function to record the results w.r.t. different alpha ###

CASCORE = function(Adj, Covariate, K, alpha = NULL, alphan = 5, itermax = 100, startn = 10){
  # Inputs:
  # 1) Adj: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) Covariate: an n by p covariate matrix
  # 3) K: a positive integer which is no larger than n. This is the predefined number of communities.
  # 4) alpha: a vector of positive numbers to tune the weight of covariate matrix
  # 5) alphan: if alpha is not given, alphan is required to find a proper weight by grid search. The 
  #            default value is 5.
  
  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed. Default value 100.
  # 2) nstart: R will try startn different random starting assignments and 
  #            then select the one with the lowest within cluster variation. Default value 10.
  
  # Outputs:
  # 1) alpha: a list of alpha selected
  # 2) est: A list of labels corresponding to all choices of alpha
  # 3) prop: A list of the within-cluster sum of variances corresponding to alpha
  
  if(!isSymmetric(Adj)) stop("Error! Adj is not symmetric!")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")
  if(dim(Adj)[1] != dim(Covariate)[1]) stop("Error! Incompatible!")
  #  if(alpha < 0) stop("Negative Alpha")
  
  #Regularity check
  estall = rep(NA, dim(Adj)[1]);
  netrow = rowSums(Adj);
  covrow = rowSums(abs(Covariate));
  ind_reg = which(netrow != 0 | covrow != 0)
  Adj = Adj[ind_reg, ind_reg];
  Covariate = Covariate[ind_reg, ];
  
  ##Algorithm
  n = dim(Adj)[1]
  d = rowSums(Adj);
  X = Adj %*% Covariate
  
  diagcomp = cbind(1, median(d)/(d+1));
  lambda = apply(diagcomp, 1, min)
  
  if(is.null(alpha)){
    d.net = sort(abs(eig(Adj)), decreasing = TRUE);
    alphaupper = d.net[1]*log(n)/mean(d);
    alphalower = max(0.05, d.net[K]/4);
    
    alpha = seq(alphalower, alphaupper, length.out = alphan);
    #    print(alpha)
    est1 = matrix(0, alphan, n); est2 = est1;
    prop1 = rep(0, alphan); prop2 = rep(0, alphan)
    
    for(i in 1:alphan){
      Newmat = X + alpha[i]*diag(lambda)%*%Covariate;
      zz = Newmat%*%t(Newmat)
      c = eigen(zz)
      vec = c$vectors
      vecn = vec[,1:K]/apply(vec[,1:K], 1, Norm);
      result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
      if (result$ifault==4) { result = kmeans(X, K,  iter.max = itermax, nstart = startn, algorithm="Lloyd"); }
      prop2[i] = result$tot.withinss;
      est2[i,] = as.factor(result$cluster);
    }
    est = est2[which.min(prop2), ];
    #    print(prop2)
  }
  else{
    alphan = length(alpha);
    est1 = matrix(0, alphan, n); est2 = est1;
    prop1 = rep(0, alphan); prop2 = rep(0, alphan)
    
    for(i in 1:alphan){
      Newmat = X + alpha[i]*diag(lambda)%*%Covariate;
      zz = Newmat%*%t(Newmat)
      c = eigen(zz)
      vec = c$vectors
      vecn = vec[,1:K]/apply(vec[,1:K], 1, Norm);
      result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
      if (result$ifault==4) { result = kmeans(X, K,  iter.max = itermax, nstart = startn, algorithm="Lloyd"); }
      prop2[i] = result$tot.withinss;
      est2[i,] = as.factor(result$cluster);
    }
    est = est2[which.min(prop2), ];
  }
  estall[ind_reg] = est;
  return(list(alpha = alpha, est = est2, prop = prop2))
}



############ Set the Network  ############
# Simulate the Network
n = 600;
K1 = 4; 

theta = 0.02 + (0.5-0.02)*(seq(1:n)/n)^2; 
Theta = diag(theta); # node degree heterogeneity
P = diag(rep(0.6, K1)) + 0.4*ones(K1);

# Set the labels
set.seed(2019)
l = sample(1:K1, n, replace=TRUE); # node labels
Pi = matrix(0, n, K1) # label matrix
for (k in 1:K1){Pi[l == k, k] = 1}

# Get the expected adjacency matrix
Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;

#Generate the adjacency matrix
set.seed(2022)
A = matrix(runif(n*n, 0, 1), nrow = n);
A = Omega - A;
A = 1*(A >= 0)
diag(A) = 0
A <- as.matrix(forceSymmetric(A))

############ Set the Covariates  ############
p = 800; # length of the covariate vector
K2 = K1 + 6; 

set.seed(2022)
Q = matrix(runif(p*K2, 0, 1), nrow = p, ncol = K2)
Q = sweep(Q,2,colSums(Q),`/`) 

set.seed(2022)
prob1 = 0.8;
W = matrix(0, nrow = n, ncol = K2); 
for(jj in 1:n) {
  pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
  if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
  else 
    W[jj, sample(K2, 1, prob = pp)] = 1;
}
W = t(W);

D0 = Q %*% W # expectation of covariate matrix 

## Generate the random matrix with multinomial covariates###
set.seed(2022)
Nrange = 50; Nmin = 70;
D = matrix(0, n, p);
for (i in 1: ncol(D0)){
  N = round(runif(n)*Nrange+Nmin);
  D[i, ] = rmultinom(1, N[i], D0[, i])
}

alphamax = 5*eig(A)[1]; alphamin = 0.05;
alpha = seq(alphamin, alphamax, length = 40);
print(alpha)

result = CASCORE(A, D, K1, alpha = alpha)
errormat = alpha*0;
for(ii in 1:length(alpha)){
  est_cascore = result$est[ii,]
  comp = table(est_cascore, l);
  err_cascore = cluster(comp)$error;
  errormat[ii] = err_cascore;
  print(ii)
}

result2 = CASCORE(A, D, K1)
errormat2 = rep(0, 5);
for(ii in 1:5){
  est_cascore = result2$est[ii,]
  comp = table(est_cascore, l);
  err_cascore = cluster(comp)$error;
  errormat2[ii] = err_cascore;
  print(ii)
}

est_cascore = result2$est[which.min(result2$prop), ];
comp = table(est_cascore, l);
err_cascore = cluster(comp)$error;
erropt = err_cascore;


pdf(file = "Tuning.pdf", height = 4, width =8)
plot(alpha, errormat, type = "l", xlab = "Tuning Parameter alpha", 
     ylab = "Error Rate")
abline(v = result2$alpha[which.min(result2$prop)],  col = 2, lty = 2)

dev.off()
