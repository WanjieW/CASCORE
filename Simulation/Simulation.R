#setwd("~/Community Detection")
#setwd("~/Project/Proj-Social Network/Final Code")
library(Matrix)
library(igraph)
library(pracma)
source("Various Functions.R")
source("Mtable.R")

############ Parameter Set  ############

# Simulate the Network
#n = 1000; # number of nodes
#K1 = 4; # number of communities revealed by adjacency matrix

# use theta from "ExpSets.R"
Theta = diag(theta); # node degree heterogeneity

# Community by Community probability matrix
# use a from "ExpSets.R"
P = diag(rep(a, K1)) + (1-a)*ones(K1);
if(caseno >= 3){
  sparse = diag(c(1/2, 1/2, 1, 1)); 
  P = sparse%*%P%*%sparse;
}

# Set the labels
set.seed(2019)
l = sample(1:K1, n, replace=TRUE); # node labels

Pi = matrix(0, n, K1) # label matrix
for (k in 1:K1){
  Pi[l == k, k] = 1
}

# Get the expected adjacency matrix
Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;

############ End Parameter Set  ############
errormat = matrix(0, nrow = repetition, ncol = 7) #store the error rate in each repitition
colnames(errormat) = c("New", "CASC", "ADMM", "SCORE, giant", "nPCA", "oPCA", "giant/all")

#Calculation of the error 
for (ii in 1: repetition){
  set.seed(ii+9999)
  
  if(caseno == 1){
    # For this case, we consider low dimensional multivariate covariates
    p = 400; # length of the covariate vector
    
    Q = matrix(runif(p*K2, 0, 1), nrow = p);
    Q = sweep(Q,2,colSums(Q),`/`) 
    
    W = matrix(0, nrow = n, ncol = K2); 
    
    # use prob1 from "ExpSets.R"
    for(jj in 1:n) {
      pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
      if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
      else 
        W[jj, sample(K2, 1, prob = pp)] = 1;
    }
    W = t(W)
    D0 = Q %*% W # expectation of covariate matrix 
  }
  
  
  if(caseno == 2){
    # For this case, we consider high dimensional multinomial covariates 
    p = 800; # length of the covariate vector
    
    Q = matrix(runif(p*K2, 0, 1), nrow = p);
    Q = sweep(Q,2,colSums(Q),`/`) 
    
    W = matrix(0, nrow = n, ncol = K2); 
    
    # use prob1 from "ExpSets.R"
    for(jj in 1:n) {
      pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
      if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
      else 
        W[jj, sample(K2, 1, prob = pp)] = 1;
    }
    W = t(W)
    D0 = Q %*% W # expectation of covariate matrix 
  }
  
  if(caseno == 3){
    # For this case, we consider high dimensional covariates and more topics than the communities
    p = 800; # length of the covariate vector
    
    mu = 0.2;
    Q = matrix(rbinom(p*K2, 1, 0.2)*mu, nrow = p);
    W = matrix(0, nrow = n, ncol = K2); 
    
    # use prob1 from "ExpSets.R"
    
    for(jj in 1:n) {
      pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
      if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
      else 
        W[jj, sample(K2, 1, prob = pp)] = 1;
    }
    W = t(W)
    
    D0 = Q %*% W # expectation of covariate matrix 
  }
  
  
  if(caseno == 4){
    # For this case, we consider high dimensional covariates and more topics than the communities
    p = 2000; # length of the covariate vector
    
    mu = 0.2;
    Q = matrix(rbinom(p*K2, 1, 0.2)*mu, nrow = p);
    W = matrix(0, nrow = n, ncol = K2); 
    
    # use prob1 from "ExpSets.R"
    
    for(jj in 1:n) {
      pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
      if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
      else 
        W[jj, sample(K2, 1, prob = pp)] = 1;
    }
    W = t(W)
    
    D0 = Q %*% W # expectation of covariate matrix 
  }
  
  D = matrix(0, n, p)
#  N = switch(caseno, rep(100, n), rep(100, n), round(runif(n)*Nrange+Nmin), round(runif(n)*Nrange+Nmin))
  #the parameter in multivariate distribution 
  for (i in 1: ncol(D0)){
#    D[i, ] = rmultinom(1, N[i], D0[, i])
    if(caseno < 3) {
      N = round(runif(n)*Nrange+Nmin);
      D[i, ] = rmultinom(1, N[i], D0[, i])
      }
    else{# consider normal noise
      D[i,] = rnorm(p, mean = D0[,i], sd = 1);
    }
  }

  A = matrix(runif(n*n, 0, 1), nrow = n);
  A = Omega - A;
  A = 1*(A >= 0)
  diag(A) = 0
  A <- as.matrix(forceSymmetric(A))
  is.igraph(A) # [1] FALSE
  ix = components(graph.adjacency(A))
  componentLabel = ix$membership
  giantLabel = which(componentLabel == which.max(ix$csize))
  Giant = A[giantLabel, giantLabel]
  errormat[ii, "giant/all"] = length(giantLabel)/n;


  # New approach: CA-SCORE
  est_cascore = CASCORE(A, D, K1)
  comp = table(est_cascore, l);
  err_cascore = cluster(comp)$error;
  errormat[ii, "New"] = err_cascore;

  # CA-clustering
  est_casc = CAclustering(A, D, K1)
  comp = table(est_casc, l);
  err_casc = cluster(comp)$error;
  errormat[ii, "CASC"] = err_casc;

  # SCORE: giant component
  est_score = SCORE(Giant, K1)
  comp = table(est_score, l[giantLabel]);
  errormat[ii, "SCORE, giant"] = cluster(comp)$error;

  # nPCA
  est_npca = nPCA(A, K1)
  comp = table(est_npca, l);
  errormat[ii, "nPCA"] = cluster(comp)$error;

  # oPCA
  est_opca = oPCA(A, K1)
  comp = table(est_opca, l);
  errormat[ii, "oPCA"] = cluster(comp)$error;

  # admm
  est_admm = admm(A, C = D, K = K1, lambda = 0.5, alpha = 2, rho = 0.5, TT = 100, tol = 100)
  comp = table(est_admm, l);
  errormat[ii, "ADMM"] = cluster(comp)$error;

  print(ii)
#  print(errormat[ii,])
}
print(colMeans(errormat))
