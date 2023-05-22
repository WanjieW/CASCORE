library(pracma)
library(Matrix)


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
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.
  
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
  
  #  lambda = max(log(n), quantile(d, probs = 0.25))/(d + max(log(n), median(d, probs = 0.25)));
  lambda = log(n)/(d + log(n));
  
  if(is.null(alpha)){
    d.net = sort(abs(eig(Adj)), decreasing = TRUE);
    alphaupper = d.net[1]*log(n)/mean(d);
    alphalower = max(0.05, d.net[K]/4);
    
    alpha = seq(alphalower, alphaupper, length.out = alphan);
    #    print(alpha)
    est1 = matrix(0, alphan, n); est2 = est1;
    prop1 = rep(0, alphan); prop2 = rep(0, alphan)
    
    for(i in 1:alphan){
      #      Newmat = X + diag(alpha[i]/(d+1))%*%Covariate;
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
    print(alpha[which.min(prop2)])
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
  return(estall)
}



