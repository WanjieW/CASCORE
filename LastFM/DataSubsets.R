### run algorithm on the countries and compare ###
allcomp <- function(A, X, label, class){
  A1 = A; X1 = X;
  
  result = rep(0,5)
  names(result) = c("CASCORE", 
                    "CASC",
                    "SDP",
                    "Net-Based", "Cov-Based")
  timeresult = result;
  NMIresult = result;
  error = result; 
  
  kk = length(unique(label))
  
  # CASCORE
  t1 = Sys.time()
  d = rowSums(A1);
  est_cascore = CASCORE(A1, X1, kk, alpha = mean(d)/2, startn = 50)
  est_cascore[is.na(est_cascore)] = which.max(summary(as.factor(est_cascore)))
  timeresult["CASCORE"] = Sys.time() - t1
  result["CASCORE"] = NMI(est_cascore, label)
  if(kk <= 7) error["CASCORE"] = cluster(table(est_cascore, label))$error; 
  

    # CAclustering
  t1 = Sys.time()
  est_casc = CAclustering(A1, X1, kk)
  timeresult["CASC"] = Sys.time() - t1
  est_casc[is.na(est_casc)] = which.max(summary(as.factor(est_casc)))
  result["CASC"] = NMI(est_casc, label)
  if(kk <= 7) error["CASC"] = cluster(table(est_casc, label))$error;
   
  # SDP
  t1 = Sys.time()
  est_sdp = admm(A1, X1, 0.5, kk, 2, 0.5, 100, 5)
  timeresult["SDP"] = Sys.time() - t1
  est_sdp[is.na(est_sdp)] = which.max(summary(as.factor(est_sdp)))
  result["SDP"] = NMI(est_sdp, label)
  if(kk <= 7) error["SDP"] = cluster(table(est_sdp, label))$error;

  
  # Net-Based (RSC)
  t1 = Sys.time()
  est_rsc = Net_based(A1, kk, startn = 50)
  timeresult["Net-Based"] = Sys.time() - t1
  result["Net-Based"] = NMI(est_rsc, label)
  if(kk <= 7) error["Net-Based"] = cluster(table(est_rsc, label))$error;
  
  # Cov-Based (SpectralGem)
  t1 = Sys.time()
  est_spec = Cov_based(X1, kk, startn = 50)
  timeresult["Cov-Based"] = Sys.time() - t1
  result["Cov-Based"] = NMI(est_spec, label)
  if(kk <= 7) error["Cov-Based"] = cluster(table(est_spec, label))$error;
  
  return(
    result = list(
      class = class,
      error = error, 
      NMI = result,
      time = timeresult,
      Sizes = dim(A1)[1],
      est = list(cascore = est_cascore, casc = est_casc, 
                 sdp =est_sdp, net = est_rsc, cov = est_spec),  
      label = label
    )
  )
}

###Select the corresponding nodes 
ind.select = which(label %in% as.vector(class_select))
Aselect = A[ind.select, ind.select]; Xselect = X[ind.select,]; 
labelselect = label[ind.select];

n = dim(Aselect)[1]; p = dim(Xselect)[2]
kk = length(unique(labelselect))

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
dselect = rowSums(Aselect)

n = dim(Aselect)[1]; p = dim(Xselect)[2]

#run results
result <- allcomp(Aselect, Xselect, labelselect, class_select)
