### run algorithm on the countries and compare ###
allcomp <- function(A, X, label, class){
  A1 = A; X1 = X;
  # ix = components(graph.adjacency(A1))
  # componentLabel = ix$membership
  # 
  # giantlabel = label[which(componentLabel == which.max(ix$csize))]
  # Giant = A1[which(componentLabel == which.max(ix$csize)), 
  #            which(componentLabel == which.max(ix$csize))]
  # giantp = length(giantlabel) / (dim(A1)[1]) 
  
  result = rep(0,4)
  names(result) = c("CASCORE", 
                    "CASC",
                    "ADMM",
                    "nPCA")
  timeresult = result;
  NMIresult = result;
  error = result; 
  
  kk = length(unique(label))
  
  # CASCORE
  t1 = Sys.time()
  est_cascore = CASCORE(A1, X1, kk)
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
   
  # ADMM
  t1 = Sys.time()
  est_admm = admm(A1, X1, 0.5, kk, 2, 0.5, 100, 5)
  timeresult["ADMM"] = Sys.time() - t1
  est_admm[is.na(est_admm)] = which.max(summary(as.factor(est_admm)))
  result["ADMM"] = NMI(est_admm, label)
  if(kk <= 7) error["ADMM"] = cluster(table(est_admm, label))$error;


  # # SCORE
  # t1 = Sys.time()
  # kkk = length(unique(giantlabel))
  # est_score = SCORE(Giant, kkk)
  # est_score[is.na(est_score)] = which.max(summary(as.factor(est_score)))
  # timeresult["SCORE"] = Sys.time() - t1
  # comp = table(est_score, giantlabel)
  # result["SCORE"] = NMI(est_score, giantlabel)
  # 
  # if(kkk <= 7) {
  #   error["SCORE"] = cluster(comp)$error;
  # }
  # 
  # # oPCA
  # t1 = Sys.time()
  # est_opca = oPCA(A1, kk)
  # timeresult["oPCA"] = Sys.time() - t1
  # result["oPCA"] = NMI(est_opca, label)
  # if(kk <= 7) error["oPCA"] = cluster(table(est_opca, label))$error;

  # nPCA
  t1 = Sys.time()
  est_npca = nPCA(A1, kk)
  timeresult["nPCA"] = Sys.time() - t1
  result["nPCA"] = NMI(est_npca, label)
  if(kk <= 7) error["nPCA"] = cluster(table(est_npca, label))$error;

  
  return(
    result = list(
      class = class,
      error = error, 
      NMI = result,
      time = timeresult,
      Sizes = dim(A1)[1],
#      GiantProp = giantp,  
      est = list(cascore = est_cascore, casc = est_casc, 
                 admm =est_admm, opca = est_opca), 
      label = label
    )
  )
}

###Select the corresponding nodes 
ind.select = which(label %in% as.vector(class_select))
Aselect = A[ind.select, ind.select]; Xselect = X[ind.select,]; 
labelselect = label[ind.select];

n = dim(Aselect)[1]; p = dim(Xselect)[2]

# Select the regional popular artists##
dartist = colSums(Xselect);
prop = dartist/dartist_all;
prob = 1 - round(min(n/2, 800)/p, 4);
Xselect = Xselect[, prop > quantile(prop, probs = prob, na.rm = TRUE)]

## Remove the users who like no artists
dnode = rowSums(Xselect)
like = which(dnode > 0)
Xselect = Xselect[like, ]
Aselect = Aselect[like, like]
labelselect = labelselect[like]

n = dim(Aselect)[1]; p = dim(Xselect)[2]

#run results
result <- allcomp(Aselect, Xselect, labelselect, class_select)
