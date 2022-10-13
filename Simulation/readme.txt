Various Functions.r includes the code of following functions:

1) With node attributes:
   CASCORE 
   CAclustering
   SDP

   Note: CASCORE and CAclustering require an n by p covariate matrix.

2) With only adjacency matrix:
   SCORE
   nPCA
   oPCA
   Note: SCORE and nPCA can only be applied on connected
         graphs.  

CASCORE, SDP, SCORE, nPCA and oPCA return a factor indicating nodes' labels. 
