##################### Overall Review ######################
# Experiment 1. explore what happens when q0 changes with O(1) 
#    Gaussian noise. Sparse + Dense network
# Exp 2a. p = 800
# Exp 2b. p = 1600

# Experiment 2: explore what happens when q0 changes with O(1) 
#    Multinomial noise
# Exp 1a. Relatively sparse case; small p
# Exp 1b. Relatively dense case; large p



##################### Gaussian Covariates ######################
# Experiment 1. Generate mu_X in the way that 20% of them are non-zeros
# Consider loosely connected and densely connected communities
# See how the mis-classification rate will impact the result
# Consider p = 800 in 2a and p = 1600 in 2b

# p = 800
caseno = 3;
n = 1000;
repetition = 50;

K1 = 4; K2 = K1 + 4; prob1seq = seq(0.4, 1, by = 0.05);
errorall = matrix(0, nrow = length(prob1seq), ncol = 7) #store the error rate in each repetition
colnames(errorall) = c("New", "CASC", "ADMM", "SCORE, giant", "nPCA", "oPCA", "giant/all")
errsd = errorall;

for(probjj in 1:length(prob1seq)){
  prob1 = prob1seq[probjj];
  theta = runif(n, min = 0.1, max = 0.4); 
  a = 0.7; 
  source('Simulation.R')
  filename = paste('Simulation Results/Exp2a_prob_', prob1*10, '.Rdata', sep = "");
  save(D0, errormat, P, Pi, Q, Theta, W, caseno, a, n, p, l, K1, K2, mu, prob1seq, repetition, theta, file = filename)
  errorall[probjj,] = colMeans(errormat); errsd[probjj,] = apply(errormat, 2, sd)
  print(paste("prob1 = ", prob1));
}
prop = 0.2; 

save(errorall, errsd, K1, K2, n, p, repetition, Theta, P, Pi, prob1seq, mu, prop,
     file = "Simulation Results/NewExp2a.Rdata")

# p = 2000
caseno = 4;
n = 1000;
repetition = 50;

K1 = 4; K2 = K1 + 4; prob1seq = seq(0.4, 1, by = 0.05);
errorall = matrix(0, nrow = length(prob1seq), ncol = 7) #store the error rate in each repetition
colnames(errorall) = c("New", "CASC", "ADMM", "SCORE, giant", "nPCA", "oPCA", "giant/all")
errsd = errorall;

for(probjj in 1:length(prob1seq)){
  prob1 = prob1seq[probjj];
  theta = runif(n, min = 0.1, max = 0.4); 
  a = 0.7; 
  source('Simulation.R')
  filename = paste('Simulation Results/Exp2b_prob_', prob1*10, '.Rdata', sep = "");
  save(D0, errormat, P, Pi, Q, Theta, W, caseno, a, n, p, l, K1, K2, mu, prob1seq, repetition, theta, file = filename)
  errorall[probjj,] = colMeans(errormat); errsd[probjj,] = apply(errormat, 2, sd)
  print(paste("prob1 = ", prob1));
}

prop = 0.2; 

save(errorall, errsd, K1, K2, n, p, repetition, Theta, P, Pi, prob1seq, mu, prop,
     file = "Simulation Results/NewExp2b.Rdata")

####### Plot Section ############
pdf(file = "Simulation2.pdf", height = 5, width =10)
par(mfrow = c(1,2))
load("Simulation Results/NewExp2a.Rdata")
prob = prob1seq;
colvec = 1:10; pchvec = 1:10;
plot(prob, errorall[ ,1], type = "b", col = colvec[1], pch = pchvec[1],
     xlab = paste("p = ", p, ", Matching Rate ", expression(a), sep = ""), ylab = "Error Rate", ylim = c(0, 0.7), 
     cex.lab = 1.5)
for(i in 2:6){
    lines(prob, errorall[, i], type = "b", col = colvec[i], pch = pchvec[i])
}
  legend("bottomleft", legend = c("CASCORE", "CASC", "SDP", "SCORE", "nPCA", "oPCA"),
         col = 1:6, pch = 1:6, lty = 1, cex = 0.8)
  
  load("Simulation Results/NewExp2b.Rdata")
  prob = prob1seq;
  colvec = 1:10; pchvec = 1:10;
  plot(prob, errorall[ ,1], type = "b", col = colvec[1], pch = pchvec[1],
       xlab = paste("p = ", p, ", Matching Rate ", expression(a), sep = ""), ylab = "Error Rate", ylim = c(0, 0.7), 
       cex.lab = 1.5)
  for(i in 2:6){
    lines(prob, errorall[, i], type = "b", col = colvec[i], pch = pchvec[i])
  }
  legend("bottomleft", legend = c("CASCORE", "CASC", "SDP", "SCORE", "nPCA", "oPCA"),
         col = 1:6, pch = 1:6, lty = 1, cex = 0.8)
dev.off()




##################### O(1) Multinomial Covariates ######################
# Experiment 2. Generate the mu_X for each class by uniform distribution, and
# then normalize to be with sum 1.
# See how the network density will impact the result


# Consider p = 400
caseno = 1;
n = 600;
repetition = 50;
K1 = 4; K2 = K1 + 6; prob1 = 0.7;
thetaseq = seq(0.01, 0.25, by = 0.01);

errorall = matrix(0, nrow = length(thetaseq), ncol = 7) #store the error rate in each repetition
colnames(errorall) = c("New", "CASC", "ADMM", "SCORE, giant", "nPCA", "oPCA", "giant/all")
errsd = errorall;
for(thetajj in 1:length(thetaseq)){
  thetaa = thetaseq[thetajj];
  Nrange = 50; Nmin = 70;
  theta = thetaa + (0.5-thetaa)*(seq(1:n)/n)^2; 
  a = 0.7; 
  source('Simulation.R')
  print(paste("thetaa = ", thetaa));
  filename = paste('Simulation Results/Exp1a_thetaa_', thetaa*10, '.Rdata', sep = "");
  save(D0, errormat, P, Pi, Q, Theta, W, caseno, a, n, p, Nrange, Nmin, l, K1, K2, prob1, repetition, theta, file = filename)
  errorall[thetajj,] = colMeans(errormat); errsd[thetajj,] = apply(errormat, 2, sd)
}


save(errorall, errsd, K1, K2, n, p, repetition, P, Pi, thetaseq, prob1, Nrange, Nmin,
     file = "Simulation Results/NewExp1a.Rdata")

caseno = 2;
n = 600;
repetition = 50;
K1 = 4; K2 = K1 + 6; prob1 = 0.7;
thetaseq = seq(0.01, 0.25, by = 0.01);

errorall = matrix(0, nrow = length(thetaseq), ncol = 7) #store the error rate in each repetition
colnames(errorall) = c("New", "CASC", "ADMM", "SCORE, giant", "nPCA", "oPCA", "giant/all")
errsd = errorall;
for(thetajj in 1:length(thetaseq)){
  thetaa = thetaseq[thetajj];
  Nrange = 50; Nmin = 70;
  theta = thetaa + (0.5-thetaa)*(seq(1:n)/n)^2; 
  a = 0.7; 
  source('Simulation.R')
  print(paste("thetaa = ", thetaa));
  filename = paste('Simulation Results/Exp1b_thetaa_', thetaa*10, '.Rdata', sep = "");
  save(D0, errormat, P, Pi, Q, Theta, W, caseno, a, n, p, Nrange, Nmin, l, K1, K2, prob1, repetition, theta, file = filename)
  errorall[thetajj,] = colMeans(errormat); errsd[thetajj,] = apply(errormat, 2, sd)
}

save(errorall, errsd, K1, K2, n, p, repetition, P, Pi, thetaseq, prob1, Nrange, Nmin,
     file = "Simulation Results/NewExp1b.Rdata")


####### Plot Section ############
pdf(file = "Simulation1.pdf", height = 5, width =10)
par(mfrow = c(1,2))
colvec = 1:10; pchvec = 1:10;

load("Simulation Results/NewExp1a.Rdata")
thetaa = thetaseq;
plot(thetaa, errorall[ ,1], type = "b", col = colvec[1], pch = pchvec[1], cex = 0.8,
     xlab = expression("Minimal"~"Theta "*theta), ylab = "Error Rate", 
     ylim = c(0, 0.6), cex.lab = 1.5)
for(i in 2:6){
  lines(thetaa, errorall[ , i], type = "b", col = colvec[i], pch = pchvec[i], cex = 0.8)
}
for(i in 1:6){
  lines(thetaa, errorall[ , i], col = colvec[i])
}
legend("topright", legend = c("CASCORE", "CASC", "SDP", "SCORE", "nPCA", "oPCA"),
       col = 1:6, pch = 1:6, lty = 1, cex = 0.8)


load("Simulation Results/NewExp1b.Rdata")
thetaa = thetaseq;
plot(thetaa, errorall[ ,1], type = "b", col = colvec[1], pch = pchvec[1], cex = 0.8,
     xlab = expression("Minimal"~"Theta "*theta), ylab = "Error Rate", ylim = c(0, 0.6), 
     cex.lab = 1.5)
for(i in 2:6){
  lines(thetaa, errorall[ , i], type = "b", col = colvec[i], pch = pchvec[i], cex = 0.8)
}
for(i in 1:6){
  lines(thetaa, errorall[ , i], col = colvec[i])
}
legend("topright", legend = c("CASCORE", "CASC", "SDP", "SCORE", "nPCA", "oPCA"),
       col = 1:6, pch = 1:6, lty = 1, cex = 0.8)
dev.off()


