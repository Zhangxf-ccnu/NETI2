#test the performance of NETI2 using simulated data
rm(list=ls())

library("NETI2")
source("evaluation.metric.R")
#source("generate.data.R")

p = 100
n = 50
K = 4
network.type = "ER"
umin = 0.5
umax = 1

rtimes = 1

lambda_list = exp(seq(log(2), log(0.05),length=10))
tau_list = c(0.3,0.5,0.7)
delta = 0.5

TP.z = matrix(0, length(lambda_list), length(tau_list))
FP.z = matrix(0, length(lambda_list), length(tau_list))
TP.y = matrix(0, length(lambda_list), length(tau_list))
FP.y = matrix(0, length(lambda_list), length(tau_list))

for (r in 1:rtimes){
  cat("r=",r)
  # generare simulated data
  dat = generate.data(p, n, K, network.type, umin, umax )

  for (i_lam in 1:length(lambda_list)){
    cat("i_lam=",i_lam)
    for (i_tau in 1:length(tau_list)){
      # run NETI2
      result = NETI2(dat$X,dat$purity,lambda_list[i_lam],tau_list[i_tau],delta)

      # compute TF Edges and FP Edges of non-cancerous network
      per.z = evaluation.metric(dat$theta.z, result$theta.z)
      TP.z[i_lam, i_tau] = TP.z[i_lam, i_tau] + per.z$TP
      FP.z[i_lam, i_tau] = FP.z[i_lam, i_tau] + per.z$FP

      # compute TF Edges and FP Edges of cancerous networks
      TP_tumor=0; TN_tumor=0; FP_tumor=0; FN_tumor=0;
      for (i_k in 1 : K){
        per.y = evaluation.metric(dat$theta.y[[i_k]], result$theta.y[[i_k]])
        TP_tumor = TP_tumor + per.y$TP
        FP_tumor = FP_tumor + per.y$FP
      }

      TP.y[i_lam, i_tau] = TP.y[i_lam, i_tau] + TP_tumor
      FP.y[i_lam, i_tau] = FP.y[i_lam, i_tau] + FP_tumor

    }
  }
}

# compute the average of precision and recall
TP.z = TP.z/rtimes
FP.z = FP.z/rtimes
TP.y = TP.y/rtimes
FP.y = FP.y/rtimes

# plot the number of true positive edges against the number of false positive edges
# The Performance of Non-cancerous network
type_list = c("o","o","o")
pch_list = c(0,1, 2)
plot(FP.z[,1], TP.z[,1], col="blue", type = type_list[1],
     pch=pch_list[1], lwd=2, xlim = c(0,600),  ylim = c(0,80),  xlab = "FP Edges", ylab = "TP Edges")
title(main = "Performance of NETI2 on ER network of non-cancerous network")

for (i in 1:3){
  points(FP.z[,i], TP.z[,i],
         col="blue",type = type_list[i], pch=pch_list[i], lwd=2)
}

method = c(expression(paste("NETI2 ", tau == 0.3)), expression(paste("NETI2 ", tau == 0.5)),
           expression(paste("NETI2 ", tau == 0.7)) )
pch=  c(0,1, 2)
col =  rep("blue",3)
legend(400,40, legend=  method, ncol=1, col = col, pch=pch, lwd=2,bty = "n")

###############
# The Performance of  cancerous networks
type_list = c("o","o","o")
pch_list = c(0,1, 2)
plot(FP.y[,1], TP.y[,1], col="blue", type = type_list[1],
     pch=pch_list[1], lwd=2, xlim = c(0,2000),  ylim = c(0,400),  xlab = "FP Edges", ylab = "TP Edges")
title(main = "Performance of NETI2 on ER network  of cancerous network")

for (i in 1:3){
  points(FP.y[,i], TP.y[,i],
         col="blue",type = type_list[i], pch=pch_list[i], lwd=2)
}

method = c(expression(paste("NETI2 ", tau == 0.3)), expression(paste("NETI2 ", tau == 0.5)),
           expression(paste("NETI2 ", tau == 0.7)) )
pch=  c(0,1, 2)
col =  rep("blue",3)

legend(1200,150, legend=  method, ncol=1, col = col, pch=pch, lwd=2,bty = "n")

