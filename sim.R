rm(list=ls())
library(pacman)
p_load(igraph, caret,SparseDC,VGAM,Matrix,pROC,kernlab,foreach,doParallel,Metrics) 
pack <- c('igraph','SparseDC','VGAM','caret','Matrix','Metrics')
# registerDoParallel(2)
setwd('D:/民大/项目/Directed network clustering/results/example of ABLR')
source("utils.R")

# nodes number
n <- 50
# dimension
v <- 2 * n
# alpha in paper
s <- matrix(runif(v,-0.1,0.1),v,1)
# beta in paper
b <- matrix(runif(v,-0.1,0.1),v,1)
# intercept b in paper
a0 <- runif(1,-0.1,0.1)
# link prob setting
a1<-0.2; b1<-sqrt(0.2) # normal parameters, group 1 b = 0.5
a2<-0.1; b2<-sqrt(0.2)
a3<-0.6; b3<-sqrt(0.2)
a4<-0.3; b4<-sqrt(0.2)
# lambda setting
lambda_0 <- 0.01
lambda_1 <- 0.01
# set the limitation of update time
count_limit <- 1 
# change 
delta <- 0.05
# sample size
n_sample_1 <- n_sample_2 <- 100

# mark the active sending and receiving nodes
node_in <- node_out <- rep(0,n)
node_out[1:(ceiling(n/4) + n/2)] <- node_in[1:(ceiling(n/4) + n/2)] <-1

# label of sample
label <- c(rep(1,n_sample_1),rep(0,n_sample_2))

## generate data
# generate adjacency matrix list
Gra <- model.weighted.sparse(rep = n_sample_2,a1,b1,a2 ,b2,a3 ,b3,a4 ,b4,delta,n)
# sample in group 1
A1 <- Gra$Gra1
# sample in group 2
A2 <- Gra$Gra2
# combine sample
A <- c(A1,A2)
# split data
parts = createDataPartition(label, p = .8, list = F)

## classification
## FSVM
# using feature in-strength
res_in <- Hesam_SVM(A1,A2,label,feature = 'inward',parts = parts)
# using feature out-strength
res_out <- Hesam_SVM(A1,A2,label,feature = 'outward',parts = parts)

# calculate AUC of in- and out-nodes
node_in_NSVM <- rep(0,n)
node_in_NSVM[res_in[["nodes_select"]]] <- 1
node_out_NSVM <- rep(0,n)
node_out_NSVM[res_out[["nodes_select"]]] <- 1
# AUC result
auc_in_NSVM <- auc(node_in, node_in_NSVM)
auc_out_NSVM <- auc(node_out, node_out_NSVM)

## ABLR
# construct A_tuta
A_tuta <- lapply(A, function(x){as.matrix(bdiag(x,t(x)))})
# training data
A_tuta_train <- A_tuta[parts]
# test data
A_tuta_test <- A_tuta[-parts]
# label
label_A_tuta_train <- label[parts]
label_A_tuta_test <- label[-parts]
# classify
res_ABLR <- Dir_LOGI_BCD_new(A_tuta_train,A_tuta_test,label_A_tuta_train,label_A_tuta_test,
                             s,b,a0,lambda_0 = lambda_0,lambda_1=lambda_1,tol = 1e-7,
                             eta = 0.5, count = 1, count_limit = count_limit)
# calculate AUC of in- and out-nodes
ABLR_s <- unlist(res_ABLR$s)
s_send <- which(ABLR_s[1:n] != 0)
s_recv <- which(ABLR_s[(n+1):(2*n)] != 0)

node_send_ABLR <- node_recv_ABLR <- rep(0,n)
node_send_ABLR[s_send] <- node_recv_ABLR[s_recv] <- 1
# AUC result
auc_send_ABLR <- auc(node_out, node_send_ABLR)
auc_recv_ABLR <- auc(node_in, node_recv_ABLR)

accuracy_in <- res_in$acc    # accuracy
AUC_in <- res_in$auc         # AUC

accuracy_out <- res_out$acc    # accuracy
AUC_out <- res_out$auc         # AUC

accuracy_ABLR <- res_ABLR$acc
AUC_ABLR <- res_ABLR$auc

results <- c(accuracy_in, accuracy_out, accuracy_ABLR,
             AUC_in, AUC_out, AUC_ABLR,
             auc_in_NSVM, auc_out_NSVM,
             auc_recv_ABLR,auc_send_ABLR)

names(results) <- c('acc_in','acc_out','acc_ABLR',
                       'AUC_in', 'AUC_out', 'AUC_ABLR',
                       'auc_in_NSVM', 'auc_out_NSVM',
                       'auc_recv_ABLR','auc_send_ABLR')
print(results)