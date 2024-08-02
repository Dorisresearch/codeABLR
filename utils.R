#### adjacency matrix ####
adj_matrix <- function(p,n,n_sample){
  A <- list()
  for (i in 1:n_sample){
    temp <- matrix(rbinom(length(p),1,p), nr=n, nc=n)
    A[[i]] <- temp
  } 
  return(A)
}

#### weighted adjacency matrix ####
model.weighted.sparse <- function(rep = m,a1,b1,a2,b2,a3,b3,a4,b4,delta=0.05,n,noloop = TRUE){
  r1 <- n/2
  r2 <- n - n/2
  Gra1 <- Gra2 <- list()
  for (i in 1:rep) {
    # normal distribution
    G1 <- matrix(rnorm(r1^2,a1,b1) ,nrow = r1)
    G4 <- matrix(rnorm(r2^2,a4,b4) ,nrow = r2)
    G3 <- matrix(rnorm(r1 * r2,a3,b3) ,nrow = r2)
    G2 <- matrix(rnorm(r1 * r2,a2,b2) ,nrow = r1)
    
    F1 <- matrix(rnorm(r1^2,a1,b1) ,nrow = r1)
    F4 <- matrix(rnorm(r2^2,a4,b4) ,nrow = r2)
    F3 <- matrix(rnorm(r1 * r2,a3,b3) ,nrow = r2)
    F2 <- matrix(rnorm(r1 * r2,a2,b2) ,nrow = r1)
    
    tmpG <- cbind(rbind(G1,G3),rbind(G2,G4))
    tmpF <- cbind(rbind(F1,F3),rbind(F2,F4))
    
    # row and col change
    d1 <- d2 <- ceiling(n/4)
    d3 <- ceiling(n/4)
    tmpF[1:(n/2),(n/2 + 1):(n/2 + d2)] <- matrix(rnorm(d2 * n / 2,a2 + delta,b2), ncol = d2)
    tmpF[(n/2 + 1):(n/2 + d3),(n/2 + 1):(n/2 + d2)] <- matrix(rnorm(d2 * d3,a4 - delta,b4), ncol = d2)
    tmpF[1:d1,1:(n/2)] <- matrix(rnorm(d1 * n/2,a1 - delta,b1) ,nrow = d1)
    
    if (noloop) {
      diag(tmpG) <- diag(tmpF) <- 0
    }
    Gra1[[i]] = tmpG
    Gra2[[i]] = tmpF
  }
  G <- list()
  G$Gra1 <- Gra1
  G$Gra2 <- Gra2
  return(G)
}

#### FSVM ####
# "Directed brain network analysis in anxious and non-anxious depression based on EEG source reconstruction and graph theory" [2023, Hesam] 
# library(igraph)
# library(caret)  data splitting
# A1,A2: two lists of network adjacency matrix
# label: label of A1 and A2
# feature: inward,outward, BC
# parts: training part id
Hesam_SVM <- function(A1,A2,label,feature,BC_1 = NA,BC_2 = NA,parts){
  n <- dim(A1[[1]])[1]
  # feature vector 1
  if (feature == 'inward'){
    inward_1 <- lapply(A1, function(x){colSums(x)})
    inward_1 <- matrix(unlist(inward_1), nrow = length(A1), byrow = TRUE)
    
    inward_2 <- lapply(A2, function(x){colSums(x)})
    inward_2 <- matrix(unlist(inward_2), nrow = length(A2), byrow = TRUE) 
    
    feature_1 <- inward_1
    feature_2 <- inward_2
  }else if (feature == 'outward'){  # feature vector 2
    outward_1 <- lapply(A1, function(x){rowSums(x)})
    outward_1 <- matrix(unlist(outward_1), nrow = length(A1), byrow = TRUE)
    
    outward_2 <- lapply(A2, function(x){rowSums(x)})
    outward_2 <- matrix(unlist(outward_2), nrow = length(A2), byrow = TRUE)
    
    feature_1 <- outward_1
    feature_2 <- outward_2
  }else  if (feature == 'BC'){  # feature vector 3: compute betweenness centrality (BC)
    # translate adjacency matrix to graph object
    g_1 <- lapply(A1, function(x){graph_from_adjacency_matrix(as.matrix(x))})
    BC_1 <- lapply(g_1, function(x){betweenness(x)})
    BC_1 <- matrix(unlist(BC_1), nrow = length(A1), byrow = TRUE)
    
    g_2 <- lapply(A2, function(x){graph_from_adjacency_matrix(as.matrix(x))})
    BC_2 <- lapply(g_2, function(x){betweenness(x)})
    BC_2 <- matrix(unlist(BC_2), nrow = length(A2), byrow = TRUE)
    
    feature_1 <- BC_1
    feature_2 <- BC_2
  }else if(feature == 'custom'){
    feature_1 <- BC_1
    feature_2 <- BC_2
  }
  
  ## Mann-Whitney U test for each feature
  # select nodes using inward/outward/BC feature vector
  sel_in <- c()
  for (i in c(1:n)){
    res <- wilcox.test(feature_1[,i],feature_2[,i])
    p_value <- res$p.value
    if (p_value < 0.05){
      sel_in <- c(sel_in, i)
    }
  }
  
  ## SVM 
  if (length(sel_in) == 0){
    acc <- F1 <- spe <- auc_value <- table <- 0
  }else{
    feature <- rbind(feature_1[,sel_in],feature_2[,sel_in])
    data <- data.frame(feature=feature,cl = as.factor(label))
    
    data_train = data[parts, ]
    data_test = data[-parts, ]
    # scaling data
    d <- dim(data)[2]
    data_train[,-d] <- scale(data_train[,-d])
    data_test[,-d] <- scale(data_test[,-d])
    
    # leave one out cross-validation method
    train_control <- trainControl(method = "LOOCV")
    # training
    model <- train(cl~., data = data_train, trControl = train_control, method = "svmLinear")
    # use model to make predictions on test data
    pred_y = predict(model, data_test)
    
    if (all(data_test$cl == 0)|all(data_test$cl == 1)){
      acc <- sum(pred_y == data_test$cl) / length(data_test$cl)
      F1 <- spe <- auc_value <- table <- NA
    }else{
      # confusion Matrix
      result <- confusionMatrix(data = pred_y, data_test$cl,mode = "everything",positive="1")
      table <- result$table
      acc <- result$overall['Accuracy']  # Accuracy
      auc_value <- auc(data_test$cl,pred_y) # calculate AUC
    }
  }
  
  # outputs
  res <- vector("list")
  res$acc <- acc
  res$auc <- auc_value
  res$table <- table
  res$nodes_select <- sel_in
  return(res)
}


#### ABLR ####
### new with intercept
loss_p_new <- function(A_tuta_train,label_A_tuta_train,b,s,a0,lambda_0,lambda_1,eta = 0.5){
  loss <- 0
  p <- c()
  n_sample <- length(A_tuta_train)
  for (i in (1:n_sample)){
    temp2 <- t(s) %*% A_tuta_train[[i]] %*% b + a0
    temp <- temp2 * label_A_tuta_train[i] - log1pexp(temp2)
    loss <- loss + temp
    
    p[i] <- 1 - 1 / (1 + exp(temp2))
  }
  loss_old <- - loss / n_sample + 
    lambda_0 * (norm(b, type = '2'))^2  + 
    lambda_1 * (1 - eta) * (norm(s, type = '2'))^2 + lambda_1 * eta * sum(abs(s)) 
  
  res <- list()
  res$loss <- loss_old
  res$p <- p
  return(res)
}

p_ud_new <- function(A_tuta_train,s,b,a0){
  p <- c()
  for (i in 1:length(A_tuta_train)){
    val1 <- t(s) %*% A_tuta_train[[i]] %*% b + a0
    p[i] <- 1 - 1 / (1 + exp(val1))
  }
  return(p)
}

Dir_LOGI_BCD_new <- function(A_tuta_train,A_tuta_test,label_A_tuta_train,label_A_tuta_test,
                             s,b,a0,lambda_0,lambda_1,eta = 0.5,
                             tol = 1e-5, count = 1, count_limit = 20000){
  # sample size of training data
  n_sample <- length(A_tuta_train)
  # loss_old 
  loss_p_res <- loss_p_new(A_tuta_train,label_A_tuta_train,b,s,a0,lambda_0,lambda_1,eta)
  change <-  loss_old <- loss_p_res$loss
  # prob vector
  p <- loss_p_res$p
  while (change > tol) {
    gamma <- t <- min(0.01,(1/(count + 1))^0.5)
    dif <- label_A_tuta_train - p
    # update s
    val1 <- lapply(A_tuta_train, function(x){x %*% b})
    val2 <- mapply("*",val1,dif,SIMPLIFY = F)
    val3 <- - Reduce("+",val2) / n_sample + 2 * lambda_1 * (1 - eta) * s
    s <- S_func(s - gamma * val3, lambda_1 * eta * t)
    
    # update b
    p <- p_ud_new(A_tuta_train,s,b,a0)
    dif <- label_A_tuta_train - p
    val1 <- lapply(A_tuta_train, function(x){t(x) %*% s})
    val2 <- mapply("*",val1,dif,SIMPLIFY = F)
    val3 <- - Reduce("+",val2) / n_sample + 2 * lambda_0 * b
    b <- b - gamma * val3
    
    # update a0
    p <- p_ud_new(A_tuta_train,s,b,a0)
    dif <- label_A_tuta_train - p
    val1 <- sum(dif) / n_sample 
    a0 <- a0 + gamma * val1
    
    loss_p_res <- loss_p_new(A_tuta_train,label_A_tuta_train,b,s,a0,lambda_0,lambda_1,eta)
    loss_new <- loss_p_res$loss
    p <- loss_p_res$p
    change <- abs(loss_new - loss_old)/loss_new
    count <- count + 1
    if (count > count_limit){
      break
    }
    loss_old <- loss_new
  }
  
  # use model to make predictions on test data
  pred_logi <- p_ud_new(A_tuta_test,s,b,a0)
  pred_logi <- ifelse(pred_logi >= 0.5, 1, 0)
  # confusion Matrix
  if (length(unique(label_A_tuta_test)) == 1){
    acc <- sum(pred_logi == label_A_tuta_test) / length(label_A_tuta_test)
    F1 <- spe <- auc_value <- table <- NA
  }else{
    result <- confusionMatrix(data = as.factor(pred_logi), as.factor(label_A_tuta_test),mode = "everything",positive="1")
    table <- result$table
    acc <- result$overall['Accuracy']  # Accuracy
    auc_value <- auc(label_A_tuta_test,pred_logi ) # calculate AUC
  }
  # outputs
  res <- vector("list")
  res$table <- table
  res$s <- s
  res$b <- b
  res$a0 <- a0
  res$rel_change <- change
  res$loss <- loss_new
  res$acc <- acc
  res$auc <- auc_value
  res$k <- count
  res$lambda_0 <- lambda_0
  res$lambda_1 <- lambda_1
  res$pre_label <- pred_logi
  return(res)
}