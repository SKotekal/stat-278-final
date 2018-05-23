## Probability Estimation Methods
# Not sure how important these methods are, given that the
# paper achieves distribution-free, finite sample validity
# via split-conformal inference. 

## k Nearest Neighbors
# Estimating p(Y|X) using kNN method with num_neighbors specified
# by the user. In particular, given y and x and the associated
# training data, returns the conditional probability.
kNN_condl_prob <- function(y, x, num_neighbors, xTrain, yTrain){
  library(FNN)
  nbhrs = get.knnx(xTrain, x, k=num_neighbors) 
  return(length(which(yTrain[nbhrs$nn.index] == y))/num_neighbors)
}

## Regularized Logistic Regression
# Estimating p(Y|X) using Logistic LASSO 
# We give the code below to copy and paste into 
# main method, since it is not efficient to keep
# fitting logistic lasso every time to evaluate test point on model

# library(glmnet)
# cvfit = cv.glmnet(Xtrain, Ytrain, family="multinomial", type.multinomial="grouped")
# cndl_probs = predict(cvfit, newx=I2 s="lambda.min", type="response")
# phat_I2 <- c()
# for (i in 1:length(I2_ind)) {
#     phat_I2 <- c(phat_I2, cndl_probs[i, I2y[i], 1])
# }

## Local Polynomial Estimator
# Estimating p(Y|X) via a local polynomial estimator (see Sadinle et al.
# and Tsyabkov, 2009 - Introduction to Nonparametric Estimation)
# We provide some kernels below
kernel_prob <- function(y, x, xTrain, yTrain, kernel, bandwith){
    cterm <- 1.0/(length(yTrain)*(bandwith**2))
    prob <- 0
    for(i in 1:length(yTrain)){
        if(kernel == "rect"){
            prob = prob + as.double(y == yTrain[i])*rect_kernel(x, xTrain[i,], bandwith)
        }
        if(kernel == "gaussian"){
            prob = prob + as.double(y == yTrain[i])*gauss_kernel(x, xTrain[i,], bandwith)
        }
    }
    return(prob/cterm)
}

# we use the Parzen-Rosenblatt estimator and (see pg 4 http://ins.sjtu.edu.cn/files/common/20121209191850_7.pdf
# in Tsykabov 2009)
rect_kernel <- function(x, x0, h){
    return(prod(0.5*as.double(abs((x-x0)/h) < 1)))
}
gauss_kernel <- function(x, x_0, h){
    return(prod(dnorm((x-x_0)/h)))
}
