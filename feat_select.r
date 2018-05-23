# Feature selection methods

## LASSO
# Use the glmnet package - contains a cross-validated method version as well
# Example code is below

# library(glmnet)
# cvfit = cv.glmnet(Xtrain, Ytrain, family="multinomial", type.multinomial="grouped")

# # Grab lambda values and associated selected variables
# lambda_vals <- cvfit$lambda
# selected_var <- predict(cvfit,type="nonzero")
# lamb_index = cvfit$lambda.1se


## Fisher Score 
# Features should assign similar values to individuals with the 
# same tumor type, and different values to individuals with different tumor types
# We compute a Fisher Score for each feature i. Then returns the top K features
# See https://pdfs.semanticscholar.org/310e/a531640728702fce6c743c1dd680a23d2ef4.pdf
fisher_features <- function(X, Y, k) {
    library(plyr)
    scores <- rep(0, dim(X)[2])
    tumors <- 1:5
    
    # compute descriptive statistics
    counts <- as.double(count(Y)[,2])
    feat_means <- as.matrix(colMeans(X))
    class_feat_means <- matrix(0, nrow=length(tumors), ncol=dim(X)[2])
    class_feat_var <- matrix(0, nrow=length(tumors), ncol=dim(X)[2])
    for(i in 1:dim(X)[2]) {
        for(t in tumors){
            class_feat_means[t, i] = mean(X[which(as.double(Y)==t), i])
            class_feat_var[t, i] = var(X[which(as.double(Y)==t), i])
        }
    }
    
    # compute fisher score 
    for(i in 1:dim(X)[2]) {
        scores[i] <- sum(counts*(class_feat_means[,i]-feat_means[i])**2)/sum(counts*(class_feat_var[,i])**2)
    }
    
    return(order(scores[1:k]))

}

## Variance Threshold
# The core idea is that low variance features do not give much information
# regarding the class type. Thus, we drop features that have variance less or equal to
# some user-specified threshold. Note that normalization is necessary as variance
# is not scale invariant. 
variance_threshold <- function(X, thold) {
    feat_vars <- apply(scale(X, center=FALSE, scale=colSums(X)), 2, var)
    return(which(feat_vars > thold))
}

## Correlation Feature Selection 
# General idea is that a good set of features has
# features that are highly correlated with the class 
# but are uncorrelated with each other. 
# To this end we calculate Spearman rho correlation of each gene with each classification type
# then pick the k features with highest correlation for each class. This results in a total of 5*k features
# There are many correlation measures that may be used, 
# but for simplicity, we will use Spearman's rho (simple, nonparametric). 
# There are many extensions to this very simple 
# feature selection method; see below.
# See https://www.cs.waikato.ac.nz/~mhall/thesis.pdf
# and https://en.wikipedia.org/wiki/Feature_selection#Correlation_feature_selection
correlation_features <- function(X, Y, k) {
    tumors <- 1:5
    feat_class_corr <- matrix(0, nrow=length(tumors), ncol=dim(X)[2])
    for(i in 1:dim(X)[2]){
        for(t in tumors){
            feat_class_corr[t, i] = cor(X[, i], as.double(Y), method=c("spearman"))
        }
    }

    features <- c()
    for(t in tumors){
        features <- c(features, order(feat_class_corr[t,])[1:k])
    }

    return(features)
}

## OMP (Greedy) 
# How does adding features progressively affect
# the size of the null regions and ambiguity? 
# We have a theoretical guarantee of (1-\alpha) coverage
# so not clear that OMP for subset selection really makes 
# sense (minimizing residuals? when we already have (1-\alpha) coverage?)
# We could implement with a different black box model (this will have to be 
# the case with any wrapper method).
# We implement this with user-specified model size (k)
