
## Read Data from csv files

labels = read.csv('./TCGA-PANCAN-HiSeq-801x20531/labels.csv')[,-1] # labels[i] = tumor type for patient #i (5 types: BRCA, COAD, KIRC, LUAD, PRAD)
data = read.csv('./TCGA-PANCAN-HiSeq-801x20531/data.csv')[,-1]
data = data[,which(colSums(data!=0)>0)] # removing genes with zero expression level across all patients
colnames(data) = NULL
# 801 patients, 20264 genes

###############
# Set up data #
###############
set.seed(278)
X = as.matrix(data) 
Y = labels
n = length(Y)
tumors <- c("BRCA", "COAD", "KIRC", "LUAD", "PRAD")
NUM_GENES <- dim(X)[2]

num_train <- 640
num_test <- n-num_train

train <- sample(1:n, num_train)
test <- setdiff(1:n, train)

Xtest = X[test,]
Ytest = Y[test]
Xtrain = X[train,]
Ytrain = Y[train]

# Split data for split-conformal inference step
ntrain <- length(train)
I1_ind <- sample(1:ntrain, floor(ntrain/2))
I2_ind <- setdiff(1:ntrain, I1_ind)
I1 <- Xtrain[I1_ind,]
I1y <-  Ytrain[I1_ind]
I2 <- Xtrain[I2_ind,]
I2y <- Ytrain[I2_ind]


num_neighbors <- 10
alpha <- 0.1


# We will use kNN to fit p(y|x)
library(FNN)

kNN_condl_prob <- function(y, x, num_neighbors, xTrain, yTrain){
  nbhrs = get.knnx(xTrain, x, k=num_neighbors) 
  return(length(which(yTrain[nbhrs$nn.index] == y))/num_neighbors)
}

 # Compute thresholds for total coverage
thlds <- function(phat, alpha) {
  S = sort(phat, decreasing=FALSE)
  return(S[ceiling((length(phat)+1)*alpha-1)])
}

# ## Greedily add genes and see ambiguity change
# model_genes <- c()
# genes_used <- rep(0,NUM_GENES)

# # First, we drop genes that have low variance in order to reduce the 
# # the total search space
# feat_vars <- apply(scale(I1, center=FALSE, scale=colSums(I1)), 2, var)
# hist(feat_vars)
# quantile(feat_vars, 0.9, na.rm=TRUE)
# high_var_ind = which(feat_vars >= quantile(feat_vars, 0.99, na.rm=TRUE))
# feat_vars <- feat_vars[which(feat_vars >= quantile(feat_vars, 0.99, na.rm=TRUE))]
# length(feat_vars)


# for(k in 1:3) {
#   min_ambig <- Inf
#   min_null <- Inf
#   best_gene <- 0
#   for(g in high_var_ind){
#     if(genes_used[g] == 1){
#       next
#     }
    
#     # Estimate the distribution of p(y|x)
#     phat_I2 <- c()
#     trial_model <- c(model_genes, g)
#     for(i in 1:length(I2_ind)) {
#       phat_I2 <- c(phat_I2, kNN_condl_prob(I2y[i], t(I2[i, trial_model]), num_neighbors, I1[,trial_model], I1y[trial_model]))
#     }
    
#     # Get thresholds
#     t = thlds(phat_I2, alpha)

#     # Classify train data
#     results = matrix(0L, ncol=5, nrow=length(train))
#     pvals = c()
#     for(i in 1:length(train)) {
#       x = Xtrain[i, trial_model]
      
#       index = 1
#       for(tum in tumors){
#         p <- kNN_condl_prob(tum, t(x), num_neighbors, I1[, trial_model], I1y[trial_model])
#         pvals <- c(pvals, p)
#         if(p >= t){
#           results[i, index] = 1
#         }
#         index <- index+1
#       }
#     }
#     # Count number of points with null classification
#     num_null <- length(which(rowSums(results) == 0))
#     #print(num_null)
  
#     # Count number of points with more than 1 classification
#     num_ambig <- length(which(rowSums(results) > 1))
#     #print(num_ambig)
    
#     if(num_null < min_null){
#       min_null <- num_null
#     }
#     if(num_ambig < min_ambig){
#       min_ambig <- num_ambig
#       best_gene <- g
#     }
#   }
#   print(best_gene)
#   print(min_ambig)
#   print(min_null)
#   model_genes <- c(model_genes, best_gene)
#   print(model_genes)
#   genes_used[best_gene] <- 1
# }

# First, we drop genes that have low variance in order to reduce the 
# the total search space
feat_vars <- apply(scale(I1, center=FALSE, scale=colSums(I1)), 2, var)
hist(feat_vars)
high_var_ind = which(feat_vars >= quantile(feat_vars, 0.99, na.rm=TRUE))
feat_vars <- feat_vars[which(feat_vars >= quantile(feat_vars, 0.99, na.rm=TRUE))]
length(feat_vars)


# Estimate the distribution of p(y|x)
phat_I2 <- c()
trial_model <- feat_vars
for(i in 1:length(I2_ind)) {
  phat_I2 <- c(phat_I2, kNN_condl_prob(I2y[i], t(I2[i, trial_model]), num_neighbors, I1[,trial_model], I1y[trial_model]))
}

# Get thresholds
t = thlds(phat_I2, alpha)

# Classify train data
results = matrix(0L, ncol=5, nrow=length(train))
pvals = c()
for(i in 1:length(train)) {
  x = Xtrain[i, trial_model]
  
  index = 1
  for(tum in tumors){
    p <- kNN_condl_prob(tum, t(x), num_neighbors, I1[, trial_model], I1y[trial_model])
    pvals <- c(pvals, p)
    if(p >= t){
      results[i, index] = 1
    }
    index <- index+1
  }
}
# Count number of points with null classification
num_null <- length(which(rowSums(results) == 0))
print(num_null)

# Count number of points with more than 1 classification
num_ambig <- length(which(rowSums(results) > 1))
print(num_ambig)
