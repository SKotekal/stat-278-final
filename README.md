# stat-278-final
Final project for STAT 278 - Spring 2018

## Introduction
This project is prediction of tumor types from gene expression levels with bounded misclassification error. The procedure developed in ["Least Ambiguous Set-Valued Classifiers with Bounded Error Levels" (Sadinle et al., 2017)](https://arxiv.org/pdf/1609.00451.pdf) is used. We investigate how different methods of feature selection, probability distribution estimation, and data splitting affect various metrics. In particular, we will evaluate the performance of the classifier, the size of null regions (points that have no classification), the size of ambigiuous regions (points that have multiple classifications), and hyperparameter effects. 


## Project Outline

### Partitioning Data
The data will need to be partitioned into a training and test set. Further, the training set will hvae to be further partitioned as part of the classifier procedure, but this will be discussed later. The breakdown of the data into a training and test set will require some thought. In particular, we must consider

* 50/50 split of test and training? 80/20 split? 
* Completely random selection of training and test sets?
* Ensure the same proportion of tumor types in training and test sets?
* Same number of all tumor types in training set?

### Feature Selection 
The gene dataset is very high-dimensional (801 patients, 20264 genes) and it is computationally intractable to construct "least ambiguous with bounded error level" (LABEL) classifiers directly from all of the data. As such, a feature selection step is necessary. Clearly any choice here will affect the final results (such as the size of ambiguous and null regions), and so we will investigate several feature selection methods. In particular, we may investigate 

* LASSO (this will require a choice of &lambda; - perhaps chosen via cross-validation?)
* Orthogonal Matching Pursuit (greedy method - iteratively add features until good fit)
* Correlate genes to tumor type and pick K most higly correlated (choose genes also uncorrelated with other genes?)
* Cluster based on correlation and pick representative from each cluster
* Principal Components Analysis and use first K components as features (Other dim reduct. methods?)
* Variance Threshold - exclude genes that have low variance across all patients (threshold is user-defined parameter)

### Estimating Conditional Probability Distribution 
The procedure laid out in the paper requires estimating the conditional probability distribution ```p(Y|X)``` where ```Y``` is the tumor type and ```X``` is the vector of features. The estimation methods that could be used include

* kNN (Nearest Neighbors - choice of k required)
* Local Polynomial Estimator
* Regularized Logistic Regression
* Kernel Regression (other SVM methods for regression?)
* Naive Bayes

I dont think there should be much of a difference in choosing different methods at this stage since there is a ***conformal inference*** step that guarantees distribution free, finite sample validity. I think the estimation procedure just needs to satisfy some technical requirements given in the paper. However, my understanding may be mistaken here, so we should investigate this point as well. 

### Split-Conformal Inference 
This split-conformal inference method first splits the data into two parts, ```I_1``` and ```I_2```. In this step of the procedure, ```I_1``` is used to estimate the conditonal probability distribution, and ```I_2``` is used to compute the classification region (which is specified by the level sets of the conditional probability distribution). 


As with splitting the data into training and test sets, thought must be put into how ```I_1``` and ```I_2``` are constructed. Again, we must consider

* Proportion of tumor type labels should be the same between ```I_1``` and ```I_2```
* Same number of tumor types (especially if one tumor type is rare in dataset)
* 50/50 split? Other split proportions? We may need more data for estimating ```p(Y|X)``` as compared to ```t```?

The classification region is given by ```{(x, y) : p_hat(y|x) > t_hat}``` where ```p_hat``` and ```t_hat``` are estimates computed from ```I_1``` and ```I_2```. See paper for more details. 

### Assessing Models
The different choices we make must be assessed in order to find an "optimal" classifier. We will compare on the following metrics

* Empirical Misclassification (error ought to be with &alpha; level by split-conformal guarantee - should verify anway)
   * Is there a trend with which ones are misclassified? The &alpha; bound is NOT conditional on a point, but on average.
* Size of null region (number of points with null set prediction)
* Size of ambiguity (number of points with 2 or more class prediction)
* Relationship of above metrics to complexity of model (number of features, etc)
* Relationship of above metrics as &alpha; changes (phase diagram?)

## Further Extensions
### Total Coverage vs Class Specific Coverage
A caveat to this procedure is that the error bound does not guarantee &alpha; error bound conditional on a point. Instead, the error bound is on average. 

Method can be adapted to have &alpha; error bound for each class. This will affect the null and ambiguity regions. We may wish to look into which methods give less ambiguous class specific coverage and which give less ambiguous total coverage. Maybe there is an interesting trade-off. 

Further, there are instances where we wish to control class error more than total error. Consider the following passage from the paper (Sadinle et. al, 2017, pg 8)

> Although it seems reasonable to work with procedures that control the total probability
of an error, in some circumstances this approach may lead to unsatisfactory classifiers.
In particular, when one of the classes is much more prevalent than the others, the
probability of properly labeling an element of the smaller classes may be quite low,
and it decreases as the probability of the largest class increases.


### Dealing with Null Regions 
The authors of the paper give methods to deal with null regions, namely by filling with a baseline classifier and a procedure termed Accretive Completion. What interesting questions are there regarding this?

### Inference
We may also wish to do some inference at the end of the classification step. Specifically, once we have identified a classifier with bounded &alpha; error that best minimizes ambiguity and the null region, we have a good model for which genes are most predictive of tumor type. With these genes, we may be able to do inference and determine which genes are expressed higher/lower given a tumor type. This means that we have turned the classification procedure into a very extensive feature selection procedure, but this may give some guarantees regarding the inference. This is seemingly (I'm probably mistaken) equivalent to "regressing out" effects of other genes to study the relationships between selected genes and tumor type. This is really an extension of the project if needed.  

This will require that we partition the data further into a training set, test set, and inference set, where we use the inference set to actually test our hypothesis. In this paradigm, we can use classical hypothesis testing.

However, if we decide to make it more complicated, we can get rid of an inference set and do inference via selective inference methods. This seems much more complicated as the selection procedure for the genes we do select is very extensive. This approach is probably intractable.