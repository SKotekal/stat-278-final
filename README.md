# stat-278-final
Final project for STAT 278 - Spring 2018

## Introduction
This project is prediction of tumor types from gene expression levels with bounded misclassification error. The procedure developed in ["Least Ambiguous Set-Valued Classifiers with Bounded Error Levels" (Sadinle et al., 2017)](https://arxiv.org/pdf/1609.00451.pdf) is used. We investigate how different methods of feature selection, probability distribution estimation, and data splitting affect various metrics. In particular, we will evaluate the performance of the classifier, the size of null regions (points that have no classification), the size of ambigiuous regions (points that have multiple classifications), and hyperparameter effects. 


## Project Outline

### Objectives/Questions we asked
* How do different feature selection methods affect the results of the classifier?
    * Feature selection can be approached from a different perspectives. In particular, we will be analyzing using LASSO, Fisher Score, and PCA to select features. 
* How does varying &alpha; (the misclassification error bound) affect the ambiguity of the classifier when different methods for feature selection are used?
* The classifier utilizes a split-conformal inference step in which the data must be split for estimating a conditional probability distribution and learning a decision rule. How does the proportion of the split affect the number of ambiguious classifications?
* For a desired ambiguity level, what proportion does our data split have to be? How does this depend on our feature selection method?
* The classifier can be modified to bound the error of misclassification conditional on tumor type. As such, the classifier must be more conservative with its predictions. How conservative does it become with various alpha and data splitting choices?
* How does this method compare with other standard classification algorithms? In particular, how does this compare with LASSO? Is LASSO prone to overfitting? The classifer discussed in the above paper guarantees an error bound, but has a tradeoff with ambiguity. 
* Which genes decrease the ambiguity the most? In particular, we will run OMP with the residual being the number of ambiguious classifications. 


### Partitioning Data
The data will need to be partitioned into a training and test set. Further, the training set will hvae to be further partitioned as part of the classifier procedure, but this will be discussed later. The breakdown of the data into a training and test set will require some thought. In particular, we must consider

* 50/50 split of test and training? 80/20 split? 
* Completely random selection of training and test sets?
* Ensure the same proportion of tumor types in training and test sets?
* Same number of all tumor types in training set?

### Feature Selection 
The gene dataset is very high-dimensional (801 patients, 20264 genes) and it is computationally intractable to construct "least ambiguous with bounded error level" (LABEL) classifiers directly from all of the data. As such, a feature selection step is necessary. Clearly any choice here will affect the final results (such as the size of ambiguous and null regions), and so we will investigate several feature selection methods. In particular, 

* LASSO (&lambda; will be chosen via cross validation) - this is a standard method for feature selection
* Fisher Score - common filter method for feature selection
* Variance Threshold - exclude genes that have low variance across all patients (threshold is user-defined parameter)
* Principal Components Analysis and use first K components as features - common dimensionality reduction method

These methods are chosen to reflect a few approaches to feature selection, and we are interested in how these approaches will differ in their results.

### Estimating Conditional Probability Distribution 
The procedure laid out in the paper requires estimating the conditional probability distribution ```p(Y|X)``` where ```Y``` is the tumor type and ```X``` is the vector of features. The estimation methods that will be considered

* kNN (Nearest Neighbors - choice of k required)
* Regularized Logistic Regression

This step is not entirely critical since there is a ***conformal inference*** step that guarantees distribution free, finite sample validity. The estimation procedure just needs to satisfy some technical requirements given in the paper.

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
* Relationship of above metrics as &alpha; changes (phase diagram)

A short remark regarding the null region. The null region is the region of space in which there is no classification. Conceptually, this is no different than a fully ambiguous classification, in which the classifier gives all 5 tumor types to a point in the null region. These null regions arise due to a technical constraint in the paper. The authors address this null region via methods to extend the classifier, so it is our opinion that is intellectually dishonest to call points that initially have no classification as ambiguous. Yet we still believe it is an interesting question to ask about how these null regions change with different parameters before handling the regions with methods described in the paper. 

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
The authors of the paper give methods to deal with null regions, namely by filling with a baseline classifier and a procedure termed Accretive Completion.