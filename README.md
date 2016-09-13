# BayesEM
R package (Bayes EM algorithms for Gaussian Graphical Models)

##Maintainer
Lingrui(Gary) Gan, University of Illinois at Urbana and Champaign

##Description
Fast Bayesian algorithms for learning the sparse structure of a Gaussian graphical model

##Example
###Simulated Samples from AR(1) Structure
```{r}
C=toeplitz(c(1,0.5,rep(0,p_n-2)))
Sigma=solve(C)
Y<-mvrnorm(n,rep(0,p_n),Sigma)
S<-cov(Y) 
```

###Learning the Sparse Structure From BayesEM Algorithms
####Hyperparameters
```{r}
a0=1
b0=p_n
alpha=1000
beta=10
alpha1=1
beta1=100
v0=1
v1=5
v11=5
w=1
l=1

maxiter=30
maxiter1=30
```
####Implement BayesEM Algorithm for the Sparse Learning
```{r}
result<-EM(S,n,p_n,a0,b0,alpha1,beta1,v0,v11,maxiter1,w,l) #full
```

####Implement Simplified BayesEM Algorithm for the Sparse Learning
```{r}
result<-EM1(S,n,p_n,a0,b0,alpha,beta,v0,v1,maxiter,w,l) #simplified version
```

####Outputs 
#####The marginal posterior probability matrix of each entries
```{r}
result$P
```
#####The MAP estimate of the precision matrix 
```{r}
result$Theta
```
##Reference
Manuscript  
