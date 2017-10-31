# BayesEM
R package (Bayes EM algorithms for Gaussian Graphical Models)

##Maintainer
Lingrui(Gary) Gan, University of Illinois at Urbana and Champaign

Feng Liang, University of Illinois at Urbana and Champaign
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
b0=1
alpha=1
beta=c(1,50,100)

  
v0=0.1
v1=c(0.3,0.4,0.5)

maxiter=100
Ra=3000
```
#### Tune the parameters
```{r}
Tune=Tune_EMLasso(S,n,p_n,a0,b0,alpha,beta,v0,maxiter,w,l,Ra)
Ra=Tune$Ra
v1=Tune$v1
beta=Tune$beta
```
####Implement BayesEM Algorithm for the Sparse Learning
```{r}
result<-EM(S,n,p_n,a0,b0,alpha,beta,v0,v1,maxiter,w,l,Ra)
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
