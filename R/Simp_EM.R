Simp_EM <-
function(S,n,p_n,a0,b0,alpha,beta,v0,v1,maxiter,w,l){
  KL1<-NULL
  Q3=NULL
  o=0
  Q1=-Inf
  Q2=100000
  sigma2=0.01
  ###initial###
  j=0
  p=0.5
  tau1=100
  tau=0
  Theta1=100
  W=S
  if(p_n>=n){
    diag(W)=diag(S)+p_n/lambda
  }
  Theta=solve(W)
  
  id=which(upper.tri(Theta)==1)
 
  Theta1=0
  sigma1=100
  tau1=100
  #Entropy1=-1000
  P=matrix(0.5,nrow=p_n,ncol=p_n)
  while(j<maxiter&(sum(abs(Theta-Theta1))>0.001|abs(sigma2-sigma1)>0.001|abs(tau1-tau)>0.001)){
    Q1=-Inf
    sigma1=sigma2
    tau1=tau
    #update posterior for r

    P=1/(1+sqrt(v1/v0)*exp(-Theta^2/(2*v0*sigma2)+Theta^2/(2*v1*sigma2))*(1-p)/p)
    
    #P[P>0.9995]=0.9995
    #P[P<0.0005]=0.0005
    #P=(dnorm(Theta,mean=0,sd=sqrt(v1*sigma2))*p)/
    # (dnorm(Theta,mean=0,sd=sqrt(v1*sigma2))*p+
    #   dnorm(Theta,mean=0,sd=sqrt(v0*sigma2))*(1-p))
    
    #Entropy=-(1-P)*log(1-P)-P*log(P)
    #if(max(abs(Entropy1-Entropy))<0.0001){
    #  break;
    #}
    #Entropy1=Entropy
    
    #P=(dnorm(Theta,mean=0,sd=sqrt(v1*sigma2))*p)/
     # (dnorm(Theta,mean=0,sd=sqrt(v1*sigma2))*p+
      #   dnorm(Theta,mean=0,sd=sqrt(v0*sigma2))*(1-p))
 
    ###update p#####   
    p=(a0-1+sum(P[id]))/(a0+b0-2+p_n*(p_n-1)/2)
    
    ##update sigma2
    sigma2=(sum((P*Theta^2/(2*v1)+(1-P)*Theta^2/(2*v0))[
      id])+w*l/2)/(p_n*(p_n-1)/4+w/2+1)
 
#####ensure not achieve extreme value######
    if(p<0.00001){
      p=0.00001
    }else if(p>0.99999){
      p=0.99999
    }
    
      
    tau1=100
    Theta1=Theta
    Theta2=1000
q=0
    while(abs(tau1-tau)>0.001&q<=1){
      
      q=q+1
      tau1=tau
      ###update tau#####
      tau=(alpha-1+p_n)/(sum(diag(Theta))+beta)
      Theta2=Theta
      for(i in 1:p_n){ 
        Theta3=Theta
        W3=W
        #ensure positive definite

        if(tau>0){
          W[i,i]=S[i,i]+2/n*tau
        }
        ###update Theta_12#####
       Theta[-i,i]=Theta[i,-i]=
          -Theta[-i,-i]%*%matrix(W[i,-i])/W[i,i]
        
    #  Theta[-i,i]=Theta[i,-i]=v*0.001+Theta[i,-i]
       
        Theta[i,i]=(1-W[i,-i]%*%matrix(Theta[i,-i]))/W[i,i]
      # W[-i,i]=W[i,-i]=
       #  -W[-i,-i]%*%matrix(Theta[i,-i])/Theta[i,i]
    Q2=Q4(Theta,tau,alpha,beta,a0,b0,v0,v1,p,P,sigma2,id)
    if(Q2<(Q1)){
      Theta=Theta3
      W=W3
      o=o+1
      Q2=Q1
    }
    Q1=Q2

      }

      
    }
    
    
    j=j+1
    print(j)
  }
P=1/(1+sqrt(v1/v0)*exp(-Theta^2/(2*v0*sigma2)+Theta^2/(2*v1*sigma2))*(1-p)/p)
  return(list(Theta=Theta,P=P,Q=Q2,tau=tau,p=p,W=W,sigma2=sigma2))
}
