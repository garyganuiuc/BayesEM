Tune_EMLasso <-
function(S,n,p_n,a0,b0,alpha,beta,v0,maxiter,w,l,Ra){
  pb <- txtProgressBar(min = 0, max = length(v0)*10*length(tau), style = 3)
  bic=NULL
  for(m in 1:length(beta)){
    for(i in 1:length(v0)){
      v1=seq(v0[i]+0.1,5*v0[i],0.5*v0[i])
      for(j in 1:length(v1)){
        for(k in 1:length(Ra)){
          ##Bayes EM####
          #v1=sqrt(n/(log(p_n)))/5
          #v0=sqrst(n/(p_n*log(p_n)))/5
          
          #tau=sqrt((log(p_n)/n))
          w=1
          l=1
          maxiter=30
          result1<-EM(S,n,p_n,a0,b0,alpha,beta[m],v0[i],v1[j],maxiter,w,l,Ra[k])
          
          bic=rbind(bic,list(v0=v0[i],v1=v1[j],Ra=Ra[k],beta=beta[m],BIC=BIC_EMLasso(result1$Theta,S,result1$P,n)))
          
          setTxtProgressBar(pb, (m-1)*length(v0)*length(v1)*length(tau)+(i-1)*length(v1)*length(tau)+(j-1)*length(tau)+k)
          
        }
      } 
    }
    
  }
  
  close(pb)
  return(bic[which.min(bic[,5]),])
}
