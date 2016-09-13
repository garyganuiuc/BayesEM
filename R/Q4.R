Q4 <-
function(Theta,tau,alpha,beta,a0,b0,v0,v1,p,P,sigma2,id){
  p_n=nrow(Theta)

    p_n*log(tau)-tau*sum(diag(Theta))+
    (a0-1)*log(p)+(b0-1)*log(1-p)+
    (alpha-1)*log(tau)-beta*tau+
    sum((-1/2*log(2*pi*v1*sigma2)-Theta[id]^2/(2*v1*sigma2))*P[id]+
          (-1/2*log(2*pi*v0*sigma2)-Theta[id]^2/(2*v0*sigma2))*(1-P[id])+
          P[id]*log(p)+(1-P[id])*log(1-p))
}
