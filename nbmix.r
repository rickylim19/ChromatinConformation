
#function to fit mixture of negative binomial (Poisson with Gamma overdispersion) distribution
# K = number of components
fit.nbmix <- function(data,K) {
   lambda <- rep(0,K-1)
   prob.init   <- seq(0.2,0.8,len=K)
   lsize.init <- log(rep(round(mean(data)),K))
   beta.init <- log(prob.init/(1-prob.init))
   par.init  <- c(lambda,beta.init,lsize.init)
    
   # fit model
   iter <- 1
   cat('Iteration ',iter,'\n',sep="")
   out <- optim(p=par.init,fn=nlogl.nbmix,data=data,K=K,method='Nelder-Mead')
   if(out$conv!=0) {
      out <- optim(p=out$p,fn=nlogl.nbmix,data=data,K=K,method='Nelder-Mead')
   }
   while(out$conv!=0) {
     iter <- iter+1
     cat('Iteration ',iter,'\n',sep="")
     out <- optim(p=out$p,fn=nlogl.nbmix,data=data,K=K,method='L-BFGS-B')
   }
   # collate results
   lambda = c(0,out$par[1:(K-1)])
   prop = exp(lambda)/sum(exp(lambda))

   prob = exp(out$par[K:(2*K-1)])/(1+exp(out$par[K:(2*K-1)]))
   size = exp(out$par[(2*K):(3*K-1)])

   AIC = 2*out$v + 2*length(out$p)
   BIC = 2*out$v + length(out$p)*log(length(data))

   result = list(prop=prop,prob=prob,size=size,AIC=AIC,BIC=BIC)

   result
}

#negative log-likelihood function for mixture of NB
nlogl.nbmix <- function(p,data,K) {
    
    lambda = c(0,p[1:(K-1)])
    prop = exp(lambda)/sum(exp(lambda))
    prob = exp(p[K:(2*K-1)])/(1+exp(p[K:(2*K-1)]))
    size = exp(p[(2*K):(3*K-1)])

    lik = t(apply(as.matrix(data),1,dnbinom,prob=prob,size=size)) %*% as.matrix(prop)
    #lik = ifelse(lik==0,1e-200,lik)
    logl = sum(log(lik))
    nlogl= - logl
nlogl
}

#class prediction
class.nbmix <- function(fit,data) {
   K = length(fit$prop)
   prop = fit$prop
   prob = fit$prob
   size = fit$size

   mar.prob = t(apply(as.matrix(data),1,dnbinom,prob=prob,size=size)) %*% as.matrix(diag(prop))
   mar.prob = mar.prob/apply(mar.prob,1,sum)

   max.index <- which(mar.prob==apply(mar.prob,1,max),arr.ind=TRUE)
   class <- rep(0,length(data))
   class[max.index[,1]] <- max.index[,2]

   class
}


## EXAMPLE
# generate data from Poisson and fit mixture of NB
y = c(rpois(100,3),rpois(50,5),rpois(25,10))
fit <- fit.nbmix(data=y,K=3)
cl  <- class.nbmix(fit,data=y)






