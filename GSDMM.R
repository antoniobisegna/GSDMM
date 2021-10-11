# i = {1,..,N_D} ---> 1 <= m <= M ---> i = 1..n
# v = {1,..,N_W} ---> 1 <= n <= N_m ---> j = 1..J_I
# k = {1,..,K}





# sampling
library(extraDistr)

sample_Pi = function(z,k,alpha){
  for (i in 1:length(z)){
    counts = colSums(outer(z,1:k,FUN="=="))
    Pi[i] = gtools::rdirichlet(1,counts+alpha)
  }
  return(Pi)
}


sample_z = function(Pi,X,k){
  z= matrix(0, nrow=nrow(X),ncol=ncol(X))
  for (i in 1:length(z)){                  ######
    for (v in 1:length(unique(as.list(X)))){
      B_X.iv= sample_B(X[i,v],z[i,v]==k,gamma)
      p.z.iv = exp(log(Pi[i]) + log(B_X.iv))
      p.z.iv = p.z.iv/ sum(p.z.iv)
      z[i,v] = rmultinom(n = 1:length(Pi) ,size = 1 ,prob = p.z.iv)
    }
  }
  return(z)
}

sample_B = function(X,z,gamma){
  for(v in 1:length(unique(as.list(X)))){
    for(i in length(z)){
      for(l in length(unique(as.list(X)))){
        if (X[i,l]== 1:v){
          if  (z[i,l] == 1:k){
            count.v = colSums(outer(z[i,l],1:k,FUN="==")) + colSums(outer(X[i,l],1:v,FUN="=="))
          }
        }
          
      }
    }
    B[k,]= rdirichlet(1,gamma+count.v)
  }
  return(B)
}


################# GIBBS ####################

GSDMM = function(G,burnin,thin,K,data,alpha,gamma){
  
  
  # Initialize storage
  
  z=matrix(NA,nrow = G+burnin,ncol = length(data))
  Pi.post=matrix(NA,nrow = K,ncol = length(data))
  B.post=matrix(NA,nrow = G+burnin,ncol = K)
  
  Iterations = burnin+(thin*G)
  
  # Initialize parameters
  
  z[1,]=sample_z(Pi = Pi.post[1,],X =data,K)
  Pi.post[1,]= rep(1/K,K)
  B.post[1,]= sample_B(X =data,z = z[1,],gamma =1 )
  
  # Initialize tmp
  z.tmp=z.post[1,]
  Pi.tmp=Pi.post[1,]
  B.tmp=B.post[1,]
  
  #Progression bar
  pb = txtProgressBar(min = 1, max = Iterations, style = 3)
  
  #Counter for saving data
  g=2
  
  for (iter in 2:Iterations){
    
    ##########################
    ######## SAMPLE z ########
    ##########################
    z.tmp= sample_z(Pi=Pi.tmp,
                    X=data,
                    k=max(unique(z.tmp)))
    if(iter<= burnin | iter%%thin == 0){
      z[g,]=z.tmp
    } 
    
    #########################
    ######## SAMPLE Pi ######
    #########################
    Pi.tmp=sample_Pi(z = z.tmp,
                     k = max(unique(z.tmp)),
                     alpha = alpha)
    if(iter<= burnin | iter%%thin == 0){
      Pi.post[g,]=Pi.tmp
    }
    
    #########################
    ######## SAMPLE B #######
    #########################
    B.tmp=sample_B(X =data,
                   z =z.tmp,
                   gamma = gamma)
    if(iter<= burnin | iter%%thin == 0){
      B.post[g,]=B.tmp
    }
    
    ######################################
    ### MOVING COUNTER FOR SAVING DATA ###
    ######################################
    if(iter<=burnin | iter%%thin == 0){
      g=g+1
    }
    
    #PROGRESS BAR 
    setTxtProgressBar(pb, iter)
  }
  
  close(pb)
  
  results      = list()
  results$B    = B.post
  results$Pi   = Pi.post
  results$z    = z
  return(results)
}

proto_GSDMM  = GSDMM(G=10000,
                     burnin = 1000,
                     thin = 1,
                     K=2,
                     data=data_dirch$samples,
                     alpha = 1,
                     gamma = 1)

View(data_dirch$samples)
