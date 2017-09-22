#=========================================================================
# program: function.R
# note   : code submission for paper "Integrative Sparse Principal Component 
#                                    Analysis of Gene Expression Data".
# In this file:
# Vmat_value,  update_V_wide,  vgl_mat  define the functions for iSPCA
# SVmat_value, Supdate_V_wide, Svgl_mat define the functions for iSPCA_M 
# SVmat_value_sign, Supdate_V_wide_sign, Svgl_mat_sign define the functions 
# for iSPCA_S  
# lVmat_value, lupdate_V_wide, vl_mat define the functions for SPCA
#


#==========================================================================
# Function: Vmat_value
# Source  : update_V_wide, vgl_mat
# Input   :
#          X         as a matrix of M datasets
#          Lambda    as a vector of the tuning parameters
#          p         as the dimension
#          M         as the number of datasets
# Output  :
#          A list of the first PC loadings at each tuning parameters for iSPCA
          
#==========================================================================
# Function: SVmat_value
# Source  : Supdate_V_wide , Svgl_mat
# Input   :
#          X         as a matrix of M datasets
#          Lambda1   as a vector of the tuning parameters for penalty 1
#          Lambda2   as a vector of the tuning parameters for penalty 2
#          p         as the dimension
#          M         as the number of datasets
# Output  :
#          A list of the first PC loadings at each tuning parameters for iSPCA_M

#==========================================================================
# Function: SVmat_value_sign
# Source  : Supdate_V_wide_sign, Svgl_mat_sign
# Input   :
#          X         as a matrix of M datasets
#          Lambda1   as a vector of the tuning parameters for penalty 1
#          Lambda2   as a vector of the tuning parameters for penalty 2
#          p         as the dimension
#          M         as the number of datasets
# Output  :
#          A list of the first PC loadings at each tuning parameters for iSPCA_S

#===========================================================================
# Function: lVmat_value
# Source  : lupdate_V_wide, vl_mat 
# Input   :
#          X         as a matrix of M datasets
#          Lambda    as a vector of the tuning parameters
#          p         as the dimension
#          M         as the number of datasets
# Output  :
#          A list of the first PC loadings at each tuning parameters for SPCA

#============================          iSPCA          ====================
Vmat_value     <- function(X,Lambda,p,M){   
  set.seed(1234)
  n            <- nrow(X)
  V_mat        <- matrix(0,nr=p,nc=M)
  U_mat        <- matrix(0,nr=n,nc=M)
  m=1
  X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
  V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d
  U_mat[,m]  <- X_m_SVD$u 
  angle <- function(x,y){
    # compute the arc-cosine between x and y
    dot.prod <- x%*%y 
    norm.x <- norm(x)
    norm.y <- norm(y)
    theta <- acos(dot.prod /(norm.x*norm.y))*180/pi
    as.numeric(theta)
  } 
  # Initialize
  for(m in 2:M){
    X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
    if(angle(t(X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))<=angle(t((-1)*X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))){
      s=2
    }else{s=1}
    V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d*(-1)^s
    U_mat[,m]  <- X_m_SVD$u*(-1)^(s)
  }
  
  V_wide_mat   <- matrix(rep(V_mat,length(Lambda)),nr=p)   #####repeat for lambda
  U_wide_mat   <- matrix(rep(U_mat,length(Lambda)),nr=n)
  pos_c        <- 1:length(Lambda)
  while (length(pos_c)>0){
    # update V_wide_mat
    V_wide_mat0                <- V_wide_mat
    pos_long_c                 <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    V_wide_mat[,pos_long_c]    <- update_V_wide(V_wide_mat[,pos_long_c,drop=FALSE],
                                                U_wide_mat[,pos_long_c,drop=FALSE],
                                                X,Lambda[pos_c])
    # update U_wide_mat
    for (m in 1:M){
      deno                       <- sqrt(colSums((as.matrix(X[,((m-1)*p+1):(m*p)])%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE]))^2))+10^(-6)
      deno_mat                   <- matrix(rep(deno,each=n),nr=n)
      U_wide_mat[,M*(pos_c-1)+m] <- as.matrix(X[,((m-1)*p+1):(m*p)])%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE])/deno_mat
      
    }
    
    pos_c                        <-unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-3))!=0)-1)%/%M+1)
    
  }
  # standardize
  V_wide_mat   <- ifelse(abs(V_wide_mat)>10^(-4),V_wide_mat,0)
  v_norm       <- sqrt(colSums(V_wide_mat^2))+0.0001
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)
  V_Dir        <- V_wide_mat/v_mat_norm
  V_Dir        <- ifelse(abs(V_Dir)>10^(-4),V_Dir,0)
  return(list(V_mat=V_Dir))
}

update_V_wide  <- function(V_wide_mat,U_wide_mat,X,Lambda){
  p            <- nrow(V_wide_mat)
  M            <- ncol(V_wide_mat)/length(Lambda)
  v_norm       <- sqrt(colSums(V_wide_mat^2))                    
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)                
  pos_r        <- rep(TRUE,p)
  pos_c        <- 1:length(Lambda)
  n            <- nrow(X)
  for(loop in 1:1){
    V_wide_mat0<- V_wide_mat
    pos_cindex <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    gl_M       <- vgl_mat(V_wide_mat[pos_r,pos_cindex,drop=FALSE],
                       U_wide_mat[,pos_cindex],
                       Lambda[pos_c],n,X) 
    
    V_wide_mat[pos_r,pos_cindex]    <- matrix(gl_M,nr=sum(pos_r))               
    
    pos_c      <- unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-4))!=0)-1)%/%M+1)
    pos_r      <- apply(abs(V_wide_mat-V_wide_mat0),1,max) > 10^(-4) 
    if(sum(pos_c)==0 ) break
    
  }
  return(V_wide_mat)
}

vgl_mat       <- function(V_wide_mat,U_wide_mat,Lambda,n,X){
  p           <- nrow(V_wide_mat)
  M           <- ncol(V_wide_mat)/length(Lambda)
  Lambda_v    <- rep(Lambda,each=M)
  Lambda_v    <- as.matrix(Lambda_v,length(Lambda_v),1)
  Vgl_mat     <- matrix(rep(Lambda_v,each=p),nr=p)
  v           <- list(0) 
  u           <- list(0)
  s           <- matrix(0,nr=p,nc=length(Lambda)*M)
  for (j in 1:p){
    j_index      <- j+p*(c(1:M)-1)
    X_j          <- as.matrix(X[,j_index]) 
    for(m in 1:M){
      X_jm       <- matrix(as.matrix(X_j[,m]),nr=n)
      U_wide_mat <- as.matrix(U_wide_mat)
      U_m        <- U_wide_mat[,m+M*(c(1:length(Lambda))-1)]
      s[j,m+M*(c(1:length(Lambda))-1)] <- (1/n)*t(X_jm)%*%U_m
    }
  }
  s_norm     <-matrix(0,nr=p,nc=length(Lambda)*M)
  for (j in 1:p){
    for(m in 1:M){
      s_norm[j,m+M*(c(1:length(Lambda))-1)]<-sqrt(diag( t(matrix(s[j,],nr=M))%*%matrix(s[j,],nr=M)))
      
    }
  }
  lambda_mat                   <- matrix(rep(Lambda_v,each=p),nr=p)
  s_norm_s                     <- as.vector(s_norm)
  s_norm_s[which(s_norm_s==0)] <- 1
  s_norm_s                     <- matrix(s_norm_s,nr=p)
  Vgl_mat                      <- (s_norm>lambda_mat)*n*s*(1-lambda_mat/s_norm_s)
  return(Vgl_mat)
}



#============================          iSPCA_M        ======================  
SVmat_value    <- function(X,Lambda1,Lambda2,p,M){   
  set.seed(1234)
  n            <- nrow(X)
  V_mat        <- matrix(0,nr=p,nc=M)
  U_mat        <- matrix(0,nr=n,nc=M)
  tuning       <- expand.grid(Lambda1 ,Lambda2)
  m            <- 1
  X_m_SVD      <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
  V_mat[,m]    <- X_m_SVD$v*X_m_SVD$d
  U_mat[,m]    <- X_m_SVD$u 
  angle        <- function(x,y){
    # compute the arc-cosine between x and y 
    dot.prod   <- x%*%y 
    norm.x     <- norm(x)
    norm.y     <- norm(y)
    theta      <- acos(dot.prod /(norm.x*norm.y))*180/pi
    as.numeric(theta)
  } 
  
  # Initialize
  for(m in 2:M){
    X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
    if(angle(t(X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))<=angle(t((-1)*X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))){
      s=2
    }else{s=1}
    V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d*(-1)^s
    U_mat[,m]  <- X_m_SVD$u*(-1)^(s)
  }
  
  V_wide_mat   <- matrix(rep(V_mat,length(tuning[,1])),nr=p)   #####repeat for lambda
  U_wide_mat   <- matrix(rep(U_mat,length(tuning[,1])),nr=n)
  
  
  pos_c        <- 1:length(tuning[,1])
  
  
  while (length(pos_c)>0){ 
    # update V_wide_mat
    V_wide_mat0                <- V_wide_mat
    pos_long_c                 <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    V_wide_mat[,pos_long_c]    <- Supdate_V_wide(V_wide_mat[,pos_long_c,drop=FALSE],
                                                 U_wide_mat[,pos_long_c,drop=FALSE],
                                                 X,tuning[,1][pos_c],tuning[,2][pos_c])
    # update U_wide_mat
    for (m in 1:M){
      deno                       <- sqrt(colSums((X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE]))^2))+10^(-6)
      deno_mat                   <- matrix(rep(deno,each=n),nr=n)
      U_wide_mat[,M*(pos_c-1)+m] <- X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE])/deno_mat
      
    }
    
    pos_c                      <-unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-3))!=0)-1)%/%M+1)
  }
  # standardize
  V_wide_mat   <- ifelse(abs(V_wide_mat)>10^(-4),V_wide_mat,0)
  v_norm       <- sqrt(colSums(V_wide_mat^2))+0.0001
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)
  V_Dir        <- V_wide_mat/v_mat_norm
  V_Dir        <- ifelse(abs(V_Dir)>10^(-4),V_Dir,0)
  return(list(V_mat=V_Dir))
}


Supdate_V_wide <- function(V_wide_mat,U_wide_mat,X,Lambda1,Lambda2){
  p            <- nrow(V_wide_mat)
  M            <- ncol(V_wide_mat)/length(Lambda1)
  v_norm       <- sqrt(colSums(V_wide_mat^2))                    
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)                
  pos_r        <- rep(TRUE,p)
  pos_c        <- 1:length(Lambda1)
  n            <- nrow(X)
  for(loop in 1:1){
    V_wide_mat0     <- V_wide_mat
    pos_cindex      <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    gl_M            <- Svgl_mat(V_wide_mat[pos_r,pos_cindex,drop=FALSE],
                        U_wide_mat[,pos_cindex],
                        Lambda1[pos_c],Lambda2[pos_c]) 
    
    V_wide_mat[pos_r,pos_cindex]    <- matrix(gl_M,nr=sum(pos_r))               
    
    pos_c       <- unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-4))!=0)-1)%/%M+1)
    pos_r       <- apply(abs(V_wide_mat-V_wide_mat0),1,max) > 10^(-4) 
    if(sum(pos_c)==0 ) break
    
  }
  return(V_wide_mat)
}


Svgl_mat      <- function(V_wide_mat,U_wide_mat,Lambda1,Lambda2){
  p           <- nrow(V_wide_mat)
  M           <- ncol(V_wide_mat)/length(Lambda1)
  Lambda_v    <- rep(Lambda1,each=M)
  Lambda_v    <- as.matrix(Lambda_v,length(Lambda_v),1)
  Lambda_v2   <- rep(Lambda2,each=M)
  Lambda_v2   <- as.matrix(Lambda_v2,length(Lambda_v2),1)
  Vgl_mat     <- matrix(rep(Lambda_v,each=p),nr=p)
  v           <- list(0) 
  u           <- list(0)
  s           <- matrix(0,nr=p,nc=length(Lambda1)*M)
  for (j in 1:p){
    j_index   <- j+p*(c(1:M)-1)
    X_j       <- X[,j_index] 
    for(m in 1:M){
      X_jm    <- matrix(X_j[,m],nr=n)
      U_m     <- U_wide_mat[,m+M*(c(1:length(Lambda1))-1)]
      s[j,m+M*(c(1:length(Lambda1))-1)] <- (1/n)*t(X_jm)%*%U_m+2*Lambda_v2[m+M*(c(1:length(Lambda2))-1)]*
                                           colSums(matrix(V_wide_mat[j,-(m+M*(c(1:length(Lambda1)-1)))],M-1,length(Lambda1)))
    }
  }
  s_norm      <- matrix(0,nr=p,nc=length(Lambda1)*M)
  for (j in 1:p){
    for(m in 1:M){
      s_norm[j,m+M*(c(1:length(Lambda1))-1)] <- sqrt(diag(t(matrix(s[j,],nr=M))%*%matrix(s[j,],nr=M)))
    }
  }
  lambda_mat1  <- matrix(rep(Lambda_v,each=p),nr=p)
  lambda_mat2  <- matrix(rep(Lambda_v2,each=p),nr=p)
  
  s_norm_s                     <- as.vector(s_norm)
  s_norm_s[which(s_norm_s==0)] <- 1
  s_norm_s                     <- matrix(s_norm_s,nr=p)
  Vgl_mat                      <- (s_norm>lambda_mat1)*s*(s_norm-lambda_mat1)/(1/n+2*(M-1)*lambda_mat2)/s_norm_s
  return(Vgl_mat)
}
#============================          iSPCA_S        ==================== 

SVmat_value_sign     <- function(X,Lambda1,Lambda2,p,M){   ####1
  set.seed(1234)
  n                  <- nrow(X)
  V_mat              <- matrix(0,nr=p,nc=M)
  U_mat              <- matrix(0,nr=n,nc=M)
  tuning             <- expand.grid(Lambda1 ,Lambda2)
  m                  <- 1
  X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
  V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d
  U_mat[,m]  <- X_m_SVD$u 
  angle <- function(x,y){
    # compute the arc-cosine between x and y
    dot.prod <- x%*%y 
    norm.x <- norm(x)
    norm.y <- norm(y)
    theta <- acos(dot.prod /(norm.x*norm.y))*180/pi
    as.numeric(theta)
  }
  # Initialize
  for(m in 2:M){
    X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
    if(angle(t(X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))<=angle(t((-1)*X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))){
      s=2
    }else{s=1}
    V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d*(-1)^s
    U_mat[,m]  <- X_m_SVD$u*(-1)^(s)
  }
  
  V_wide_mat   <- matrix(rep(V_mat,length(tuning[,1])),nr=p)   #####repeat for lambda
  U_wide_mat   <- matrix(rep(U_mat,length(tuning[,1])),nr=n)
  
  pos_c        <- 1:length(tuning[,1])
  while (length(pos_c)>0){
    # update V_wide_mat
    V_wide_mat0                <- V_wide_mat
    pos_long_c                 <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    V_wide_mat[,pos_long_c]    <- Supdate_V_wide_sign(V_wide_mat[,pos_long_c,drop=FALSE],
                                                      U_wide_mat[,pos_long_c,drop=FALSE],
                                                      X,tuning[,1][pos_c],tuning[,2][pos_c])
    # update U_wide_mat
    for (m in 1:M){
      deno                       <- sqrt(colSums((X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE]))^2))+10^(-6)
      deno_mat                   <- matrix(rep(deno,each=n),nr=n)
      U_wide_mat[,M*(pos_c-1)+m] <- X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE])/deno_mat
      
    }
    
    pos_c                        <- unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-4))!=0)-1)%/%M+1)
    
  }
  # standardize
  V_wide_mat   <- ifelse(abs(V_wide_mat)>10^(-5),V_wide_mat,0)
  v_norm       <- sqrt(colSums(V_wide_mat^2))+0.0001
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)
  V_Dir        <- V_wide_mat/v_mat_norm
  V_Dir        <- ifelse(abs(V_Dir)>10^(-4),V_Dir,0)
  return(list(V_mat=V_Dir))
}

Supdate_V_wide_sign  <- function(V_wide_mat,U_wide_mat,X,Lambda1,Lambda2){
  p            <- nrow(V_wide_mat)
  M            <- ncol(V_wide_mat)/length(Lambda1)
  v_norm       <- sqrt(colSums(V_wide_mat^2))                    
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)                
  pos_r        <- rep(TRUE,p)
  pos_c        <- 1:length(Lambda1)
  n            <- nrow(X)
  
  for(loop in 1:1){
    V_wide_mat0     <- V_wide_mat
    pos_cindex      <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    
    gl_M_sign       <- Svgl_mat_sign(V_wide_mat[pos_r,pos_cindex,drop=FALSE],
                                  U_wide_mat[,pos_cindex],
                                  Lambda1[pos_c],Lambda2[pos_c]) 
    
    V_wide_mat[pos_r,pos_cindex]    <- matrix(gl_M_sign,nr=sum(pos_r))               
    
    pos_c          <- unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-4))!=0)-1)%/%M+1)
    pos_r          <- apply(abs(V_wide_mat-V_wide_mat0),1,max) > 10^(-4) 
    if(sum(pos_c)==0 ) break
    
  }
  return(V_wide_mat)
}

Svgl_mat_sign <- function(V_wide_mat,U_wide_mat,Lambda1,Lambda2){
  p           <- nrow(V_wide_mat)
  M           <- ncol(V_wide_mat)/length(Lambda1)
  Lambda_v    <- rep(Lambda1,each=M)
  Lambda_v    <- as.matrix(Lambda_v,length(Lambda_v),1)
  Lambda_v2   <- rep(Lambda2,each=M)
  Lambda_v2   <- as.matrix(Lambda_v2,length(Lambda_v2),1)
  
  v           <- list(0) 
  u           <- list(0)
  tau         <- 0.1
  s           <- matrix(0,nr=p,nc=length(Lambda1)*M)
  for (j in 1:p){
    j_index   <- j+p*(c(1:M)-1)
    X_j       <- X[,j_index] 
    for(m in 1:M){
      X_jm    <- matrix(X_j[,m],nr=n)
      U_m     <- U_wide_mat[,m+M*(c(1:length(Lambda1))-1)]
      matrix_m<- V_wide_mat[j,-(m+M*(c(1:length(Lambda1)-1)))]/
                 sqrt(tau^2+V_wide_mat[j,-(m+M*(c(1:length(Lambda1)-1)))]^2)
      s[j,m+M*(c(1:length(Lambda1))-1)]     <- (1/n)*t(X_jm)%*%U_m+2*Lambda_v2[m+M*(c(1:length(Lambda2))-1)]*
      colSums(matrix(matrix_m,M-1,length(Lambda1)))/sqrt(tau^2+V_wide_mat[j,m+M*(c(1:length(Lambda1))-1)]^2)
    }
  }
  s_norm     <-matrix(0,nr=p,nc=length(Lambda1)*M)
  for (j in 1:p){
    for(m in 1:M){
      s_norm[j,m+M*(c(1:length(Lambda1))-1)] <- sqrt(diag(t(matrix(s[j,],nr=M))%*%matrix(s[j,],nr=M)))
    }
  }
  lambda_mat1  <- matrix(rep(Lambda_v,each=p),nr=p)
  lambda_mat2  <- matrix(rep(Lambda_v2,each=p),nr=p)
  
  s_norm_s                     <- as.vector(s_norm)
  s_norm_s[which(s_norm_s==0)] <- 1
  s_norm_s                     <- matrix(s_norm_s,nr=p)
  Vgl_mat_sign                 <- (s_norm>lambda_mat1)*s*(s_norm-lambda_mat1)/(1/n+2*(M-1)*lambda_mat2/(V_wide_mat^2+tau^2))/s_norm_s
  return(Vgl_mat_sign)
}


#============================          SPCA           =====================
lVmat_value    <- function(X,Lambda,p,M){   ####1
  set.seed(1234)
  n            <- nrow(X)
  V_mat        <- matrix(0,nr=p,nc=M)
  U_mat        <- matrix(0,nr=n,nc=M)
  m            <- 1
  X_m_SVD      <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
  V_mat[,m]    <- X_m_SVD$v*X_m_SVD$d
  U_mat[,m]    <- X_m_SVD$u 
  angle        <- function(x,y){
    
    dot.prod   <- x%*%y 
    norm.x     <- norm(x)
    norm.y     <- norm(y)
    theta      <- acos(dot.prod /(norm.x*norm.y))*180/pi
    as.numeric(theta)
  } 
  if(M>1){
    for(m in 2:M){
      X_m_SVD    <- irlba(as.matrix(X[,((m-1)*p+1):(m*p)]),nv=1)
      if(angle(t(X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))<=angle(t((-1)*X_m_SVD$v*X_m_SVD$d),as.matrix(V_mat[,1]))){
        s=2
      }else{s=1}
      V_mat[,m]  <- X_m_SVD$v*X_m_SVD$d*(-1)^s
      U_mat[,m]  <- X_m_SVD$u*(-1)^(s)
    }
  }
  V_wide_mat   <- matrix(rep(V_mat,length(Lambda)),nr=p)   #####repeat for lambda
  U_wide_mat   <- matrix(rep(U_mat,length(Lambda)),nr=n)
  pos_c        <- 1:length(Lambda)
  for(loop  in 1:200){ 
    V_wide_mat0                <- V_wide_mat
    pos_long_c                 <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
    V_wide_mat[,pos_long_c]    <- lupdate_V_wide(V_wide_mat[,pos_long_c,drop=FALSE],
                                                 U_wide_mat[,pos_long_c,drop=FALSE],
                                                 X,Lambda[pos_c])
    for (m in 1:M){
      deno                       <- sqrt(colSums((X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE]))^2))+10^(-6)
      deno_mat                   <- matrix(rep(deno,each=n),nr=n)
      U_wide_mat[,M*(pos_c-1)+m] <- X[,((m-1)*p+1):(m*p)]%*%(V_wide_mat[,M*(pos_c-1)+m,drop=FALSE])/deno_mat
      
    }
    
    pos_c                        <-unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-3))!=0)-1)%/%M+1)
    if(length(pos_c)==0) break
  }
  
  V_wide_mat   <- ifelse(abs(V_wide_mat)>10^(-4),V_wide_mat,0)
  v_norm       <- sqrt(colSums(V_wide_mat^2))+0.0001
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)
  
  V_Dir        <- V_wide_mat/v_mat_norm
  V_Dir        <- ifelse(abs(V_Dir)>10^(-3),V_Dir,0)
  return(list(V_mat=V_Dir))
}
lupdate_V_wide <- function(V_wide_mat,U_wide_mat,X,Lambda){
  p            <- nrow(V_wide_mat)
  M            <- ncol(V_wide_mat)/length(Lambda)
  v_norm       <- sqrt(colSums(V_wide_mat^2))                    
  v_mat_norm   <- matrix(rep(v_norm,each=p),nr=p)                
  pos_r        <- rep(TRUE,p)
  pos_c        <- 1:length(Lambda)
  n            <- nrow(X)
  
  V_wide_mat0  <- V_wide_mat
  pos_cindex   <- rep(M*(pos_c-1),each=M)+rep(1:M,length(pos_c))
  l_M          <- vl_mat(V_wide_mat[pos_r,pos_cindex,drop=FALSE],
                   matrix(U_wide_mat[,pos_cindex],nc=length(pos_cindex)),
                   Lambda[pos_c]) 
  
  V_wide_mat[pos_r,pos_cindex]    <- matrix(l_M,nr=sum(pos_r))               
  
  pos_c          <- unique((which(colSums(abs(V_wide_mat0-V_wide_mat)>10^(-4))!=0)-1)%/%M+1)
  pos_r          <- apply(abs(V_wide_mat-V_wide_mat0),1,max) > 10^(-4) 
  
  return(V_wide_mat)
}
vl_mat         <- function(V_wide_mat,U_wide_mat,Lambda){
  p            <- nrow(V_wide_mat)
  M            <- ncol(V_wide_mat)/length(Lambda)
  Lambda_v     <- rep(Lambda,each=M)
  Lambda_v     <- as.matrix(Lambda_v,length(Lambda_v),1)
  Vl_mat       <- matrix(rep(Lambda_v,each=p),nr=p)
  v            <- list(0) 
  u            <- list(0)
  s            <- matrix(0,nr=p,nc=length(Lambda)*M)
  lambda_mat   <- matrix(rep(Lambda_v,each=p),nr=p)
  for (j in 1:p){
    j_index    <- j+p*(c(1:M)-1)
    X_j        <- matrix(X[,j_index],nc=M) 
    for(m in 1:M){
      X_jm     <- matrix(X_j[,m],nr=n,nc=1)
      U_m      <- matrix(U_wide_mat[,m+M*(c(1:length(Lambda))-1)],nc=length(m+M*(c(1:length(Lambda))-1)))
      s[j,m+M*(c(1:length(Lambda))-1)] <- (1/n)*t(X_jm)%*%U_m
      v1       <- ((1/n)*t(X_jm)%*%U_m-Lambda>0)*n*((1/n)*t(X_jm)%*%U_m-Lambda)
      v2       <- ((1/n)*t(X_jm)%*%U_m+Lambda<=0)*n*((1/n)*t(X_jm)%*%U_m+Lambda)
      Vl_mat[j,m+M*(c(1:length(Lambda))-1)] <- v1+v2
    }
  }
  return(Vl_mat)
}




