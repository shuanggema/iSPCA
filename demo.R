#============================================================================
# program: demo.R
# note   : code submission for paper "Integrative Sparse Principal Component 
#                                    Analysis of Gene Expression Data".
# purpose: illustrate how to implement the proposed method 
# steps:
#      (1) source R files to load pre-defined functions
#      (2) generate datasets
#      (3) specify tunings
#      (4) implement the algorithms
#      (5) choose tuning

#============================================================================

rm(list=ls()) 

# (1) source R files to load pre-defined functions
library(MASS)
library(Matrix)
library(irlba)
source('function.R')

# (2) generate datasets
loop    <- 1
set.seed(loop)
M       <- 6        # the number of datasets
n       <- 25       # the number of observations
p       <- 500      # the dimension

alpha_v          <- rep(0.5,M)
beta_v           <- rep(0.3,M)
Type             <- c('ss1') # ss1 denotes Scenario 1; ss2 denotes Scenario 2


V               <- matrix(rep(diag(0,p,p),M),nr=p)
for(m in 1:M){
  beta          <- beta_v[m]
  alpha         <- alpha_v[m]
  l             <- 100
  set.seed(m*l+loop*200)
  
  if(Type=='ss1'){ ##### V for Scenario 1
    V[,(m-1)*p+1]            <- c(rep(1,floor(p^beta)),rep(0,p-floor(p^beta)))
    V[,((m-1)*p+2):(m*p)]    <- matrix(runif((p-1)*p,1,2),nr=p,nc=p-1)
    V[,((m-1)*p+1):(m*p)]    <- qr.Q(qr(V[,((m-1)*p+1):(m*p)] ))
  }
  
  
  if(Type=='ss2'){ ##### V for Scenario 2
    V_sig                    <- rnorm(floor(p^beta),mean=(floor(p^beta):1)^1.5,sd=floor(p^beta):1/10)
    V[,(m-1)*p+1]            <- c(sample(V_sig,floor(p^beta)),rep(0,p-floor(p^beta)))
    V[,((m-1)*p+2):(m*p)]    <- matrix(runif((p-1)*p,1,2),nr=p,nc=p-1)
    V[,((m-1)*p+1):(m*p)]    <- qr.Q(qr(V[,((m-1)*p+1):(m*p)]))
  }
  
  if(Type=='dsp'){ ##### V for Scenario 3
    
    V[,(m-1)*p+1]            <- rep(0,p)
    if(beta<0.4){
      pos_1                  <- (1+(m-1)*floor(1*(p^beta)/6)):(m*floor(1*p^beta/6))
      V[pos_1,(m-1)*p+1]     <- rnorm(floor(1*p^beta/6), mean=3, sd=0.2) 
      pos_2                  <-  (M*floor(1/6*p^beta)+1):(M*floor(1/6*p^beta)+floor(p^beta)-floor(1/6*p^beta))
      V[pos_2,(m-1)*p+1]     <- runif(floor(p^beta)-floor(1*p^beta/6),0.5,1)  
    }else{
      pos_1                  <- (1+(m-1)*floor(1*(p^beta)/10)):(m*floor(1*p^beta/10))
      pos_2                  <-  (M*floor(1/10*p^beta)+1):(M*floor(1/10*p^beta)+floor(p^beta)-floor(1/10*p^beta))
      V[pos_1,(m-1)*p+1]     <- rnorm(floor(1*p^beta/10), mean=3, sd=0.2) 
      V[pos_2,(m-1)*p+1]     <- runif(floor(p^beta)-floor(1*p^beta/10),0.5,1)  
    }
    
    
    
    
    V[,((m-1)*p+2):(m*p)]    <- matrix(runif((p-1)*p,1,2),nr=p,nc=p-1)
    V[,((m-1)*p+1):(m*p)]    <- qr.Q(qr(V[,((m-1)*p+1):(m*p)]))
  }
  
  
  V_org            <- V[,((m-1)*p+1):(m*p)]  
  eig_v            <- c(p^alpha,rep(1,p-1))
  assign(paste0('Sigma_',m), V_org%*%diag(eig_v)%*%t(V_org))		
}


set.seed(loop*300)

X                <- NULL
X_c              <- NULL
for (m in 1:M){
  assign(paste0('X',m), mvrnorm(n, rep(0,p), Sigma=get(paste0('Sigma_',m))))
  X              <- cbind(X,get(paste0('X',m)))
  X_c            <- rbind(X_c,get(paste0('X',m)))
}


# (3) specify tunings 
lLambda1        <- seq(0,0.5,length=100)

# (4) implement the algorithm for a single merged data set
n               <- nrow(X_c)
comb_lSSVD_X    <- lVmat_value(X_c,lLambda1,p,M=1)

# (4) implement the algorithm for SPCA
n               <- nrow(X)
lSSVD_X         <- lVmat_value(X,lLambda1,p,M)


# (3) specify tunings
Lambda1         <- seq(0,2,length=15)      # for iSPCA
Lambda21        <- seq(0,2,length=15)      # for iSPCA_M
Lambda22        <- seq(0,0.05,length=15)   # for iSPCA_M
Lambda31        <- seq(0,2,length=15)      # for iSPCA_S
Lambda32        <- seq(0,0.04,length=20)   # for iSPCA_S
pLambda1        <- rep(0,1)                # for PCA


# (4) implement the algorithm for iSPCA
SSVD_X              <- Vmat_value(X,Lambda1,p,M)

# (4) implement the algorithm for iSPCA_M
SSVD_X_prom         <- SVmat_value(X,Lambda21,Lambda22,p,M)

# (4) implement the algorithm for iSPC_S
SSVD_X_prom_sign    <- SVmat_value_sign(X,Lambda31,Lambda32,p,M)

# (4) implement the algorithm for PCA
pSSVD_X             <- Vmat_value(X,pLambda1,p,M)









#=========================       choose tuning        ======================= 



set.seed(378378)
X_train            <- NULL
X_train_c          <- NULL

for (m in 1:M){
  assign(paste0('X',m), mvrnorm(2*n, rep(0,p), Sigma=get(paste0('Sigma_',m))))
  X_train          <- cbind(X_train,get(paste0('X',m)))
  X_train_c        <- rbind(X_train_c,get(paste0('X',m)))
}


###### for iSPCA; iSPCA_M; iSPCA_S; PCA
tuning             <- expand.grid(Lambda21 ,Lambda22)
tuning2            <- expand.grid(Lambda31 ,Lambda32)

CV_error           <- vector() # for iSPCA
CV_error_prom      <- vector() # for iSCPA_M
CV_error_prom_sign <- vector() # for iSCPA_S

for(i in 1:length(Lambda1)){
  CV_error[i]    <-  0
  for(m in 1:M){
    CV_error[i]  <-  CV_error[i]+(norm(X_train[,((m-1)*p+1):(m*p)]-
                                         X_train[,((m-1)*p+1):(m*p)]%*%SSVD_X$V_mat[,M*(i-1)+m]%*%t( SSVD_X$V_mat[,M*(i-1)+m]),type='F'))
    
  }
}

for(i in 1:length(tuning[,1])){
  
  CV_error_prom[i]    <-  0
  for(m in 1:M){
    CV_error_prom[i]  <-  CV_error_prom[i]+(norm(X_train[,((m-1)*p+1):(m*p)]-
                                                 X_train[,((m-1)*p+1):(m*p)]%*%SSVD_X_prom$V_mat[,M*(i-1)+m]%*%t(SSVD_X_prom$V_mat[,M*(i-1)+m]),type='F'))
    
  }
}

for(i in 1:length(tuning2[,1])){
  
  CV_error_prom_sign[i]    <-  0
  for(m in 1:M){
    CV_error_prom_sign[i]  <-  CV_error_prom_sign[i]+(norm(X_train[,((m-1)*p+1):(m*p)]-
                                                           X_train[,((m-1)*p+1):(m*p)]%*%SSVD_X_prom_sign$V_mat[,M*(i-1)+m]%*%t(SSVD_X_prom_sign$V_mat[,M*(i-1)+m]),type='F'))
    
  }
}
i_th            <- which.min(CV_error)
i_th_prom       <- which.min(CV_error_prom)
i_th_prom_sign  <- which.min(CV_error_prom_sign)
Lambda          <- Lambda1[i_th ]
lambda_1        <- tuning[,1][i_th_prom ]
lambda_2        <- tuning[,2][i_th_prom ]
lambda_1_sign   <- tuning2[,1][i_th_prom_sign ]
lambda_2_sign   <- tuning2[,2][i_th_prom_sign ]

v_1_mat         <- SSVD_X$V_mat[,M*(i_th-1)+(1:M)]           # for iSPCA
v_1_mat_prom    <- SSVD_X_prom$V_mat[,M*(i_th_prom-1)+(1:M)] # for iSPCA_M
v_1_mat_prom_sign   <- SSVD_X_prom_sign$V_mat[,M*(i_th_prom_sign-1)+(1:M)] # for iSPCA_S

pv_1_mat        <- pSSVD_X$V_mat[,1:M]                       # PCA

###### for SPCA
CV_error         <- vector()
for(i in 1:length(lLambda1)){
  CV_error[i]    <-  0
  for(m in 1:M){
    CV_error[i]  <-  CV_error[i]+(norm(X_train[,((m-1)*p+1):(m*p)]-
                                         X_train[,((m-1)*p+1):(m*p)]%*%lSSVD_X$V_mat[,M*(i-1)+m]%*%t(lSSVD_X$V_mat[,M*(i-1)+m]),type='F'))
    
  }
}
i_th             <- which.min(CV_error)
Lambda_lasso     <- lLambda1[i_th]
v_l_mat          <- lSSVD_X$V_mat[,M*(i_th-1)+(1:M)]       # for SPCA


#######  for a single merged data set
CV_error         <- vector()
for(i in 1:length(lLambda1)){
  CV_error[i]    <-  0
  
  CV_error[i]    <-  CV_error[i]+(norm(X_train_c-
                                       X_train_c%*%comb_lSSVD_X$V_mat[,i]%*%t(comb_lSSVD_X$V_mat[,i]),type='F'))
  
}

i_th            <- which.min(CV_error)
Lambda_lasso    <-lLambda1[i_th]
v_comb_l_mat    <- matrix(comb_lSSVD_X$V_mat[,i_th],nc=1) # for a single merged data set

#=========================       End  file      ============================= 