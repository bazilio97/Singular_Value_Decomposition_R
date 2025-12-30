# Singular Value Decomposition for the math exam

A <- c(2,5,-4,6,3,0);
A <- c(3,2,2,3,2,-2);
A <- matrix(A,nrow=2,ncol=3)
A <- matrix(1:9,nrow=3,ncol=3)  

ATA <- t(A)%*%A;
ATA_e <- eigen(ATA)
u_mat <- ATA_e$vectors
lambda <- ATA_e$values;
lambda[lambda < 10^-12] <- 0;

sigma_indx <- sqrt(lambda)>0;
r <- diag(sqrt(lambda)[sigma_indx])
if (ncol(r)>nrow(u_mat)) r <- r[,-ncol(r)];
if (ncol(r)<nrow(u_mat)) r <- cbind(r,rep(0,nrow(r)));

sigma <- function(v)
{
   matrix <- diag(sqrt(v));
   
   k <- length(v);
   if (v[k] == 0)
   {
      matrix <- matrix[,-k];
   }
   
}

v_mat <- mapply( function(sigma,u) A%*%u/sigma, sqrt(lambda)[sigma_indx], asplit(u_mat[,sigma_indx],2) )
svd.matrix <- round(v_mat %*% r %*% t(u_mat));
