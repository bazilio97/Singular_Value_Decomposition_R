############################################################
# Copyright (c) June 2023.
# Dmytro Kotov, mail to dmytro.kotov97@gmail.com
# ================================================================
# README for SVD-based Calculation in R
# ================================================================
#
# This script computes a singular value decomposition-like approximation
# using eigen decomposition on the matrix ATA, followed by adjustment
# based on small eigenvalues and the construction of the 'sigma' matrix
# for the final approximation.
#
# The following steps outline the key operations and their purpose:
#
# 1. **ATA Calculation:**
#    - Compute ATA where A is the input matrix.
#    - ATA is a square matrix that encodes important properties of A.
#
# 2. **Eigen Decomposition of ATA:**
#    - Perform eigen decomposition to obtain eigenvalues and eigenvectors.
#    - `ATA_e$values` contains eigenvalues (λ), and `ATA_e$vectors` contains eigenvectors (U).
#
# 3. **Thresholding of Eigenvalues:**
#    - Eigenvalues less than 10^-12 are considered as numerical noise and are set to zero.
#
# 4. **Sigma Matrix Construction:**
#    - Only positive square roots of eigenvalues are used to form a diagonal matrix `r` (the "sigma" matrix).
#    - A conditional check ensures proper alignment of the matrix dimensions by adjusting its columns if necessary.
#
# 5. **Sigma Function Definition:**
#    - A helper function `sigma(v)` creates a diagonal matrix using the square root of non-zero eigenvalues.
#    - If the last eigenvalue is zero, the last column of the matrix is removed.
#
# 6. **Calculation of Approximation:**
#    - `v_mat` computes the matrix approximation by dividing the product of matrix A and eigenvectors by the corresponding singular values.
#    - Finally, `svd.matrix` is obtained by multiplying `v_mat` with the matrix `r` and the transpose of the eigenvector matrix.
#
# 7. **Output:**
#    - The result `svd.matrix` is the matrix that approximates the singular value decomposition of the original matrix A.
#
# ================================================================
# Detailed Code Walkthrough:
# ================================================================
#
# # Step 1: Compute ATA (A Transpose A)  
ATA <- t(A) %*% A;

# # Step 2: Eigen decomposition of ATA  
ATA_e <- eigen(ATA);

# Extract eigenvalues (λ) and eigenvectors (U)
u_mat <- ATA_e$vectors;
lambda <- ATA_e$values;

# # Step 3: Threshold small eigenvalues  
# Set eigenvalues smaller than 10^-12 to zero (numerical noise threshold)
lambda[lambda < 10^-12] <- 0;

# # Step 4: Sigma matrix construction  
# Determine which eigenvalues correspond to non-zero singular values
sigma_indx <- sqrt(lambda) > 0;

# Construct the "sigma" matrix (diagonal matrix of non-zero singular values)
r <- diag(sqrt(lambda)[sigma_indx]);

# Ensure correct matrix dimensions by adjusting rows and columns
if (ncol(r) > nrow(u_mat)) r <- r[,-ncol(r)];
if (ncol(r) < nrow(u_mat)) r <- cbind(r, rep(0, nrow(r)));

# # Step 5: Define sigma function  
# This function constructs a diagonal matrix for singular values
sigma <- function(v)
{
   matrix <- diag(sqrt(v));
   
   k <- length(v);
   if (v[k] == 0)
   {
      matrix <- matrix[,-k];  # Remove the last column if eigenvalue is zero
   }
   
}

# # Step 6: Calculate matrix approximation  
# Use mapply to apply the sigma function across the relevant eigenvalues and eigenvectors
v_mat <- mapply(function(sigma, u) A %*% u / sigma, sqrt(lambda)[sigma_indx], asplit(u_mat[,sigma_indx], 2));

# # Final approximation of the matrix (SVD-like result)
svd.matrix <- round(v_mat %*% r %*% t(u_mat));

# Output the approximation matrix  
# svd.matrix is the result of the operation, which approximates A
#
# ================================================================
# How to Use:
# ================================================================
#
# 1. Replace the matrix `A` with your data matrix before running the code.
# 2. The final result will be stored in the variable `svd.matrix`, which approximates the singular value decomposition of A.
# 3. You can visualize or analyze `svd.matrix` depending on your needs.
#
# ================================================================
# Notes:
# ================================================================
#
# - The precision of the approximation depends on the choice of the threshold (10^-12).
# - This script provides a way to approximate SVD when dealing with large or ill-conditioned matrices.
#
# ================================================================
# END OF DOCUMENTATION
# ================================================================
