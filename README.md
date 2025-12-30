Singular Value Decomposition (SVD) 
=======================================================

This project demonstrates computing the Singular Value Decomposition (SVD) of a matrix in R. 
SVD is a fundamental technique in linear algebra used for matrix factorization, dimensionality reduction, 
and solving systems of linear equations.

---

Workflow Overview
-----------------

1. Define the Matrix
--------------------
- The matrix A is defined either manually or using sequences.
- Examples:
  - A <- matrix(c(2,5,-4,6,3,0), nrow=2, ncol=3)
  - A <- matrix(1:9, nrow=3, ncol=3)

2. Compute A^T * A
------------------
- ATA <- t(A) %*% A
- This is used to compute eigenvalues and eigenvectors, which are needed for SVD.

3. Eigen Decomposition
---------------------
- ATA_e <- eigen(ATA)
- u_mat <- ATA_e$vectors
- lambda <- ATA_e$values
- Very small eigenvalues (less than 1e-12) are set to zero to avoid numerical instability.

4. Construct Sigma Matrix
------------------------
- sigma_indx <- sqrt(lambda) > 0
- r <- diag(sqrt(lambda)[sigma_indx])
- Adjustments are made if the dimensions of r do not match u_mat.
- sigma() function handles zero singular values by removing corresponding columns.

5. Compute V Matrix
------------------
- v_mat <- mapply(function(sigma, u) A %*% u / sigma, sqrt(lambda)[sigma_indx], asplit(u_mat[, sigma_indx], 2))
- This calculates the left singular vectors (V) from A, eigenvectors (U), and singular values (Sigma).

6. Reconstruct Original Matrix
------------------------------
- svd.matrix <- round(v_mat %*% r %*% t(u_mat))
- Reconstructed matrix approximates the original matrix using the SVD components.

---

Key Functions and Concepts
--------------------------
- t(A): Transpose of matrix A
- eigen(): Computes eigenvalues and eigenvectors
- diag(): Creates diagonal matrix for singular values
- mapply(): Applies a function over multiple arguments, used here to compute left singular vectors
- Sigma (Σ): Diagonal matrix of singular values
- U and V matrices: Orthogonal matrices from SVD
- Reconstruction: A ≈ V * Σ * U^T

---

Notes
-----
- Small eigenvalues (< 1e-12) are considered zero to handle numerical errors.
- The code demonstrates SVD computation manually without using R’s built-in svd() function.
- Rounding is applied for easier comparison and interpretation of reconstructed matrices.

---

Usage
-----
1. Define the matrix A.
2. Run the code to compute ATA, eigenvalues, eigenvectors, singular values, and left singular vectors.
3. Use the reconstruction step to verify that the original matrix is approximately recovered.
