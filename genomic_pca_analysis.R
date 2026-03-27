# Genomic PCA Analysis in R
# This script preprocesses high-dimensional genomic data,
# applies PCA, visualizes principal component structure,
# and verifies the connection between PCA, SVD, and eigen-decomposition.

library(data.table)
library(ggplot2)

# Set working directory and load data
genome_data <- fread("genome.txt", header = FALSE)

# Extract metadata
id <- genome_data[[1]]
sex <- genome_data[[2]]
population <- genome_data[[3]]

# Extract genomic nucleobase data (ignore first 3 metadata columns)
genomic_data <- genome_data[, -(1:3)]
genomic_matrix <- as.matrix(genomic_data)

# Convert nucleobase data into a binary matrix:
# 0 = matches the mode nucleobase in that column
# 1 = differs from the mode nucleobase
n_rows <- nrow(genomic_matrix)
n_cols <- ncol(genomic_matrix)

mode_col <- character(n_cols)

for (j in 1:n_cols) {
  column_counts <- table(genomic_matrix[, j])
  mode_col[j] <- names(column_counts)[which.max(column_counts)]
}

mode_matrix <- matrix(mode_col, nrow = n_rows, ncol = n_cols, byrow = TRUE)
X <- 1L * (genomic_matrix != mode_matrix)

# Dimension of principal component vectors:
# since X has 10101 genomic positions, each principal component has dimension 10101

# Run PCA on the binary matrix
pca_result <- prcomp(X, center = TRUE, scale. = FALSE)

# Store PCA scores for plotting
scores <- as.data.frame(pca_result$x)
scores$Population <- population
scores$Sex <- sex

# Plot individuals on the first two principal components, colored by population
plot_pc1_pc2 <- ggplot(scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point() +
  labs(
    title = "Projection of Individuals onto the First Two Principal Components",
    x = "PC1",
    y = "PC2",
    color = "Population"
  )

print(plot_pc1_pc2)

# Plot individuals on the first and third principal components, colored by sex
plot_pc1_pc3 <- ggplot(scores, aes(x = PC1, y = PC3, color = Sex)) +
  geom_point() +
  labs(
    title = "Projection of Individuals onto PC1 and PC3",
    x = "PC1",
    y = "PC3",
    color = "Sex"
  )

print(plot_pc1_pc3)

# Inspect the absolute loadings of the third principal component
pc3_loadings <- pca_result$rotation[, 3]

loading_df <- data.frame(
  index = 1:length(pc3_loadings),
  abs_loading = abs(pc3_loadings)
)

plot_pc3_loadings <- ggplot(loading_df, aes(x = index, y = abs_loading)) +
  geom_line() +
  labs(
    title = "Nucleobase Index vs Absolute Loading of PC3",
    x = "Nucleobase Index",
    y = "|PC3 Loading|"
  )

print(plot_pc3_loadings)

# PCA / SVD verification on a small synthetic dataset

set.seed(123)

N <- 5
D <- 3
toy_matrix <- matrix(rnorm(N * D), nrow = N, ncol = D)

# Center the data
toy_matrix_centered <- scale(toy_matrix, center = TRUE, scale = FALSE)

# Singular Value Decomposition
svd_result <- svd(toy_matrix_centered)
U <- svd_result$u
sigma <- svd_result$d
V <- svd_result$v

# Sample covariance matrix
S <- (t(toy_matrix_centered) %*% toy_matrix_centered) / N

# Eigen-decomposition of covariance matrix
eigen_result <- eigen(S)
Q <- eigen_result$vectors
lambda <- eigen_result$values

# Compare covariance eigenvalues with sigma^2 / N
comparison_table <- cbind(
  eigenvalues_S = lambda,
  sigma_sq_over_N = sigma^2 / N
)

print(toy_matrix)
print(toy_matrix_centered)
print(svd_result)
print(S)
print(lambda)
print(Q)
print(comparison_table)

# Compare eigenvectors from covariance matrix with right singular vectors
print(abs(cor(V, Q)))
