library(RSpectra)
library(biganalytics)
library(kmed)
library(aricode)

# Build graph Laplacian
Build_Gaussian_kernel = function(data, kNN = 10){
  n = nrow(data)
  kNN = min(kNN, n)
  
  # make dist sparse (only keep kNN dist)
  knn_list <- FNN::get.knn(data, k = kNN, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  knn_dist = knn_list$nn.dist
  
  # Gaussian kernel (adaptive bandwidth)
  ep = 1
  epi = ep * (apply(knn_dist,1,median)^2)
  knn_dist = knn_dist^2 / epi
  knn_dist = exp(-knn_dist)
  
  W = Matrix::sparseMatrix(i = rep(1:nrow(Idx),each = kNN), 
                           j = c(t(Idx)),
                           dims = c(n, n),
                           x = c(t(knn_dist)))
  W = as.matrix(W)
  
  # Symmetrize 
  diag(W) = 1
  W = (W + t(W)) / 2
  
  return(W)
  
}
# Calculate EigenVal & EigenVec
ncut = function(W,nbEigenValues,ifnormalized = TRUE,methodgraph = "beltrami"){
  # input: W = symmetric similarity matrix
  #        nbEigenValues =  number of Ncut eigenvectors computed
  #        ifnormalized: # TRUE for eigenfunctions of D^{-1} W, otherwise orthonormal basis of D^{-1/2}WD^{-1/2}
  # output: Eigenvectors= continuouse Ncut eigenvectos, dim = nrow(W) x nbEigenValues
  #         Eigenvalues= Ncut eigenvalues, dim = 1 x nbEigenValues
  # methodgraph: beltrami
  
  n = dim(W)[1]
  nbEigenValues = min(nbEigenValues,n)
  
  # from symmetric kernel to symmetric laplacian (beltrami method)
  if(methodgraph == "beltrami"){
    valeurMin = 1e-8
    W =W * (W>valeurMin)
    d = as.matrix(rowSums(W) + .Machine$double.eps)
    P = W / (d %*% t(d) )
    d2 = as.matrix(sqrt(rowSums(P)) + .Machine$double.eps)
    P = P / (d2 %*% t(d2) )
    
    P = (P + t(P)) / 2 # force symmetry
  }
  
  
  # obtain eigen vector & eigen values
  if(nbEigenValues < dim(P)[1]){
    res = RSpectra::eigs_sym(P, k = nbEigenValues, which = "LM")
  }else{
    res = eigen(P)
  }
  Eigenvalues = res$values
  Eigenvectors = res$vectors
  
  # reorder eigvals & eigvecs
  sub = order(Eigenvalues,decreasing = TRUE)
  Eigenvalues = Eigenvalues[sub]
  Eigenvectors = Eigenvectors[,sub]
  
  # remove negative eigvals & corresponding eigvecs
  sub = which(Eigenvalues > 0)
  Eigenvalues = Eigenvalues[sub]
  Eigenvectors = Eigenvectors[,sub]
  
  d = Eigenvectors[,1]; d[d == 0] = 1
  if(ifnormalized == TRUE){
    for(j in 1:ncol(Eigenvectors)){
      Eigenvectors[,j] = Eigenvectors[,j] / d
      if(Eigenvectors[1,j] != 0){
        Eigenvectors[,j] = Eigenvectors[,j] * sign(Eigenvectors[1,j])
      }
    }
  }
  
  # index = order(abs(diff(Vals)),decreasing = TRUE)[1]
  # Eigenvectors = Eigenvectors[,1:index]
  # Eigenvalues = Eigenvalues[1:index]
  return(list(Eigenvectors,Eigenvalues))
}
# Spectral Clustering
Spectral_clustering = function(sim_mat, k, nbEigenValues = 50){
  res = ncut(W = sim_mat,nbEigenValues = nbEigenValues)
  Vecs = res[[1]]; Vals = res[[2]]
  tmp_Vecs = Vecs[,1:k,drop = FALSE]
  
  # QR rotation (https://github.com/asdamle/QR-spectral-clustering/blob/master/clusterQR.m)
  piv = qr(t(tmp_Vecs),LAPACK = TRUE)
  piv = piv$pivot[1:ncol(tmp_Vecs)]
  
  res <- svd(t(tmp_Vecs[piv,,drop = FALSE]))
  tmp_Vecs = tmp_Vecs %*% (res$u %*% t(res$v))
  
  res_cluster = apply(abs(tmp_Vecs),1,function(x){which(x == max(x))})
  return(res_cluster)
}

# Diffusion Distance
Calculate_diffusiondist = function(sim_mat, nbEigenValues = 100, time_scale, delta = 0.1){
  res = ncut(W = sim_mat,nbEigenValues = nbEigenValues)
  Vecs = res[[1]]; Vals = res[[2]]
  
  Vecs = Vecs[,-1,drop = FALSE]
  Vals = Vals[-1]
  
  vec_sub = min(which(Vals ^ time_scale < delta),length(Vals))
  Vecs = Vecs[,1:vec_sub,drop = FALSE]; Vals = Vals[1:vec_sub]
  
  cat("diffusion time scaler : ", time_scale, "\n", 
      "# of Eigen for graph : ", vec_sub, "\n")
  
  coord = t(diag(x = Vals ^ time_scale, nrow = length(Vals)) %*% t(Vecs))
  diffusion_dist = stats::dist(coord, method = "euclidean")
  
  return(diffusion_dist)
}

# Hierachical Clustering
Hierarchical_clustering <- function(diffusion_dist,height = NULL,
                                    k = NULL){
  hc = hclust(diffusion_dist, method = "complete")
  if(!is.null(height)){
    cluster = cutree(hc, h = height)
  }else if(!is.null(k)){
    cluster = cutree(hc, k = k)
  }else{
    cluster = NULL
  }
  res = list(cluster, hc)
  names(res) = c("partition","htree")
  return(res)
}

kmedoid_clustering <- function(diffusion_dist,k = NULL){
  res = kmed::skm(diffusion_dist, ncluster = k, seeding = 10, iterate = 100)
  return(res$cluster)
}

# Twin-sample cross validation (using Rand Index as cross validation measure)
Twin_sample_metric = function(diffuse_dist, seed){
  diffuse_dist = as.matrix(diffuse_dist)
  n = nrow(diffuse_dist)
  if(n/2 < 10){return(1)}
  
  # split sample into two twin samples
  set.seed(seed)
  sample1 = sample(n)[1:floor(n/2)]
  sample2 = setdiff(1:n,sample1)
  
  sample1_nn_idex = apply(diffuse_dist[sample1,sample2],2,function(x){which(x == min(x))})
  
  # Clustering (a series of k)
  # k_ls = c(seq(floor(length(sample1)/10),1,-2),1)
  k_ls = seq(10,1,-1)
  RI_ls = unlist(lapply(k_ls,function(k){
    sample1_cl = kmedoid_clustering(
      as.dist(diffuse_dist[sample1,sample1]),
      k = k)
    sample2_cl = kmedoid_clustering(
      as.dist(diffuse_dist[sample2,sample2]),
      k = k)
    sample1_nn_cl = sample2_cl[sample1_nn_idex]
    RI(sample1_cl,sample1_nn_cl)
  }))
  
  # # Clustering (a series of height)
  # height = Hierarchical_clustering(as.dist(diffuse_dist[sample1,sample1]))[[2]]
  # height_ls = height$height[height$height > median(height$height)]
  # 
  # RI_ls = unlist(lapply(height_ls,function(h){
  #   sample1_cl = Hierarchical_clustering(
  #     as.dist(diffuse_dist[sample1,sample1]),
  #     height = h)[[1]]
  #   sample2_cl = Hierarchical_clustering(
  #     as.dist(diffuse_dist[sample2,sample2]),
  #     height = h)[[1]]
  #   sample1_nn_cl = sample2_cl[sample1_nn_idex]
  #   RI(sample1_cl,sample1_nn_cl)
  # }))

  return(median(RI_ls))
  
}
