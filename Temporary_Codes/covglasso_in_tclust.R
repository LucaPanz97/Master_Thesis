### TESI

library(mclust)
library(tclust)
library(ellipse)
library(MASS)
library(robustbase)
library(covglasso)


##	calculates the initial cluster assignment and parameter
.InitClusters <- function (X, iter, pa)
{
  iter$assig <- sample(x = pa$K,size = pa$n,replace = T)
  iter$assig[sample(1:pa$n, pa$trimm)] <- 0
  iter$csize <- tabulate(iter$assig, nbins = pa$K)
  iter$cw <- iter$csize / pa$no.trim
  assignment <- iter$assig[iter$assig != 0]
  
  for (l in 1:pa$K) {
    assign(paste0("matrice_", l), matrix(nrow=0, ncol=ncol(X)))
  }
  for (i in 1:pa$no.trim) {
    index <- assignment[i]
    matrixx <- get(paste0("matrice_", index))
    matrixx <- rbind(matrixx, X[i,])
    assign(paste0("matrice_", index), matrixx)
  }
  
  lista_matrix <- list()
  for (l in 1:pa$K) {
    matricee <- get(paste0("matrice_", l))
    lista_matrix[[l]] <- matricee
  }
  
  for (k in seq_along(lista_matrix)) {
    n_covlas <- dim(lista_matrix[[k]])[1]  
    p_covlas <- dim(lista_matrix[[k]])[2]
    iter$center[k,] <- colMeans(lista_matrix[[k]])
    iter$S[,,k] <- cov(lista_matrix[[k]])*(n_covlas-1)/n_covlas
    if (p_covlas >= n_covlas)
      iter$S[,,k] <- iter$S[,,k] + 1e-6 * diag(p_covlas)
    covfit <- covglasso(data = lista_matrix[[k]], S = iter$S[,,k], n = n_covlas, lambda = pa$lambda, rho = NULL,
                        duplicated = TRUE, L = 10, crit = c("bic", "ebic"), gamma = 0.5, penalize.diag = FALSE,
                        start = NULL, ctrl = control(), path = FALSE)
    iter$sigma[,,k] <- covfit$sigma
  }
  return (iter)
}


##  Multivariate normal density
.dmnorm <- function (X, mu, sigma)
{
  ((2 * pi)^(-length(mu) / 2)) * (det(sigma)^(-1/2)) * exp (-0.5 * .evmaha (X, mu, sigma))
}


.evmaha <- function (X, mu, sigma)    ##  calculate mahalanobis distances 
{                                     
  v <- eigen (sigma)
  Zt <- t (v$vectors) %*% (t (X) - mu)
  colSums ((Zt * v$values^(-.5))^2)
}


##	get a matrix object out of the sigma - tensor
.ssclmat <- function (x, k) as.matrix (x[,,k])  


##	calculates the obj. functions value
.calcobj <- function (X, iter, pa)
{	
  iter$obj <- 0
  for (k in 1:pa$K)
  {
    if (det(iter$sigma[,,k]) != 0){
      log_det_sigma <- log(det(iter$sigma[,,k]))
      trace_term <- sum(diag(solve(iter$sigma[,,k]) %*% iter$S[,,k]))
      sparsity_term <- pa$lambda * sum(abs(pa$P * iter$sigma[,,k]))
      
      iter$obj <- iter$obj + log_det_sigma + trace_term + sparsity_term
    }
    else{
      iter$obj <- Inf
      iter$esc_init <- TRUE
      break
    }
  }
  return (iter)
}


##	finds the current cluster assignment based on the given center and sigma
.findClustAssig <- function (X, iter, pa)
{														
  ll = matrix (NA, pa$n, pa$K)
  for (k in 1:pa$K)			
    ll[,k] <- iter$cw[k] * .dmnorm(X, iter$center[k,], .ssclmat (iter$sigma,k))   ##	calculating the likelihood for each observation in cluster k --> D_k(X;theta) = p_k*phi(X;m_k,S_k) --> inner part of (2.1)
  
  old.assig <- iter$assig
  iter$assig <- apply(ll,1,which.max)						
  disc <- ll [cbind (1:pa$n, iter$assig)]   #gli elementi di disc sono D(x_1;theta),D(x_2;theta),...,D(x_n;theta) dove D(x_i;theta)=max{D_1(x_i;theta),...,D_K(x_i;theta)} (vedi (3.1))  
  idx.out <- which (rank (disc, ties.method = c ("random"))<=pa$trim)	
  iter$assig [idx.out] <- 0								
  iter$code <- all (old.assig == iter$assig)				
  iter$csize <- tabulate (iter$assig, pa$K)
  iter$z_ij <- matrix (0, ncol = pa$K, nrow = pa$n)		
  iter$z_ij[cbind (1:pa$n, iter$assig)[-idx.out, ]] <- 1	
  
  if (!pa$equal.weights)									
    iter$cw <- iter$csize / sum(iter$csize)
  return (iter)
}


.estimClustPar <- function (X, iter, pa)   #vedi step 2.2. dell'algoritmo
{	
  for (l in 1:pa$K) {
    assign(paste0("matrice_", l), matrix(nrow=0, ncol=ncol(X)))
  }
  for (i in 1:pa$no.trim) {
    index <- (mclust::map(iter$z_ij))[i]
    matrixx <- get(paste0("matrice_", index))
    matrixx <- rbind(matrixx, X[i,])
    assign(paste0("matrice_", index), matrixx)
  }
  
  lista_matrix <- list()
  for (l in 1:pa$K) {
    matricee <- get(paste0("matrice_", l))
    lista_matrix[[l]] <- matricee
  }
  
  for (k in seq_along(lista_matrix)) {
    if (length(lista_matrix[[k]]) != 0){
      if (dim(lista_matrix[[k]])[1] > 2){
        n_covlas <- dim(lista_matrix[[k]])[1]  
        p_covlas <- dim(lista_matrix[[k]])[2]
        iter$center[k,] <- colMeans(lista_matrix[[k]])
        iter$S[,,k] <- cov(lista_matrix[[k]])*(n_covlas-1)/n_covlas
        if (p_covlas >= n_covlas)
          iter$S[,,k] <- iter$S[,,k] + 1e-6 * diag(p_covlas)
        covfit <- covglasso(data = lista_matrix[[k]], S = iter$S[,,k], n = n_covlas, lambda = pa$lambda, rho = NULL,
                            duplicated = TRUE, L = 10, crit = c("bic", "ebic"), gamma = 0.5, penalize.diag = FALSE,
                            start = NULL, ctrl = control(), path = FALSE)
        iter$sigma[,,k] <- covfit$sigma
      }
      else{ 
        iter$sigma[,,k] <- 0
        iter$esc_init <- TRUE
      }
    }
    else{
      iter$sigma[,,k] <- 0
      iter$esc_init <- TRUE
    }
  }
  return (iter)
}


##	converts the output of function .tclust.R to a "tclust" - object
.Parsetclust.Res <- function (x, iter, parlist)
{	
  idx.clust <- order (iter$csize, decreasing = TRUE)   #Crea un vettore contenente gli indici dei cluster ordinati per cluster size decrescente
  idx.nz <- iter$csize[idx.clust] != 0   #TRUE se cluster size != 0 e FALSE altrimenti
  idx.clust <- idx.clust [idx.nz]   #Tiene gli indici di idx.clust per cui la cluster size è != 0
  id.clust <- 0
  id.clust [1 + idx.clust] <- 1:length (idx.clust)   #ipotizziamo che il cluster n.3 sia quello con cluster size piu grande, e quindi idx.clust[1]=3; allora id.clust[3+1]=id.clust[4]=1
  #ipotizziamo che il cluster n.1 sia quello con quarta cluster size piu grande, e quindi idx.clust[4]=1; allora id.clust[1+1]=id.clust[2]=4
  #ipotizziamo che il cluster n.2 sia quello con terza cluster size piu grande, e quindi idx.clust[3]=2; allora id.clust[2+1]=id.clust[3]=3
  #inoltre id.clust[1]=0 sempre
  ret <- list (
    centers = t (iter$center[idx.clust, , drop = FALSE]),
    cov = iter$sigma [,, idx.clust, drop = FALSE],
    cluster = id.clust [iter$assig + 1],   #ipotizziamo che la prima osservazione stia nel cluster n.2, quindi iter$assig[1]=2, quindi iter$assig[1]+1=2+1=3, quindi id.clust[iter$assig[1]+1]=id.clust[3]=3, quindi cluster[1]=id.clust[iter$assig[1]+1]=3, ovvero la prima osservazione sta nel terzo cluster più grande per cluster size
    par = parlist,
    k = length (idx.clust),
    obj = iter$obj,
    loglik_vec = iter$loglik_vec,
    num_init = iter$num_init,
    size = iter$csize [idx.clust],
    weights = iter$cw [idx.clust],
    ret.orig = iter
  )
  class (ret) <- "tclust"
  ret
}


## funzione per printare i risultati
print.tclust <- function (x, ...)
{
  cat ("* Results for TCLUST algorithm: *\n")
  cat ("trim = ", x$par$alpha, ", k = ", x$k, "\n", sep = "")
  
  cat ("Classification (trimmed points are indicated by 0", "):\n")
  print (x$cluster)
  
  cat ("Means:\n")
  print (x$centers)
  
  cat ("Covariances:\n")
  print (x$cov)
  
  if (x$obj < (-1e+20))
    warning ("The solution is not reliable. More iterations are probably needed.")
  
  cat ("\nInitialization with best results: ", x$num_init, "\n")
  
  cat ("\nTrimmed objective function: ", x$obj, "\n")
  
  cat ("\nLoglikelihood vector: ", x$loglik_vec, "\n")
  
  invisible(x)
}


.covglasso_in_tclust <- function (x, k, alpha, nstart = 50, iter.max = 20, equal.weights = FALSE, zero.tol = 1e-16, lambda, P)
{
  if (is.data.frame(x))
    x <- data.matrix(x)
  else if (!is.matrix (x))
    x <- matrix(x, nrow = length(x), ncol = 1, dimnames = list (names(x), deparse(substitute(x))))
  if (!is.numeric (x))
    stop ("parameter x: numeric matrix/vector expected")
  
  parlist <- list (k = k, alpha = alpha, nstart = nstart, iter.max = iter.max, equal.weights = equal.weights, zero.tol = zero.tol, lambda = lambda, P = P)
  
  n <- nrow (x)
  p <- ncol (x)
  no.trim <- floor(n * (1 - alpha))
  
  # preparing lot's of lists	
  pa <- list (					       ##	these are variables which all the iterations have in common, and will never change (pa for "params")
    n = n,													 ##	number of observations
    p = p,													 ##	number of dimensions
    no.trim = no.trim,							 ##	number of observations which are considered as to be not outlying
    trimm = n-no.trim,							 ##	number of observations which are considered as to be outlying
    alpha = alpha,                   ## trimming level
    K = k,													 ##	number of clusters to be searched for
    equal.weights = equal.weights,   ##	wether equal weights shall be assumed for all clusters
    zero.tol = zero.tol,					   ##	zero tolerance
    lambda = lambda,
    P = P
  )
  
  iter <- list (							## these variables change each iteration - this object is passed to all the functions, modified and returned by them
    obj = Inf,												       ##	current objective value
    loglik_vec = c(), 
    assig = array (0, n),							       ##	cluster assignment
    csize = array (NA, k),						       ##	cluster sizes
    cw = rep (NA, k),									       ##	cluster weights
    sigma = array (NA, c (p, p, k)),         ##	cluster's sigmas
    center = array (NA, c(k, p)),			       ##	cluster's centers
    S = array (NA, c (p, p, k)),             ## cluster's sample covariance matrices
    code = FALSE,												     ##	this is a return code supplied by functions like .findClustAssig
    esc_init = FALSE,                        ## if TRUE, it ends that initialization
    num_init = 0,
    z_ij = matrix (0, nrow = n, ncol = k )   ##	z_ij ## -> what was it's old name - "ind" - right?
  )
  
  best.iter <- iter   # initializing the best.iter - structure
  
  for (j in 1:nstart)
  {
    cat ("\nINIZIALIZZ NUM: ", j, "\n")
    iter$loglik_vec <- c()
    iter$esc_init <- FALSE
    iter$num_init <- j
    iter <- .InitClusters (x, iter, pa)
    
    for (i in 0:iter.max)
    {
      temp <- i
      cat ("\nITERAZ NUM: ", i, "\n")
      cat ("\nPRE ASSIG: ", iter$assig, "\n")
      iter <- .findClustAssig (x, iter, pa)
      cat ("\nPOS ASSIG: ", iter$assig, "\n")
      
      if (iter$code || i == iter.max)									
        break										
      
      iter <- .calcobj(x, iter, pa)
      if (iter$esc_init)
        break
      
      iter$loglik_vec[temp+1] <- iter$obj
      
      iter <- .estimClustPar (x, iter, pa)				## estimates the cluster's parameters (cov, center) based on the current iter$assig
      if (iter$esc_init)
        break
    }
    
    cat ("\nESC INIT: ", iter$esc_init, "\n")
    
    if (!iter$esc_init){
      iter <- .calcobj(x, iter, pa)
      iter$loglik_vec[temp+1] <- iter$obj
      cat ("\nFINAL OBJ: ", iter$obj, "\n")
      iter$code = as.numeric (i == iter.max)   #iter$code=1 se l'algoritmo ha raggiunto iter.max, cioè non è arrivato a convergenza prima, mentre	iter$code=0 se l'algoritmo è arrivato a convergenza prima di raggiungere iter.max 	
      if (j == 1 || iter$obj < best.iter$obj)
        best.iter = iter
    }
  }
  return (.Parsetclust.Res (x, best.iter, parlist))
}


#GENERAZIONE DATI

#set.seed(547)     # to replicate example

p <- 50

# mixture parameters
mu1 <- c(1:50)
mu2 <- c(3,4,11:28,21:50)
mu3 <- c(5,6,21:38,21:50)
sigma1 <- diag(seq(0.1, 1, length.out = 50))
sigma2 <- diag(seq(0.1, 1, length.out = 50))
sigma3 <- diag(seq(0.1, 1, length.out = 50))
tau <- c(0.4, 0.3, 0.3)

mu_out <- c(9,10,81:128)
sigma_out <- diag(seq(0.5, 1.5, length.out = 50))

# 3-component mixture distribution
K <- 3
N_no_out <- 95
Nk <- rowSums( rmultinom(N_no_out, 1, tau) )
group <- factor(c(rep(1:K, Nk),rep("out",5)))   # group membership vector

# generate data
x1 <- mvrnorm(Nk[1], mu1, sigma1)
x2 <- mvrnorm(Nk[2], mu2, sigma2)
x3 <- mvrnorm(Nk[3], mu3, sigma3)
x_out <- mvrnorm(5, mu_out, sigma_out)
x <- rbind(x1, x2, x3, x_out)
N <- nrow(x)
plot(x[,1:2])
plot(x[,1:2], col = group)


#RICHIAMO IL .covglasso_in_tclust()

tclust <- .covglasso_in_tclust(x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, equal.weights = FALSE, zero.tol = 1e-16, lambda = 8, P = matrix(1, nrow = p, ncol = p))

plot(x[,1:2], col=tclust$cluster+1,pch=tclust$cluster+1, asp=1)

for(k in 1:3) {
  lines(ellipse(x = tclust$cov[,,k], centre = tclust$centers[,k]), lty = 3, col=k+1)
}

#Con 10 colonne diverse (delta uguale a 2) e 40 uguali la clusterizzazione trovata è mooolto brutta (però sembra esserci comunque un trend crescente del loglik vector)
#Uguale con 20 colonne diverse e 30 uguali (delta uguale a 2)
#Con 20 colonne diverse e 30 uguali (delta uguale a 10) clusterizzazione, locations e scatters sono molto buoni e il loglik vector ha un trend crescente




