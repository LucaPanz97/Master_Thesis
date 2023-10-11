### TESI

library(mclust)
library(tclust)
library(ellipse)
library(MASS)
library(robustbase)


#FUNZIONI UTILIZZATE IN .tclust.R()

##	calculates the initial cluster sizes
.getini <- function (K, no.trim)
{
  if (K == 1)
    return (no.trim)
  
  pi.ini  <- runif(K)   #Genera K numeri casuali compresi tra 0 e 1
  ni.ini <- sample(x = K, size = no.trim, replace = T, prob = pi.ini / sum (pi.ini))   #Seleziona no.trim numeri casuali (con replacement) dalla sequenza di numeri interi da 1 a K, i cui pesi sono determinati dalla distribuzione di probabilità 'pi.ini/sum(pi.ini)'      
  return (tabulate(ni.ini, nbins = K))   #Ritorna un vettore di lunghezza K che che conta il numero di volte che ciascun valore compreso tra 1 e K è stato selezionato dalla funzione 'sample'
}


##	calculates the initial cluster assignment and parameter
.InitClusters <- function (X, iter, pa)
{
  for (k in 1:pa$K)
  {
    idx <- sample (1:pa$n, pa$p+1)   #Genera un campione casuale di p+1 numeri interi compresi tra 1 e n
    X.ini = X [drop = F,idx,]   #Estrae le righe della matrice X corrispondenti ai numeri interi presenti in idx, mantenendo tutte le colonne
    iter$center [k,] <- colMeans (X.ini)   #Calcola la media di ogni colonna della matrice X.ini e assegna i risultati alla k-esima riga della matrice iter$center = m_k_0		
    iter$sigma[,,k] <- (pa$p/(pa$p+1))*cov (X.ini)   #Calcola la matrice di covarianza del subset di dati X.ini, moltiplicandola per il fattore correttivo 'p/(p+1)' e assegna i risultati alla k-esima "fetta" della matrice iter$sigma = S_k_0
  }
  
  if (pa$equal.weights)   #	if we're considering equal weights, cw is set here **ONLY ONCE AND NEVER CHANGED**
    iter$csize <- rep (pa$no.trim / pa$K, pa$K)
  else
    iter$csize = .getini (pa$K, pa$no.trim)
  iter$cw <- iter$csize / pa$no.trim   #p_0
  # if we're considering different weights, calculate them, and they're gonna be recalculated every time in .estimClustPar
  return (iter)
}


## restricts the clusters' covariance matrices by averaging them
.restr.avgcov <- function (iter, pa)
{
  s.all <- matrix (0, pa$p, pa$p)
  for (k in 1:pa$K)
    s.all <- s.all + iter$sigma[,,k] * iter$csize[k]/sum(iter$csize)
  iter$sigma[,,] <- s.all   #non sta assegnando un array bidimensionale a un  array tridimensionale???	
  iter$code = sum(diag (s.all)) > pa$zero.tol
  return (iter)
}


.restr2_eigenv <- function (autovalues, ni.ini, restr.fact, zero.tol)
{
  ev <- autovalues
  if (!is.matrix (ev))												
    if (is.atomic (ev))												
      ev <- t (as.numeric (ev))									
  else															
    ev <- as.matrix (ev)										
  stopifnot (ncol (ev) == length (ni.ini))   # check wether the matrix autovalues and ni.ini (i.e. vector with current sample size of the clusters) have the right dimension  
  d <- t (ev)   #matrice che ha per elementi d_jl, con j=1,...,K e l=1,...p
  p <- nrow (ev)
  K <- ncol (ev)	
  n <- sum(ni.ini)
  nis <- matrix(data=ni.ini,nrow=K,ncol=p)   
  idx.nis.gr.0 <- nis > zero.tol								
  used.ev <- ni.ini > zero.tol   #indici delle colonne, cioè dei cluster, che hanno current sample size > zero.tol										
  ev.nz <- ev[,used.ev]   # non-zero eigenvalues
  if ((max (ev.nz) <= zero.tol))										
    return (matrix (0, nrow = p, ncol = K))							
  if (max (ev.nz) / min (ev.nz) <= restr.fact)					
  {																
    ev[,!used.ev] <- mean (ev.nz)									
    return (ev)														
  }																
  
  d_ <- sort (c (ev, ev / restr.fact))   #d_ contiene i valori e_1 <= e_2 <= ... <= e_2kp (vedi Proposition 3.2.)								
  dim <- length (d_)   #2kp											
  d_1 <- d_
  d_1[dim+1] <- d_[dim] * 2
  d_2 <- c (0, d_)
  ed <- (d_1 + d_2) / 2   #ed contiene i valori f_1,f_2,...,f_2kp+1 t.c. f_1 <= e_1 <= f_2 <= e_2 <= ... <= f_2kp <= e_2kp <= f_2kp+1
  dim <- dim + 1;   #2kp+1
  
  t <- s <- r <- array(0, c(K, dim))
  sol <- sal <- array(0, c(dim))
  for (mp_ in 1:dim)
  {
    for (i in 1:K)
    {
      r[i,mp_] <- sum ((d[i,] < ed[mp_])) + sum((d[i,] > ed[mp_]*restr.fact))   #inner part a denominatore di m_i (vedi Proposition 3.2.)
      s[i,mp_] <- sum (d[i,]*(d[i,] < ed[mp_]))   #primo pezzo dell'inner part a numeratore di m_i
      t[i,mp_] <- sum (d[i,]*(d[i,] > ed[mp_] * restr.fact))   #secondo pezzo dell'inner part a numeratore di m_i
    }
    sol[mp_] <- sum (ni.ini / n * (s[,mp_] + t[,mp_] / restr.fact)) / (sum(ni.ini / n * (r[, mp_])))   #m_i
    e <-	sol[mp_] * (d < sol[mp_]) + d * (d >= sol[mp_]) * (d <= restr.fact * sol[mp_]) + (restr.fact*sol[mp_]) * (d > restr.fact * sol[mp_])   #d_jl_m_i (vedi (3.2))
    o <- -1/2 * nis / n * (log(e) + d / e)   #objective function in (3.3) senza la sommatoria iniziale  #perchè moltiplicato per (-1/2) e diviso per n (come sopra nel calcolo di sol) ???
    sal[mp_] <- sum(o)   ##objective function in (3.3)
  }
  m <- sol[which.max (sal)]   #m_opt  #nel paper voglio minimizzare l'objective function in (3.3), ma qui questa viene moltiplicata per (-1/2) e quindi voglio massimizzarla					
  t (m * (d < m) + d * (d >= m) * (d <= restr.fact * m) + (restr.fact * m) * (d > restr.fact * m))   #matrice che ha per elementi d_lj_m_opt, con l=1,...p e j=1,...,K
}


.multbyrow <- function (a, b) t (t(a) * as.numeric (b))


.HandleSmallEv <- function (autovalues, zero.tol)							
{	
  K <- nrow (autovalues)   #K non è il numero di colonne???													
  autovalues[autovalues <= zero.tol] <- zero.tol							
  mi <- apply(autovalues,2,min)											
  ma <- apply(autovalues,2,max)											
  idx.iter <- which (mi/ma <= zero.tol)   #indice dei cluster per cui min_autoval / max_autoval <= zero.tol									
  for (i in idx.iter)														
    autovalues[autovalues[,i] > mi[i] / zero.tol, i] <- mi[i] / zero.tol
  det = apply(autovalues, 2, prod)											
  p <- nrow (autovalues)													
  autovalues_det <- .multbyrow (autovalues, det^(-1/p))						
  return (autovalues_det)
}


.restr2_deter_ <- function (autovalues, ni.ini, restr.fact, zero.tol)
{
  p = nrow (autovalues)
  if (p == 1)
    return  (.restr2_eigenv (autovalues, ni.ini, restr.fact, zero.tol))
  K = ncol (autovalues)
  es = apply (autovalues, 2, prod)
  idx.ni.ini.gr.0 <- ni.ini > zero.tol										
  if (max(es[idx.ni.ini.gr.0]) <= zero.tol)	
    return (matrix (0, p, K))												
  d = t(es)	
  autovalues_det <- .HandleSmallEv (autovalues, zero.tol)					
  if (max (d[idx.ni.ini.gr.0]) / min (d[idx.ni.ini.gr.0]) <= restr.fact)		
  {
    d [!idx.ni.ini.gr.0] <-  mean (d[idx.ni.ini.gr.0])						
    dfin <- d^(1/p)
  }
  else
    dfin <- .restr2_eigenv (d^(1/p), ni.ini, restr.fact^(1/p), zero.tol)   #???
  .multbyrow (autovalues_det, dfin)											
}


##	restricts the clusters' covariance matrices by restricting it's eigenvalues or determinants without restricts directions 
.restr.diffax <- function (iter, pa, f.restr.eigen = .restr2_eigenv, restr.deter, restr.fact = 12)   #restr.fact è c nel paper
{
  if (pa$p == 1)						#	one - dimensional
    if (restr.fact == 1)				#	all variances should be equal -> use the simpler .restr.avgcov function instead
      return (.restr.avgcov (iter, pa))	
  else									
    restr.deter <- FALSE				#	else - if p == 1 always use the simpler eigen - restriction 
  
  if (!missing (restr.deter))					#	evaluating the appropriate restriction function
    f.restr.eigen <- if (restr.deter) .restr2_deter_ else .restr2_eigenv							
  
  u <- array (NA, c(pa$p, pa$p, pa$K))
  d <- array (NA, c(pa$p, pa$K))
  for (k in 1:pa$K)
  {
    ev <- eigen (iter$sigma[,,k])
    u [,,k] <- ev$vectors
    d [,k] <- ev$values
  }
  d [d < 0] <- 0		##	all eigenvalue < 0 are restricted to 0 BECAUSE: min (d) < 0 -> max (d) / min (d) < 0 which would always fit our eigenvalue - criterion!
  d <- f.restr.eigen (d, iter$csize, restr.fact, pa$zero.tol)
  iter$code = max(d) > pa$zero.tol
  if (!iter$code)
    return (iter)
  for (k in 1:pa$K)	
    iter$sigma[,,k] <- u[,,k] %*% diag (d[,k], nrow = pa$p) %*% t(u[,,k])   #S_j_l+1 = (U_j)'(D_j*)(U_j)  con  D_j* = diag(d_lj_m_opt)
  return (iter)
}


##  Multivariate normal density
.dmnorm <- function (X, mu, sigma)
{
  ((2 * pi)^(-length(mu) / 2)) * (det(sigma)^(-1/ 2)) * exp (-0.5 * .evmaha (X, mu, sigma))
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
.calcobj <- function (X, iter, pa)   #???
{			
  iter$obj <- ifelse (pa$equal.weights, 0, sum ( (iter$csize * log(iter$cw)  )[iter$csize > pa$zero.tol] ) )   #iter$obj = 0 se pa$equal.weights = TRUE, altrimenti iter$obj è uguale alla somma dei valori (per cui iter$csize > pa$zero.tol) del vettore prodotto tra iter$csize e log(iter$cw)
  for (k in 1:pa$K)
  {
    w <- .dmnorm(X,iter$center[k,],.ssclmat (iter$sigma,k))
    if (sum(iter$z_ij[  ,k]) >  pa$zero.tol)
      if (sum(iter$z_ij[w <= 0, k]) <= pa$zero.tol)
        iter$obj <- iter$obj + sum(iter$z_ij[w > 0,k] * log(w[w > 0]))
    else
    {
      iter$obj <- iter$obj -Inf   # we can already return -Inf, no need to do any further calculation
      return (iter)					
    }
  }
  return (iter)
}


.TreatSingularity <- function (iter, pa)
{	
  warning ("After trimming, all points in the data set are concentrated in k subspaces.")   # a single point is a subspace too
  iter$code <- 2	 # indicating the data's concentration in either k subspaces or points
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
  for (k in 1:pa$K)
  {
    if (iter$csize[k] > pa$zero.tol)   ##	this cluster's size is > 0
    {
      iter$center[k,] = (t(iter$z_ij[,k]) %*% X) / iter$csize[k]   #m_k_l+1   #iter$z_ij[,k] è un vettore di lunghezza n che indica l'appartenenza dei punti al cluster k (cioè se il punto i appartiene al cluster k, iter$z_ij[i,k] sarà uguale a 1, altrimenti sarà uguale a 0   
      X.c <- (X - matrix (iter$center[k,], ncol = pa$p, nrow = pa$n, byrow = TRUE))
      iter$sigma[,,k] <- (t(X.c * iter$z_ij[,k]) %*% X.c) / iter$csize[k]   #T_k
    }
    else   ##	this cluster's size has decreased to 0
      iter$sigma[,,k] <- 0					
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
  int <- list (                                      #ipotizziamo che il cluster n.1 sia quello con quarta cluster size piu grande, e quindi idx.clust[4]=1; allora id.clust[1+1]=id.clust[2]=4
    iter.successful = 0,                             #ipotizziamo che il cluster n.2 sia quello con terza cluster size piu grande, e quindi idx.clust[3]=2; allora id.clust[2+1]=id.clust[3]=3
    iter.converged = 0,                              #inoltre id.clust[1]=0 sempre
    dim = dim (x)
  )
  ret <- list (
    centers = t (iter$center[idx.clust, , drop = FALSE]),
    cov = iter$sigma [,, idx.clust, drop = FALSE],
    cluster = id.clust [iter$assig + 1],   #ipotizziamo che la prima osservazione stia nel cluster n.2, quindi iter$assig[1]=2, quindi iter$assig[1]+1=2+1=3, quindi id.clust[iter$assig[1]+1]=id.clust[3]=3, quindi cluster[1]=id.clust[iter$assig[1]+1]=3, ovvero la prima osservazione sta nel terzo cluster più grande per cluster size
    par = parlist,
    k = length (idx.clust),
    obj = iter$obj,
    loglik_vec = iter$loglik_vec,
    size = iter$csize [idx.clust],
    weights = iter$cw [idx.clust],
    ret.orig = iter,
    int = int
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
  
  cat ("\nTrimmed objective function: ", x$obj, "\n")
  
  cat ("\nLoglikelihood vector: ", x$loglik_vec, "\n")
  
  if (!is.null (x$restr.fact))
    cat ("Selected restriction factor:", x$restr.fact, "\n")
  
  invisible(x)
}


.tclust.R <- function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, f.restr = .restr.diffax, equal.weights = FALSE, zero.tol = 1e-16)
{
  if (is.data.frame(x))
    x <- data.matrix(x)
  else if (!is.matrix (x))
    x <- matrix(x, nrow = length(x), ncol = 1, dimnames = list (names(x), deparse(substitute(x))))
  if (!is.numeric (x))
    stop ("parameter x: numeric matrix/vector expected")
  
  parlist <- list (k = k, alpha = alpha, nstart = nstart, iter.max = iter.max, f.restr = f.restr, equal.weights = equal.weights, zero.tol = zero.tol)
  
  n <- nrow (x)
  p <- ncol (x)
  no.trim <- floor(n * (1 - alpha))
  
  # preparing lot's of lists	
  pa <- list (					       ##	these are variables which all the iterations have in common, and will never change (pa for "params")
    n = n,													 ##	number of observations
    p = p,													 ##	number of dimensions
    no.trim = no.trim,							 ##	number of observations which are considered as to be not outlying
    trimm = n-no.trim,							 ##	number of observations which are considered as to be outlying
    K = k,													 ##	number of clusters to be searched for
    equal.weights = equal.weights,   ##	wether equal weights shall be assumed for all clusters
    zero.tol = zero.tol						   ##	zero tolerance				
  )
  
  iter <- list (							## these variables change each iteration - this object is passed to all the functions, modified and returned by them
    obj = -Inf,												       ##	current objective value
    loglik_vec = c(), 
    assig = array (0, n),							       ##	cluster assignment
    csize = array (NA, k),						       ##	cluster sizes
    cw = rep (NA, k),									       ##	cluster weights
    sigma = array (NA, c (p, p, k)),         ##	cluster's sigmas
    center = array (NA, c(k, p)),			       ##	cluster's centers
    code = NA,												       ##	this is a return code supplied by functions like .findClustAssig
    z_ij = matrix (0, nrow = n, ncol = k )   ##	z_ij ## -> what was it's old name - "ind" - right?
  )
  
  best.iter <- iter   # initializing the best.iter - structure
  
  for (j in 1:nstart)
  {
    iter <- .InitClusters (x, iter, pa)
    lastobj <- -Inf
    iter$loglik_vec = c()
    for (i in 0:iter.max)
    {
      temp <- i
      iter <- f.restr (iter = iter, pa = pa)   # restricting the clusters' scatter structure
      if (!iter$code)
      {								
        if (i)
          return (.Parsetclust.Res (x, .TreatSingularity (.calcobj (x, iter, pa), pa), parlist))
        else
          iter$sigma[,,] = diag (pa$p)
      }
      iter <- .findClustAssig (x, iter, pa)				
      if (iter$code || i == iter.max)									
        break										
      iter <- .calcobj(x, iter, pa)
      iter$loglik_vec[temp+1] <- iter$obj
      iter <- .estimClustPar (x, iter, pa)				## estimates the cluster's parameters (cov, center) based on the current iter$assig
    }
    iter <- .calcobj(x, iter, pa)
    iter$loglik_vec[temp+1] <- iter$obj
    iter$code = as.numeric (i == iter.max)   #iter$code=1 se l'algoritmo ha raggiunto iter.max, cioè non è arrivato a convergenza prima, mentre	iter$code=0 se l'algoritmo è arrivato a convergenza prima di raggiungere iter.max 	
    if (j == 1 || iter$obj > best.iter$obj)
      best.iter = iter
  }
  return (.Parsetclust.Res (x, best.iter, parlist))
}


#GENERAZIONE DATI

#set.seed(547)     # to replicate example

p <- 2

# mixture parameters
mu1 <- c(-5,6)
mu2 <- c(3,2)
mu3 <- c(-3,-4)

sigma1 <- diag(0.4, 2)
sigma2 <- diag(2, 2)
sigma3 <- diag(1, 2)
tau <- c(0.4, 0.3, 0.3)

# 3-component mixture distribution
K <- 3
N_no_out <- 95
Nk <- rowSums( rmultinom(N_no_out, 1, tau) )
group <- factor(c(rep(1:K, Nk),rep("out",5)))   # group membership vector

# generate data
x1 <- mvrnorm(Nk[1], mu1, sigma1)
x2 <- mvrnorm(Nk[2], mu2, sigma2)
x3 <- mvrnorm(Nk[3], mu3, sigma3)
x_out <- rbind(cbind(runif(n = 5,min = 8,max = 10),runif(n = 5,min = -6,max = -3))) # add some outliers
x <- rbind(x1, x2, x3, x_out)
N <- nrow(x)
plot(x)
plot(x, col = group)


#RICHIAMO IL .tclust.R()

tclust <- .tclust.R(x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, f.restr = .restr.diffax, equal.weights = FALSE, zero.tol = 1e-16)

plot(x, col=tclust$cluster+1,pch=tclust$cluster+1, asp=1)

for(k in 1:3) {
  lines(ellipse(x = tclust$cov[,,k], centre = tclust$centers[,k]), lty = 3, col=k+1)
}



